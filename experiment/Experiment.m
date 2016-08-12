%% EXPERIMENT
% *Summary of this class goes here*
%
% Detailed explanation goes here
%
%% Copyright
% * *2016 - TNO*
% * *Author*: Coen van Leeuwen
% * *Since*: August 12, 2016
%
%% See also:
%

%% Class definition:
classdef Experiment < handle
    
    %% Public properties
    properties
        % The type of solver to use for initialization
        initSolverType@char matrix;
        
        % The type of solver to use for iteration
        iterSolverType@char matrix;
    end
    
    %% Final protected properties (settings)
    properties (SetAccess = immutable, GetAccess = protected)
        % The type of constraints to use between variables
        constraintType@char matrix;
        
        % Arguments to initialize the constraints with
        constraintArgs@cell;
        
        % Properties to provide to the agents
        agentProps@struct scalar;
        
        % Number of colors (domain size)
        nColors@uint16 scalar;
        
        % The number of iterations that the algorithm is allowed to improve
        % on its last best result. A value of zero indicates no limit
        nStableIterations@uint16 scalar;
        
        % The total maximum number of iterations (overrules
        % nStableIterations). A value of zero indicates no limit
        nMaxIterations@uint16 scalar;
        
        % The maximum running time in seconds for non-iterative solvers
        maxtime@double scalar;
        
        % The delay between checks to see if the algorithm has converged
        % for non-iterative solvers
        waittime@double scalar;
    end
    
    %% Publicly obtainable properties
    properties (SetAccess = protected, GetAccess = public)
        % A cell array containing all variables of the graph
        variable@cell;
        
        % A cell array containing all constraints of the graph
        constraint@cell;
        
        % A cell array containing all variable agents the graph
        agent@cell;
        
        % A cell array containing all constraint agents of the graph
        constraintAgent@cell;
        
        % A struct containing properties of the graph (constraints,
        % density, size)
        graph@struct scalar;
        
        % A struct contaning the results of the experiment
        results@struct scalar;
    end
    
    %% Protected properties
    properties (SetAccess = protected, GetAccess = protected)
        % Only for internal use, the time at which the experiment starts
        t_start@uint64 scalar; 
    end
    
    %% Public methods
    methods
        
        %% EXPERIMENT - Object constructor
        function obj = Experiment(edges, options)
            % Parse all options
            obj.initSolverType = getSubOption('nl.coenvl.sam.solvers.CoCoASolver', 'char', options, 'initSolverType');
            obj.iterSolverType = getSubOption('', 'char', options, 'iterSolverType');
            obj.constraintType = getSubOption('nl.coenvl.sam.constraints.InequalityConstraint', 'char', options, 'constraint', 'type');
            obj.constraintArgs = getSubOption({}, 'cell', options, 'constraint', 'arguments');
            obj.agentProps = getSubOption(struct, 'struct', options, 'agentProperties');
            obj.nColors = getSubOption(uint16(3), 'uint16', options, 'ncolors');
            obj.nStableIterations = getSubOption(uint16(50), 'uint16', options, 'nStableIterations');
            obj.nMaxIterations = getSubOption(uint16(0), 'uint16', options, 'nMaxIterations');
            obj.maxtime = getSubOption(180, 'double', options, 'maxTime');
            obj.waittime = getSubOption(.01, 'double', options, 'waitTime');
            
            assert(~isempty(obj.initSolverType) || ~isempty(obj.iterSolverType), ...
                'EXPERIMENT:INIT:NOSOLVERTYPE', ...
                'Either set initSolverType, iterSolverType or both');
            assert(isGraph(edges), 'EXPERIMENT:INITGRAPH:INVALIDGRAPH', ...
                'Input must be a Nx2 matrix of connected nodes');
            assert(graphIsConnected(edges), ...
                'EXPERIMENT:INITGRAPH:NOTCONNECTED', ...
                'Graph must be connected');
            
            % Make sure all referenced java objects are cleared
            obj.reset();
            
            obj.graph.edges = edges;
            obj.graph.density = graphDensity(edges);
            obj.graph.size = graphSize(edges);
            
            obj.initVariables();
            obj.assignAgentProperties();
            obj.initConstraints();
            
            if ~isempty(strfind(obj.iterSolverType, 'MaxSum'))
                obj.initConstraintAgents();
            end
        end % EXPERIMENT
        
        %% SETINITSOLVERTYPE - Sets the solver type and resets experiment
        function set.initSolverType(obj, type)
            obj.initSolverType = type;
            obj.reset();
        end % SETINITSOLVERTYPE
        
        %% SETITERSOLVERTYPE - Sets the solver type and resets experiment
        function set.iterSolverType(obj, type)
            obj.iterSolverType = type;
            obj.reset();
        end % SETITERSOLVERTYPE
        
        %% RUN - Run the experiment
        function run(obj)
            if ~isempty(obj.initSolverType)
                obj.runInitSolver();
            end
            
            if ~isempty(obj.iterSolverType)
                obj.runIterSolver();
            end
            
            % Gather results to return
            obj.results.time = toc(obj.t_start);
            obj.results.cost = obj.getCost();
            obj.results.evals = nl.coenvl.sam.ExperimentControl.getNumberEvals();
            obj.results.msgs = nl.coenvl.sam.MailMan.getTotalSentMessages();
        end % RUN
        
        %% RESET - Resets the experiment so it can be run again
        function reset(obj)
            obj.results = struct();
            cellfun(@(x) x.reset(), obj.agent);
            cellfun(@(x) x.clear(), obj.variable);
            nl.coenvl.sam.ExperimentControl.ResetExperiment();
        end % RESET
        
        %% DELETE - Object constructor
        function delete(obj)
            % Experiment destructor
            obj.reset();
        end % DELETE
        
    end
    
    methods (Sealed, Access = protected)
        %% PAUSEUNTILVARIABLESARESET
        function pauseUntilVariablesAreSet(obj)
            while toc(obj.t_start) < obj.maxtime
                if all(cellfun(@(x) x.isSet(), obj.variable));
                    return
                else
                    pause(obj.waittime);
                end
            end
            
            % At this moment, we exceeded maximum time, stop experiment
            error('Solver did not terminate in %d seconds, increase maximum wait time if it is expected to terminate.', ...
                obj.maxtime);
        end % PAUSEUNTILVARIABLESARESET
        
        %% DOSTOP - Decide if the experiment should stop
        function bool = doStop(obj, numIters, countDown)
            bool = false;
            % First decide based on MAX number of iterations
            if (obj.nMaxIterations > 0) && (numIters >= obj.nMaxIterations)
                bool = true;
            end
            
            % Only then look at number of stable iterations
            if ~bool && (obj.nStableIterations > 0) && (countDown <= 0)
                bool = true;
            end
        end % DOSTOP
        
        %% GETCOST - Obtain current experiment cost
        function cost = getCost(obj)
            cost = sum(cellfun(@(x) x.getExternalCost(), obj.constraint));
        end % GETCOST
        
        %% GETSOLVERCOUNTERPART - Get the constraint (MAXSUM) solver type 
        function type = getSolverCounterPart(obj)
            dummyVariable = nl.coenvl.sam.variables.IntegerVariable(int32(1), int32(1));
            dummyAgent = nl.coenvl.sam.agents.VariableAgent(dummyVariable, 'dummy');
            dummySolver = feval(obj.iterSolverType, dummyAgent);
            
            assert(isa(dummySolver, 'nl.coenvl.sam.solvers.MaxSumVariableSolver'), ...
                'EXPERIMENT:initConstraintAgents:INVALIDSOLVERTYPE', ...
                'Unexpected solver type, constraint agents only apply to MaxSum');
            
            type = dummySolver.getCounterPart().getCanonicalName();
        end % GETSOLVERCOUNTERPART
    end
    
    %% Protected methods
    methods (Access = protected)

        %% INITVARIABLES - Initialize variables and agents
        function initVariables(obj)
            nagents = obj.graph.size;
            for i = 1:nagents
                varName = sprintf('variable%05d', i);
                agentName = sprintf('agent%05d', i);
                
                obj.variable{i} = nl.coenvl.sam.variables.IntegerVariable(int32(1), int32(obj.nColors), varName);
                obj.agent{i} = nl.coenvl.sam.agents.VariableAgent(obj.variable{i}, agentName);
                
                if ~isempty(obj.initSolverType)
                    obj.agent{i}.setInitSolver(feval(obj.initSolverType, obj.agent{i}));
                end
                if ~isempty(obj.iterSolverType)
                    obj.agent{i}.setIterativeSolver(feval(obj.iterSolverType, obj.agent{i}));
                end
            end
        end % INITVARIABLES
        
        %% ASSIGNAGENTPROPERTIES
        function assignAgentProperties(obj)
            fields = fieldnames(obj.agentProps);
            nagents = numel(obj.agent);
            
            for i = 1:nagents
                for f = fields'
                    prop = f{:};
                    if numel(agentPropts) == nagents
                        obj.agent{i}.set(prop, obj.agentProps(i).(prop));
                    elseif numel(agentPropts) == 1 && numel(obj.agentProps.(prop)) == 1
                        obj.agent{i}.set(prop, obj.agentProps.(prop));
                    elseif numel(agentPropts) == 1 && numel(obj.agentProps.(prop)) >= nagents
                        obj.agent{i}.set(prop, obj.agentProps.(prop)(i));
                    else
                        error('DOEXPERIMENT:INCORRECTPROPERTYCOUNT', ...
                            'Incorrect number of properties, must be either 1 or number of agents (%d)', ...
                            nagents);
                    end
                end
            end
        end % ASSIGNAGENTPROPERTIES
        
        %% INITCONSTRAINTS - Add constraints to the graph
        function initConstraints(obj)
            for i = 1:size(obj.graph.edges,1)
                a = obj.graph.edges(i,1);
                b = obj.graph.edges(i,2);
                
                if numel(obj.constraintArgs) <= 1
                    % Create a constraint always with same argument
                    obj.constraint{i} = feval(obj.constraintType, obj.variable{a}, obj.variable{b}, obj.constraintArgs{:});
                elseif numel(obj.constraintArgs) == size(edges,1)
                    % Create a constraint for every argument
                    obj.constraint{i} = feval(obj.constraintType, obj.variable{a}, obj.variable{b}, obj.constraintArgs{i});
                elseif mod(numel(obj.constraintArgs), size(edges,1)) == 0
                    error('Deprecated style of constraint argumentations!');
                else
                    error('DOEXPERIMENT:INCORRECTARGUMENTCOUNT', ...
                        'Incorrect number of constraint arguments, must be 0, 1 or number of edges (%d)', ...
                        size(edges,1));
                end
                
                obj.agent{a}.addConstraint(obj.constraint{i});
                obj.agent{b}.addConstraint(obj.constraint{i});
            end
            
        end % INITCONSTRAINTS
        
        %% INITCONSTRAINTAGENTS - Add constraint agents (MAXSUM only)
        function initConstraintAgents(obj)
            functionSolverType = obj.getSolverCounterPart();
            
            for i = 1:size(obj.graph.edges,1)
                a = obj.graph.edges(i,1);
                b = obj.graph.edges(i,2);
                
                % Create constraint agent
                agentName = sprintf('constraint%05d', i);
                obj.constraintAgent{i} = nl.coenvl.sam.agents.ConstraintAgent(agentName, obj.constraint{i}, obj.variable{a}, obj.variable{b});
                obj.constraintAgent{i}.setSolver(feval(functionSolverType, obj.constraintAgent{i}));
                
                % Set constraint agent address as targets
                obj.agent{a}.addFunctionAddress(obj.constraintAgent{i}.getID());
                obj.agent{b}.addFunctionAddress(obj.constraintAgent{i}.getID());
            end
        end % INITCONSTRAINTAGENT
        
        %% RUNINITSOLVER - Run initialization part of algorithm
        function runInitSolver(obj)
            a = obj.agent{randi(obj.graph.size)};
            a.set(nl.coenvl.sam.solvers.CoCoSolver.ROOTNAME_PROPERTY, true);
            
            obj.t_start = tic; % start the clock
            
            cellfun(@(x) x.init(), obj.agent);
            cellfun(@(x) x.init(), obj.constraintAgent);
            obj.pauseUntilVariablesAreSet();
        end % RUNINITSOLVER
        
        %% RUNITERSOLVER - Run iterative part of algorithm
        function runIterSolver(obj)
            % Do the iterations
            obj.results.costList = [];
            obj.results.evalList = nl.coenvl.sam.ExperimentControl.getNumberEvals();
            obj.results.msgList = nl.coenvl.sam.MailMan.getTotalSentMessages();
            obj.results.timeList = toc(obj.t_start);
            
            % Start iterations
            i = 0;
            bestSolution = inf;
            countDown = obj.nStableIterations;
            fprintf('Iteration: ');
            while ~obj.doStop(i, countDown)
                countDown = countDown - 1;
                i = i + 1;
                if mod(i, 25) == 0
                    fprintf(' %d', i);
                end
                
                cellfun(@(x) x.tick(), obj.agent);
                cellfun(@(x) x.tick(), obj.constraintAgent);
                
                cost = obj.getCost();
                obj.results.costList(i) = cost;
                obj.results.evalList(i) = nl.coenvl.sam.ExperimentControl.getNumberEvals();
                obj.results.msgList(i) = nl.coenvl.sam.MailMan.getTotalSentMessages();
                obj.results.timeList(i) = toc(obj.t_start);
                
                % If a better solution is found, reset countDown
                if cost < bestSolution
                    countDown = obj.nStableIterations;
                    bestSolution = cost;
                end
            end
            fprintf(' done\n')
            
            % The solver is not iterative, but may take a while to complete
            obj.pauseUntilVariablesAreSet();
        end % RUNITERSOLVER
    end % Protected methods
end
