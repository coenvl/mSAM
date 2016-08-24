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
classdef (Abstract) Experiment < handle
    
    %% Public properties
    properties
        % The type of solver to use for initialization
        initSolverType@char matrix;
        
        % The type of solver to use for iteration
        iterSolverType@char matrix;
        
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
    
    %% Final protected properties (settings)
    properties (SetAccess = immutable, GetAccess = protected)
        % The type of constraints to use between variables
        constraintType@char matrix;
        
        % Arguments to initialize the constraints with
        constraintArgs@cell;
        
        % Properties to provide to the agents
        agentProps@struct;
        
        % Number of colors (domain size)
        nColors@uint16 scalar;
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
           
        end % EXPERIMENT
        
        %% RUN - Run the experiment
        function run(obj)
            % Run implementation specific initialization           
            obj.init();
            
            obj.runInitSolver();
            
            if ~isempty(obj.iterSolverType)
                obj.runIterSolver();
            else
                obj.pauseUntilVariablesAreSet();
                obj.results.time = toc(obj.t_start); 
                obj.results.numIters = 1;
                obj.results.cost = obj.getCost();
                obj.results.evals = nl.coenvl.sam.ExperimentControl.getNumberEvals();
                obj.results.msgs = nl.coenvl.sam.MailMan.getTotalSentMessages();
            end
        end % RUN
        
        %% RESET - Resets the experiment so it can be run again
        function reset(obj)
            obj.results = struct();
            cellfun(@(x) x.reset(), obj.agent);
            cellfun(@(x) x.reset(), obj.constraintAgent);
            cellfun(@(x) x.clear(), obj.variable);
            obj.agent = {};
            obj.constraintAgent = {};
            obj.variable = {};
            obj.constraint = {};
            nl.coenvl.sam.ExperimentControl.ResetExperiment();
        end % RESET
        
        %% DELETE - Object constructor
        function delete(obj)
            % Experiment destructor
            obj.reset();
        end % DELETE
        
        %% SETINITSOLVERTYPE
        function set.initSolverType(obj, type)
            fprintf('Updating initialization solver to %s\n', type);
            obj.initSolverType = type;
            obj.reset();
        end % SETINITSOLVERTYPE
        
        %% SETITERSOLVERTYPE
        function set.iterSolverType(obj, type)
            fprintf('Updating iterative solver to %s\n', type);
            obj.iterSolverType = type;
            obj.reset();
        end % SETITERSOLVERTYPE
        
    end
    
    methods (Abstract)
        init(obj);
    end
    
    methods (Sealed, Access = protected)
        %% PAUSEUNTILVARIABLESARESET
        function pauseUntilVariablesAreSet(obj)
            t_puvas = tic;
            while toc(t_puvas) < obj.maxtime
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
    end
    
    %% Protected methods
    methods (Sealed, Access = protected)       
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
        
        %% RUNINITSOLVER - Run initialization part of algorithm
        function runInitSolver(obj)
            a = obj.agent{randi(numel(obj.agent))};
            a.set(nl.coenvl.sam.solvers.CoCoSolver.ROOTNAME_PROPERTY, true);
            
            obj.t_start = tic; % start the clock
            
            cellfun(@(x) x.init(), obj.agent);
            cellfun(@(x) x.init(), obj.constraintAgent);
            
            if ~isempty(obj.initSolverType)
                obj.pauseUntilVariablesAreSet();
            end
        end % RUNINITSOLVER
        
        %% RUNITERSOLVER - Run iterative part of algorithm
        function runIterSolver(obj)          
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
                if mod(i, 500) == 0
                    fprintf('\n');
                end
                
                cellfun(@(x) x.tick(), obj.agent);
                cellfun(@(x) x.tick(), obj.constraintAgent);
                
                cost = obj.getCost();
                obj.results.cost(i) = cost;
                obj.results.evals(i) = nl.coenvl.sam.ExperimentControl.getNumberEvals();
                obj.results.msgs(i) = nl.coenvl.sam.MailMan.getTotalSentMessages();
                obj.results.time(i) = toc(obj.t_start);
                
                % If a better solution is found, reset countDown
                if cost < bestSolution
                    countDown = obj.nStableIterations;
                    bestSolution = cost;
                end
            end
            fprintf(' done\n')
            
            obj.results.numIters = i;
            
            % The solver is not iterative, but may take a while to complete
            obj.pauseUntilVariablesAreSet();
        end % RUNITERSOLVER
    end % Protected methods
end
