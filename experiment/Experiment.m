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
        
        % Whether the solvers to use are rooted or not (must be properly
        % adhered to by implementing classes during init())
        useRootedSolvers@logical scalar;
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
        
        % Whether to show progress during iterations
        debug@logical scalar;
        
        % Whether to keep track of solution space exploration
        ssetrack@logical scalar;
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
    
    %% Private properties
    properties (SetAccess = private, GetAccess = private)
        % Graphics handles to show debugging
        hdebug;
        hadebug;
        
        % Solution space exploration status
        sse_status
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
            obj.useRootedSolvers = getSubOption(false, 'logical', options, 'useRootedSolvers');
            obj.debug = getSubOption(false, 'logical', options, 'debug');
            obj.ssetrack = getSubOption(false, 'logical', options, 'ssetrack');
            
            assert(~isempty(obj.initSolverType) || ~isempty(obj.iterSolverType), ...
                'EXPERIMENT:INIT:NOSOLVERTYPE', ...
                'Either set initSolverType, iterSolverType or both');
%             assert(isGraph(edges), 'EXPERIMENT:INITGRAPH:INVALIDGRAPH', ...
%                 'Input must be a Nx2 matrix of connected nodes');
%             assert(graphIsConnected(edges), ...
%                 'EXPERIMENT:INITGRAPH:NOTCONNECTED', ...
%                 'Graph must be connected');
            
            % Make sure all referenced java objects are cleared
            obj.reset();
            
            obj.graph.edges = edges;
            obj.graph.density = graphDensity(edges);
            obj.graph.size = graphSize(edges);
        end % EXPERIMENT
        
        %% RUN - Run the experiment
        function run(obj)
            % Run implementation specific initialization           
            if isempty(obj.agent)
                obj.init();
            end
            
            if isempty(obj.t_start) || obj.t_start == 0
                obj.t_start = tic; % start the clock
            end
            
            if (~obj.variable{1}.isSet())
                obj.runInitSolver();
            end
            
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
            
            if (obj.ssetrack)
                obj.results.sse_explored = size(unique(obj.sse_status, 'rows'), 1);
                obj.results.sse_size = obj.nColors ^ obj.graph.size;
            end
        end % RUN
        
        %% RUNI - Run the simulation for given number of iterations
        function runI(obj, numIters)
            ExceptOn(isempty(obj.iterSolverType), ...
                'EXPERIMENT:RUNI:NOITERSOLVER', ...
                'Cannot use this method without iterative solver');
            
            if isempty(obj.agent)
                obj.init();
            end
            
            if (~obj.variable{1}.isSet())
                obj.runInitSolver();
            end
            
            if isempty(obj.t_start) || obj.t_start == 0
                obj.t_start = tic; % start the clock
            end
            
            for i = 1:numIters
                obj.tick();
            end
        end % RUN
        
        %% RESET - Resets the experiment so it can be run again
        function reset(obj)
            obj.results = struct('numIters', 0);
            cellfun(@(x) x.reset(), obj.agent);
            cellfun(@(x) x.reset(), obj.constraintAgent);
            cellfun(@(x) x.clear(), obj.variable);
            obj.agent = {};
            obj.constraintAgent = {};
            obj.variable = {};
            obj.constraint = {};
            obj.sse_status = [];
            obj.t_start = 0;
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
        % This function should set the variable, agent, constraintAgent and
        % constraint properties of the experiment into 1xN cell arrays in
        % order to run the experiment.
        init(obj);
    end
    
    methods (Sealed, Access = protected)
        %% PAUSEUNTILVARIABLESARESET
        function pauseUntilVariablesAreSet(obj)
            t_puvas = tic;
            while toc(t_puvas) < obj.maxtime
                if all(cellfun(@(x) x.isSet(), obj.variable))
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
        function bool = doStop(obj, countDown)
            bool = false;
            % First decide based on MAX number of iterations
            if (obj.nMaxIterations > 0) && (obj.results.numIters >= obj.nMaxIterations)
                bool = true;
            end
            
            % Only then look at number of stable iterations
            if ~bool && (obj.nStableIterations > 0) && (countDown <= 0)
                bool = true;
            end
        end % DOSTOP
    end
    
    methods (Sealed)
        %% GETCOST - Obtain current experiment cost
        function cost = getCost(obj)
            cost = sum(cellfun(@(x) x.getExternalCost(), obj.constraint));
        end % GETCOST
        
        %% GETOPTIMALITY
        function k = getOptimality(obj)
            for i = 1:5
                if obj.isKOptimal(i)
                    k = i - 1;
                    return;
                end
            end
        end % GETOPTIMALITY
        
        %% ISKOPTIMAL
        function [isOpt, results] = isKOptimal(obj, k, costToBeat, tryset, valueset, results)
            if nargin < 4
                tryset = [];
                valueset = [];
                results = [];
            end
            
            if nargin < 3
                costToBeat = obj.getCost();
            end
            
            if (k == 0)
                isOpt = obj.getCost() < costToBeat;
                if isOpt
                    % Store results
                    results(end+1).idx = tryset;
                    results(end).values = valueset;
                end
                return;
            end
            
           
            if ~isempty(tryset)
                startidx = tryset(end) + 1;
            else
                startidx = 1;
            end
            
            isOpt = false;
            for i = startidx:numel(obj.variable)
                originalValue = int32(obj.variable{i}.getValue());
                iter = obj.variable{i}.iterator();
                while iter.hasNext()
                    tryValue = int32(iter.next());
                    if (originalValue == tryValue)
                        continue;
                    end
                    obj.variable{i}.setValue(tryValue);
                    if (k > 1)
                        fprintf('%s%s\n', repmat('    ', 1, k), obj.variable{i}.toString());
                    end
                    [isOpt, results] = obj.isKOptimal(k - 1, costToBeat, [tryset  i], [valueset tryValue], results);
                    if isOpt
                        break;
                    end
                end
                obj.variable{i}.setValue(originalValue);
            end
            
        end % ISKOPTIMAL
        
        %% PRINT - Show the graph as an image, or print it to file
        function print(obj, filename)
            if nargin == 1
                filename = fullfile(pwd, 'graph.png');
            end
            
            [path, ~, ext] = fileparts(filename);
            ext = ext(2:end);
            if isempty(path)
                path = pwd;
                filename = fullfile(pwd, filename);
            end
            
            if ~exist(path,'dir')
                error('No such path: %s', path);
            end
    
            colornames = {'red', 'green', 'cyan', 'blue', 'yellow', 'magenta', 'brown', 'coral', 'gold', 'black', 'azure'};
            
            if obj.nColors > numel(colornames)
                error('EXPERIMENT:TOOMANYCOLORS', ...
                    'Unable to print graphs with more than %d colors', numel(colornames));
            end
            
            gvfile = fullfile(pwd, 'test.gv');
            fid = fopen(gvfile, 'w');
            fprintf(fid, 'graph {\ngraph [K=0.03, overlap=false];\nnode [style=filled];\n');

            for i = 1:obj.graph.size
                var = obj.variable{i};
                if var.isSet()
                    fprintf(fid, ' %d [fillcolor=%s, width=.2, height=.1];\n', i, colornames{double(var.getValue())});
                else
                    fprintf(fid, ' %d [fillcolor=white, width=.2, height=.1];\n', i);
                end
            end

            for i = 1:size(obj.graph.edges,1)
                fprintf(fid, ' %d -- %d;\n', obj.graph.edges(i,1), obj.graph.edges(i,2));
            end
            fprintf(fid, '}\n');
            fclose(fid);

            searchpath = 'c:\Progra~1';
            files = dir(fullfile(searchpath, 'Graphviz*'));
            if (numel(files) == 0)
                % Also search (x86)
                searchpath = 'c:\Progra~2';
                files = dir(fullfile(searchpath, 'Graphviz*'));
            end

            if (numel(files) == 0)
                error('Unable to complete printing graph, could not find Graphviz installation');
            end
            graphvizpath = fullfile(searchpath, files(1).name, 'bin');

            % layout = fullfile(graphvizpath, 'fdp.exe');
            layout = fullfile(graphvizpath, 'sfdp.exe');
            % layout = fullfile(graphvizgraphvizpathpath, 'neato.exe');
            % layout = fullfile(graphvizpath, 'twopi.exe');

            [a,b] = system(sprintf('%s %s -T %s -Gdpi=300 -o %s', layout, gvfile, ext, filename));

            if nargin == 1
                if isempty(obj.hdebug) || ~ishandle(obj.hdebug)
                    obj.hdebug = figure('name', 'Debug graph');
                    obj.hadebug = gca(obj.hdebug);
                end
                A = imread(filename);
                imshow(A, 'Parent', obj.hadebug);
            end
        end % PRINT
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

            cellfun(@(x) x.init(), obj.agent);
            cellfun(@(x) x.init(), obj.constraintAgent);
            
            if ~isempty(obj.initSolverType)
                obj.pauseUntilVariablesAreSet();
            end
            
            if obj.debug
                obj.print();
                waitforbuttonpress;
            end
            
            % Keep track of sse_status if required
            if obj.ssetrack
                obj.sse_status(end + 1,:) = cellfun(@(x) double(x.getValue()), obj.variable);
            end
        end % RUNINITSOLVER
        
        %% RUNITERSOLVER - Run iterative part of algorithm
        function runIterSolver(obj)          
            % Start iterations
            bestSolution = inf;
            countDown = obj.nStableIterations;
            %fprintf('Iteration: ');
            while ~obj.doStop(countDown)
                countDown = countDown - 1;
                cost = obj.tick();
                
                % If a better solution is found, reset countDown
                if cost < bestSolution
                    countDown = obj.nStableIterations;
                    bestSolution = cost;
                end
            end
            %fprintf(' done\n')
            
            % The solver is not iterative, but may take a while to complete
            obj.pauseUntilVariablesAreSet();
        end % RUNITERSOLVER
        
        %% TICK - Do one iteration
        function cost = tick(obj)
            % Get the number of *this* iteration
            i = obj.results.numIters + 1;
            obj.results.numIters = i;
            
            % Display the progress
            if mod(i, 25) == 0
                fprintf(' %d', i);
            end
            if mod(i, 500) == 0
                fprintf('\n');
            end

            % Do the actual iteration
            cellfun(@(x) x.tick(), obj.agent);
            cellfun(@(x) x.tick(), obj.constraintAgent);

            % Store the intermediate results
            cost = obj.getCost();
            obj.results.cost(i) = cost;
            obj.results.evals(i) = nl.coenvl.sam.ExperimentControl.getNumberEvals();
            obj.results.msgs(i) = nl.coenvl.sam.MailMan.getTotalSentMessages();
            obj.results.time(i) = toc(obj.t_start);
            
            if obj.debug
                obj.print();
                waitforbuttonpress;
            end
            
            % Keep track of sse_status if required
            if obj.ssetrack
                obj.sse_status(end + 1,:) = cellfun(@(x) double(x.getValue()), obj.variable);
            end
        end % TICK
    end % Protected methods
end
