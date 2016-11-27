%% GraphColoringExperiment
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
classdef GraphColoringExperiment < Experiment
    
    %% Public methods
    methods
        
        %% GRAPHCOLORINGEXPERIMENT - Object constructor
        function obj = GraphColoringExperiment(edges, options)
            % Parse all options
            obj = obj@Experiment(edges, options);
        end % GRAPHCOLORINGEXPERIMENT
        
        %% INIT - build the variables and agents
        function init(obj)
            obj.initVariables();
            obj.initConstraints();
            
            if ~isempty(strfind(obj.iterSolverType, 'MaxSum'))
                obj.initConstraintAgents();
            end
            
            obj.assignAgentProperties();
        end % INIT

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
            functionSolverType = getSolverCounterPart(obj.iterSolverType);
            
            for i = 1:size(obj.graph.edges,1)
                a = obj.graph.edges(i,1);
                b = obj.graph.edges(i,2);
                
                % Create constraint agent
                agentName = sprintf('constraint%05d', i);
                obj.constraintAgent{i} = nl.coenvl.sam.agents.BinaryConstraintAgent(agentName, obj.constraint{i}, obj.variable{a}, obj.variable{b});
                obj.constraintAgent{i}.setSolver(feval(functionSolverType, obj.constraintAgent{i}));
                
                % Set constraint agent address as targets
                obj.agent{a}.addFunctionAddress(obj.constraintAgent{i}.getID());
                obj.agent{b}.addFunctionAddress(obj.constraintAgent{i}.getID());
            end
        end % INITCONSTRAINTAGENT
    end % Protected methods
end
