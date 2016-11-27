%% SENSORNETWORKEXPERIMENT
% *Summary of this class goes here*
%
% Detailed explanation goes here
%
%% Copyright
% * *2016 - TNO*
% * *Author*: Coen van Leeuwen
% * *Since*: August 26, 2016
%
%% See also:
%

%% Class definition:
classdef SensorNetworkExperiment < Experiment
       
    properties (Access = private)
        % Variables to store the battery limit (fixed)
        batteryVar@cell;
        
        % Constraints to constrain the battery lifetime
        batteryConstraint@cell;
        
        % Only for maxSum
        batteryConstraintAgent@cell;
    end
    
    methods
        %% SENSORNETWORKEXPERIMENT - Experiment constructor
        function obj = SensorNetworkExperiment(edges, options)
            obj = obj@Experiment(edges, options);
        end % SENSORNETWORKEXPERIMENT
        
        %% INIT - Experiment initialization
        function init(obj)
            obj.initVariables();
            
            obj.initStateEstimationConstraints();
            obj.initBatteryLifeConstraints();
            
            if ~isempty(strfind(obj.iterSolverType, 'MaxSum'))
                obj.initSEConstraintAgents();
                obj.initBLConstraintAgents();
            end
            
            % Concatenate battery-related properties into properties
            obj.variable = [obj.variable obj.batteryVar];
            obj.constraint = [obj.constraint obj.batteryConstraint];
            obj.constraintAgent = [obj.constraintAgent obj.batteryConstraintAgent];
        end % INIT
    end
    
    methods (Access = private)
        %% INITVARIABLES - Initialize variables and agents
        function initVariables(obj)
            for i = 1: obj.graph.size
                varName = sprintf('variable%05d', i);
                batteryVarName = sprintf('batteryVariable%05d', i);
                agentName = sprintf('agent%05d', i);
                
                obj.variable{i} = nl.coenvl.sam.variables.IntegerVariable(int32(1), int32(obj.nColors), varName);
                obj.batteryVar{i} = nl.coenvl.sam.variables.IntegerVariable(int32(1), int32(obj.nColors), batteryVarName);
                obj.agent{i} = nl.coenvl.sam.agents.VariableAgent(obj.variable{i}, agentName);
                
                if ~isempty(obj.initSolverType)
                    obj.agent{i}.setInitSolver(feval(obj.initSolverType, obj.agent{i}));
                end
                if ~isempty(obj.iterSolverType)
                    obj.agent{i}.setIterativeSolver(feval(obj.iterSolverType, obj.agent{i}));
                end
                
                % Per default make it so it HAS to be less than this value
                obj.batteryVar{i}.setValue(int32(obj.nColors * .5));
            end
        end % INITVARIABLES
        
        %% INITSTATEESTIMATIONCONSTRAINTS - Add constraints to the graph
        function initStateEstimationConstraints(obj)
            for i = 1:size(obj.graph.edges,1)
                a = obj.graph.edges(i,1);
                b = obj.graph.edges(i,2);

                ExceptOn(numel(obj.constraintArgs) ~= 1, ...
                    'SENSORNETWORKEXPERIMENT:INIT:INCORRECTARGUMENTCOUNT', ...
                    'Incorrect number of constraint arguments, must be exactly 1 NxN array');
                    
                obj.constraint{i} = feval(obj.constraintType, obj.variable{a}, obj.variable{b}, ...
                    obj.constraintArgs{1}, obj.constraintArgs{1});
 
                obj.agent{a}.addConstraint(obj.constraint{i});
                obj.agent{b}.addConstraint(obj.constraint{i});
            end
            
        end % INITCONSTRAINTS
        
        %% INITSECONSTRAINTAGENTS - Add constraint agents (MAXSUM only)
        function initSEConstraintAgents(obj)
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
        end % INITSECONSTRAINTAGENT
        
        %% INITBATTERYLIFETIMECONSTRAINTS - Add constraints to the graph
        function initBatteryLifeConstraints(obj)
            constraintType = 'nl.coenvl.sam.constraints.LessThanConstraint';
            constraintCost = 1e6;
            for i = 1:obj.graph.size
                obj.batteryConstraint{i} = feval(constraintType, obj.variable{i}, ...
                    obj.batteryVar{i}, constraintCost);
 
                obj.agent{i}.addConstraint(obj.batteryConstraint{i});
            end
        end % INITCONSTRAINTS
        
        %% INITBLCONSTRAINTAGENTS - Add constraint agents (MAXSUM only)
        function initBLConstraintAgents(obj)
            functionSolverType = getSolverCounterPart(obj.iterSolverType);
            
            for i = 1:obj.graph.size                
                % Create constraint agent
                agentName = sprintf('lifetimeConstraint%05d', i);
                obj.batteryConstraintAgent{i} = nl.coenvl.sam.agents.BinaryConstraintAgent(agentName, obj.batteryConstraint{i}, obj.variable{i}, obj.batteryVar{i});
                obj.batteryConstraintAgent{i}.setSolver(feval(functionSolverType, obj.batteryConstraintAgent{i}));
                
                % Set constraint agent address as targets
                obj.agent{i}.addFunctionAddress(obj.batteryConstraintAgent{i}.getID());
            end
        end % INITSECONSTRAINTAGENT
    end % Private methods
    
end
