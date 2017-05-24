%% WPTEXPERIMENT
% *Wireless Power Transfer Experiment*
%
% In the wireless power transfer situation the variables are the 
% transmitter powers and the regular constraints represent the receivers. 
% Their cost function represents the received power as a negative number. 
% The cost of the sensors are represented as "hard" constraints in this 
% class as additional constraints (and/or) constraint agents.
%
%% Copyright
% * *2016 - TNO*
% * *Author*: Coen van Leeuwen
% * *Since*: August 26, 2016
%
%% See also:
%

%% Class definition:
classdef WPTExperiment < Experiment
       
    properties (Access = private)
        % Arguments to initialize the sensor constraints with
        sensorConstraintArgs@cell;        
        
        % Constraints to constrain the battery lifetime
        sensorConstraint@cell;
        
        % Only for maxSum
        sensorConstraintAgent@cell;
        
        % Measurement noise defines the magnitude of the error of
        % measurements
        measurementNoise@double scalar;
    end
    
    methods
        %% WPTEXPERIMENT - Experiment constructor
        function obj = WPTExperiment(edges, options)
            obj = obj@Experiment(edges{1}, options);
            obj.sensorConstraintArgs = getSubOption({}, 'cell', options, 'sensorConstraint', 'arguments');
            obj.measurementNoise = getSubOption(0.0, 'double', options, 'measurementNoise');
                        
            obj.graph.sensorEdges = edges{2};           
        end % WPTEXPERIMENT
        
        %% INIT - Experiment initialization
        function init(obj)
            obj.initVariables();

            obj.initReceiverConstraints();
            obj.initSensorConstraints();
            if ~isempty(strfind(obj.iterSolverType, 'MaxSum'))
                obj.initReceiverConstraintAgents();
                obj.initSensorConstraintAgents();
            end
            
            % Concatenate wpt properties into properties
            obj.constraint = [obj.constraint obj.sensorConstraint];
            obj.constraintAgent = [obj.constraintAgent obj.sensorConstraintAgent];
        end % INIT
    end
    
    methods (Access = private)
        %% INITVARIABLES - Initialize receiver variables and agents
        function initVariables(obj)
            for i = 1:max(max([obj.graph.edges{:} obj.graph.sensorEdges{:}]))
                varName = sprintf('transmitter%05d', i);
                agentName = sprintf('agent%05d', i);
                
                obj.variable{i} = nl.coenvl.sam.variables.FixedPrecisionVariable(0, 10, 10 / double(obj.nColors), varName);
                obj.agent{i} = nl.coenvl.sam.agents.VariableAgent(obj.variable{i}, agentName, nl.coenvl.sam.agents.SolverAgent.MULTI_THREADED);
                
                obj.variable{i}.set('position', obj.agentProps(i).position);
                
                if ~isempty(obj.initSolverType)
                    obj.agent{i}.setInitSolver(feval(obj.initSolverType, obj.agent{i}));
                end
                if ~isempty(obj.iterSolverType)
                    obj.agent{i}.setIterativeSolver(feval(obj.iterSolverType, obj.agent{i}));
                end
            end
        end % INITVARIABLES
        
        %% INITRECEIVERCONSTRAINTS - Add constraints to the graph
        function initReceiverConstraints(obj)
            constraintType = 'nl.coenvl.sam.wpt.WPTReceiverConstraint';
            for i = 1:numel(obj.graph.edges)
                edge = obj.graph.edges{i};
                
                if (~isempty(obj.measurementNoise))
                    obj.constraint{i} = feval(constraintType, obj.constraintArgs{i}, 1 + (randn * obj.measurementNoise));
                else
                    obj.constraint{i} = feval(constraintType, obj.constraintArgs{i});
                end
                
                for j = 1:numel(edge)
                    a = edge(j);
                    obj.constraint{i}.addVariable(obj.variable{a});
                    obj.agent{a}.addConstraint(obj.constraint{i});
                end
            end
        end % INITCONSTRAINTS
        
        %% INITRECEIVERCONSTRAINTAGENTS - Add receiver constraint agents
        function initReceiverConstraintAgents(obj)
            functionSolverType = getSolverCounterPart(obj.iterSolverType);
            
            for i = 1:numel(obj.graph.edges)
                edge = obj.graph.edges{i};
                involvedVars = [obj.variable{edge}];
                
                % Create constraint agent
                agentName = sprintf('constraint%05d', i);
                
                obj.constraintAgent{i} = nl.coenvl.sam.agents.HigherOrderConstraintAgent(agentName, obj.constraint{i}, involvedVars);
                obj.constraintAgent{i}.setSolver(feval(functionSolverType, obj.constraintAgent{i}));

                % Set constraint agent address as targets
                for a = edge(:)'
                     obj.agent{a}.addFunctionAddress(obj.constraintAgent{i}.getID());
                end
            end
        end % INITSECONSTRAINTAGENT
        
        %% INITSENSORCONSTRAINTS - Add sensor constraints to the graph
        function initSensorConstraints(obj)
            constraintType = 'nl.coenvl.sam.wpt.WPTSensorConstraint';
            for i = 1:numel(obj.graph.sensorEdges)
                edge = obj.graph.sensorEdges{i};
                
                if (~isempty(obj.measurementNoise))
                    obj.sensorConstraint{i} = feval(constraintType, obj.sensorConstraintArgs{i}, 1 + (randn * obj.measurementNoise));
                else
                    obj.sensorConstraint{i} = feval(constraintType, obj.sensorConstraintArgs{i});
                end
                
                for j = 1:numel(edge) 
                    a = edge(j);
                    obj.sensorConstraint{i}.addVariable(obj.variable{a});
                    obj.agent{a}.addConstraint(obj.sensorConstraint{i});
                end
            end
        end % INITSENSORCONSTRAINTS
        
        %% INITSENSORCONSTRAINTAGENTS - Add sensor constraint agents
        function initSensorConstraintAgents(obj)
            functionSolverType = getSolverCounterPart(obj.iterSolverType);
            
            for i = 1:numel(obj.graph.sensorEdges)
                edge = obj.graph.sensorEdges{i};
                involvedVars = [obj.variable{edge}];
                
                % Create constraint agent
                agentName = sprintf('sensorConstraint%05d', i);
                obj.sensorConstraintAgent{i} = nl.coenvl.sam.agents.HigherOrderConstraintAgent(agentName, obj.sensorConstraint{i}, involvedVars);
                obj.sensorConstraintAgent{i}.setSolver(feval(functionSolverType, obj.sensorConstraintAgent{i}));
                
                % Set constraint agent address as targets
                for a = edge(:)'
                    obj.agent{a}.addFunctionAddress(obj.sensorConstraintAgent{i}.getID());
                end
            end
        end % INITSENSORCONSTRAINTAGENTS
    end % Private methods
    
end
