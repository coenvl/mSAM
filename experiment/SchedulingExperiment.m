%% SCHEDULINGEXPERIMENT
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
classdef SchedulingExperiment < Experiment
    
    properties (Constant)
        preference_property = 'preference';
    end
    
    methods
        %% SCHEDULINGEXPERIMENT - Experiment constructor
        function obj = SchedulingExperiment(edges, options)
            obj = obj@Experiment(edges, options);
            obj.graph.nmeetings = numel(edges);
        end % SCHEDULINGEXPERIMENT
        
        %% INIT - Experiment initialization
        function init(obj)
            obj.initVariables();
            obj.initEqConstraints();
            obj.initNeqConstraints();
            
            reduce = @(x) x(~cellfun(@isempty, x));
            obj.variable = reduce(obj.variable);
            obj.agent = reduce(obj.agent);
            obj.constraintAgent = reduce(obj.constraintAgent);
            obj.constraint = reduce(obj.constraint);
        end % INIT
    end
    
    methods (Access = private)
        %% INITVARIABLES - Initialize variables and agents
        function initVariables(obj)
            idx = 0;
            for m = 1:obj.graph.nmeetings
                graph = obj.graph.edges{m};
                for i = unique(graph(:))'
                    idx = idx + 1;
                    varName = sprintf('variable%05d', idx);
                    agentName = sprintf('agent%05d', idx);
                    
                    obj.variable{m,i} = nl.coenvl.sam.variables.IntegerVariable(int32(1), int32(obj.nColors), varName);
                    newagent = nl.coenvl.sam.agents.VariableAgent(obj.variable{m,i}, agentName);
                    obj.agent{m,i} = newagent;
                    
                    if ~isempty(obj.initSolverType)
                        newagent.setInitSolver(feval(obj.initSolverType, newagent));
                    end
                    if ~isempty(obj.iterSolverType)
                        newagent.setIterativeSolver(feval(obj.iterSolverType, newagent));
                    end
                    
                    % Set preference property
                    newagent.set(SchedulingExperiment.preference_property, obj.agentProps(i).(SchedulingExperiment.preference_property));
                end
            end
        end % INITVARIABLES
        
        %% INITEQCONSTRAINTS
        % Add constraints saying one meeting must have the same time for
        % all agents, but agents may have different preferences as to when
        % to have it.
        function initEqConstraints(obj)
            peqc = 'nl.coenvl.sam.constraints.PreferentialEqualityConstraint';
            constraintCost = 1e3;
            
            idx = 0;
            % First add constraints within meetings
            for m = 1:obj.graph.nmeetings
                graph = obj.graph.edges{m};
                for i = 1:size(graph,1)
                    a = graph(i,1);
                    b = graph(i,2);
                    idx = idx + 1;
                    
                    obj.constraint{idx} = feval(peqc, obj.variable{m,a}, obj.variable{m,b}, obj.agent{m,a}.get(SchedulingExperiment.preference_property), obj.agent{m,b}.get(SchedulingExperiment.preference_property), constraintCost);
                    
                    obj.agent{m,a}.addConstraint(obj.constraint{idx});
                    obj.agent{m,b}.addConstraint(obj.constraint{idx});
                    
                    if ~isempty(strfind(obj.iterSolverType, 'MaxSum'))
                        functionSolverType = getSolverCounterPart(obj.iterSolverType);
            
                        % Create constraint agent
                        agentName = sprintf('constraint%05d', i);
                        obj.constraintAgent{idx} = nl.coenvl.sam.agents.ConstraintAgent(agentName, obj.constraint{idx}, obj.variable{m,a}, obj.variable{m,b});
                        obj.constraintAgent{idx}.setSolver(feval(functionSolverType, obj.constraintAgent{idx}));
                        
                        % Set constraint agent address as targets
                        obj.agent{m,a}.addFunctionAddress(obj.constraintAgent{idx}.getID());
                        obj.agent{m,b}.addFunctionAddress(obj.constraintAgent{idx}.getID());
                    end
                end
            end
        end % INITEQCONSTRAINTS
        
        %% INITNEQCONSTRAINTS
        % Add constraints saying that one agent cannot be in two
        % meetings at the same time.
        function initNeqConstraints(obj)
            neqc = 'nl.coenvl.sam.constraints.InequalityConstraint';
            constraintCost = 1e9;
            
            idx = numel(obj.constraint);
            % Next add constraints between (same) agents
            for a = 1:size(obj.agent,2)
                for m1 = 1:obj.graph.nmeetings
                    if ~isempty(obj.agent{m1,a})
                        for m2 = (m1 + 1):obj.graph.nmeetings
                            if ~isempty(obj.agent{m2,a})
                                % add constraints between agent {m1,a} and {m2,a}
                                idx = idx + 1;
                                
                                obj.constraint{idx} = feval(neqc, obj.variable{m1,a}, obj.variable{m2,a}, constraintCost);
                                
                                obj.agent{m1,a}.addConstraint(obj.constraint{idx});
                                obj.agent{m2,a}.addConstraint(obj.constraint{idx});
                                
                                if ~isempty(strfind(obj.iterSolverType, 'MaxSum'))
                                    % Create constraint agent
                                    functionSolverType = getSolverCounterPart(obj.iterSolverType);
                                    
                                    agentName = sprintf('constraint%05d', idx);
                                    obj.constraintAgent{idx} = nl.coenvl.sam.agents.ConstraintAgent(agentName, obj.constraint{idx}, obj.variable{m1,a}, obj.variable{m2,a});
                                    obj.constraintAgent{idx}.setSolver(feval(functionSolverType, obj.constraintAgent{idx}));
                                    
                                    % Set constraint agent address as targets
                                    obj.agent{m1,a}.addFunctionAddress(obj.constraintAgent{idx}.getID());
                                    obj.agent{m2,a}.addFunctionAddress(obj.constraintAgent{idx}.getID());
                                end
                            end
                        end
                    end
                end
            end
        end % INITNEQCONSTRAINTS
    end % Private methods
    
end
