%% SwitchingSolverGCE
% *Summary of this class goes here*
%
% Detailed explanation goes here
%
%% Copyright
% * *2016 - TNO*
% * *Author*: Coen van Leeuwen
% * *Since*: October 7, 2016 (*30*)
%
%% See also:
%

%% Class definition:
classdef SwitchingSolverGCE < GraphColoringExperiment
    
    %% Public methods
    methods
        
        %% SWITCHINGSOLVERGCE - Object constructor
        function obj = SwitchingSolverGCE(edges, options)
            % Parse all options
            obj = obj@GraphColoringExperiment(edges, options);
            
            ExceptOn(~isempty(strfind(obj.iterSolverType, 'MaxSum')), ...
                'SWITCHINGSOLVERGCE:UNDEFINED', ...
                'This experiment is not suited for constraint agents');
        end % SWITCHINGSOLVERGCE
        
        %% INIT - build the variables and agents
        function init(obj)
            obj.initVariables();
            obj.initConstraints();
            
            ExceptOn(~isempty(strfind(obj.iterSolverType, 'MaxSum')), ...
                'SWITCHINGSOLVERGCE:INIT:UNDEFINED', ...
                'This experiment is not suited for constraint agents');
            
            obj.assignAgentProperties();
        end % INIT
        
        %% SWITCHSOLVER - Switch solver
        function switchSolver(obj, newIterSolverType)
            fprintf('Switching to %s', newIterSolverType); 
            for i = 1:numel(obj.agent)
                obj.agent{i}.setIterativeSolver(feval(newIterSolverType, obj.agent{i}));
            end
        end % SWITCHSOLVER
    end
    
    %% Protected methods
    methods (Access = protected)

        %% INITCONSTRAINTAGENTS - Add constraint agents (MAXSUM only)
        function initConstraintAgents(~)
            error('SWITCHINGSOLVERGCE:INITCONSTRAINTAGENTS:UNDEFINED', ...
                'This experiment is not suited for constraint agents');
        end % INITCONSTRAINTAGENT
    end % Private methods
end
