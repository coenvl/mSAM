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
            
            
        end
    end
    
end
