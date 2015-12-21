%% GRAPHDENSITY
% *Summary of this function goes here*
%
% Detailed explanation goes here
%
%% Copyright
% * *2015 - TNO*
% * *Author*: Coen van Leeuwen
% * *Since*: July 09, 2015
% 
%% See also:
%

%% Function Definition
function density = graphDensity( edges )

% Get the number of nodes
n = numel(unique(edges));

% Calculate the maximum possible number of edges
maxEdges = 0.5 * n * (n - 1);

% Calculate the density of the graph
density = size(edges, 1) / maxEdges;