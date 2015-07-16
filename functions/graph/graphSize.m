%% GRAPHSIZE
% *Summary of this function goes here*
%
% Detailed explanation goes here
%
%% Copyright
% * *2015 - TNO*
% * *Author*: Coen van Leeuwen
% * *Since*: July 10, 2015
% 
%% See also:
%

%% Function Definition
function [ size ] = graphSize( edges )

size = numel(unique(edges(:)));