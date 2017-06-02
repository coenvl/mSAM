%% GRAPHISCONNECTED
% *Summary of this function goes here*
%
% Detailed explanation goes here
%
%% Copyright
% * *2016 - TNO*
% * *Author*: Coen van Leeuwen
% * *Since*: July 29, 2016
% 
%% See also:
%

%% Function Definition
function [ isConnected, missing ] = higherOrderGraphIsConnected(edges, knownGraphSize)

% First find out what nodes exist...
if (nargin > 1)
    nodes = 1:knownGraphSize;
else
    nodes = 1:graphSize(edges);
end

edges = [edges{:}];
size = cellfun(@(x) numel(x), edges);

% Throw away edges of size 1 (they don't contribute)
edges = edges(size > 1);

% Keep a list of the connected nodes
connectedNodes = edges{1,1}(:);

% Keep a list of the "newly" connected nodes
newNodes = connectedNodes(:);

while ~isempty(newNodes)
    newNode = newNodes(end);
    newNodes(end) = [];
    
    if any(connectedNodes == newNode)
        % If newNode is connected, then so are all nodes that are connected to it
        k = cellfun(@(x) any(x == newNode), edges);
        connectK = setdiff([edges{k}], connectedNodes);
        
        newNodes = [newNodes; connectK(:)];
        connectedNodes = [connectedNodes; connectK(:)];
    end
    
end

missing = setdiff(nodes, connectedNodes);

isConnected = isempty(missing);