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
function [ isConnected, missing ] = graphIsConnected( edges )

if iscell(edges)
    edges = vertcat(edges{:});
end

% Make sure the edges are sorted
edges = sortrows(edges);

% Keep a list of the connected nodes
connectedNodes = edges(1,1);

% Keep a list of the "newly" connected nodes
newNodes = connectedNodes;

% nodes = unique(edges(:));
nodes = 1:graphSize(edges);
while ~isempty(newNodes)
    newNode = newNodes(end);
    newNodes(end) = [];
    
    if any(connectedNodes == newNode);
        % If i is connected, then so are all nodes that are connected to i
        k = edges(:,1) == newNode;
        l = edges(:,2) == newNode;
        connectK = setdiff(edges(k,2), connectedNodes);
        connectL = setdiff(edges(l,1), connectedNodes);
        
        newNodes = [newNodes; connectK(:); connectL(:)];
        connectedNodes = [connectedNodes; connectK(:); connectL(:)];
    end
    
end

missing = setdiff(nodes, connectedNodes);

isConnected = isempty(missing);