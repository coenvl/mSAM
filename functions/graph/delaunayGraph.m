function e = delaunayGraph(options)
%DCGGRAPH Summary of this function goes here
%   Detailed explanation goes here
% 
% dcg = DynamicColorGraph(nagents);
% 
% e = dcg.neqConstraints;

nagents = getSubOption(uint16(10), 'uint16', options, 'nAgents');

DT = delaunayTriangulation(rand(nagents, 1),rand(nagents, 1));

e = DT.edges;

end

