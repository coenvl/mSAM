function e = delaunayGraph(options)
%DCGGRAPH Summary of this function goes here
%   Detailed explanation goes here
% 
% dcg = DynamicColorGraph(nagents);
% 
% e = dcg.neqConstraints;

nagents = getSubOption(uint16(10), 'uint16', options, 'nAgents');

if verLessThan('matlab', '8.4')
    DT = DelaunayTri(rand(nagents, 1),rand(nagents, 1));    %#ok<DDELTRI>
else
    DT = delaunayTriangulation(rand(nagents, 1),rand(nagents, 1));
end

e = DT.edges;

end

