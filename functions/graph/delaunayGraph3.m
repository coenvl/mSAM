function [ e, pos ] = delaunayGraph3(options)
%DCGGRAPH Summary of this function goes here
%   Detailed explanation goes here
% 
% dcg = DynamicColorGraph(nagents);
% 
% e = dcg.neqConstraints;

nAgents = getSubOption(uint16(10), 'uint16', options, 'nAgents');
samplemethod = getSubOption('random', 'char', options, 'sampleMethod');

clear pos;
switch samplemethod
    case 'poisson'
        pos = poissonSample3(nAgents);
    case 'random'
        pos = rand(nAgents,3);
    otherwise
        pos = rand(nAgents,3);
end

if verLessThan('matlab', '8.4')
    DT = DelaunayTri(pos(:,1),pos(:,2), pos(:,3));    %#ok<DDELTRI>
else
    DT = delaunayTriangulation(pos(:,1),pos(:,2), pos(:,3));
end

e = DT.edges;

end

