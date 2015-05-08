function e = delaunayGraph(options)
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
        pos = poissonSample(nAgents);
    case 'random'
        pos = rand(nAgents,2);
    otherwise
        pos = rand(nAgents,2);
end

if verLessThan('matlab', '8.4')
    DT = DelaunayTri(pos(:,1),pos(:,2));    %#ok<DDELTRI>
else
    DT = delaunayTriangulation(pos(:,1),pos(:,2));
end

e = DT.edges;

end

