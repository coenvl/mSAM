function [ e, pos ] = distanceGraph3(options)

nAgents = getSubOption(uint16(10), 'uint16', options, 'nAgents');
samplemethod = getSubOption('random', 'char', options, 'sampleMethod');
scale = getSubOption(10, 'double', options, 'scale');
distance = getSubOption(7.5, 'double', options, 'maxDist');

clear pos;
switch samplemethod
    case 'poisson'
        pos = poissonSample3(nAgents);
    case 'random'
        pos = rand(nAgents,3);
    otherwise
        pos = rand(nAgents,3);
end

pos = pos - min(min(pos));
pos = pos ./ max(max(pos));
pos = pos * scale;

D = squareform(pdist(pos));
D(D == 0) = Inf;
[I,J] = find(D < distance);

e = sortrows([I(I < J) J(I < J)]);

end

