function [ e ] = scalefreeGraph(options)
%SCALEFREEGRAPH Summary of this function goes here
%   Detailed explanation goes here

nagents = getSubOption(uint16(10), 'uint16', options, 'nAgents');
maxlinks = getSubOption(uint16(3), 'uint16', options, 'maxLinks');
initialsize = getSubOption(uint16(3), 'uint16', options, 'initialsize');

net = SFNG(nagents,maxlinks,not(eye(initialsize))); %[0 1 1;1 0 1; 1 1 0]);
[i,j] = find(net);
k = i < j;
e = sortrows([i(k) j(k)]);

end

