function [ e ] = randomGraph(options)
%RANDOMGRAPH Summary of this function goes here
%   Detailed explanation goes here

nagents = getSubOption(uint16(10), 'uint16', options, 'nAgents');
density = getSubOption(.1, 'double', options, 'density');

% compute ALL possible connections
[A,B] = meshgrid(1:nagents);
e_all = [A(:) B(:)];
e_all = e_all(e_all(:,1) < e_all(:,2),:); %because undirected

% Now select a random set
k = randperm(size(e_all,1));
e = e_all(k(1:ceil(density*size(e_all,1))),:);

% Add any potentially left out nodes
not_in = setdiff(1:nagents, e);
not_partners = e(randi(numel(e),1,numel(not_in)));
e = [e; [not_in' not_partners']];

% And nicely sort them
e = sortrows(e);

end

