function e = gridGraph(options)
%GRIDGRAPH Summary of this function goes here
%   Detailed explanation goes here
% 

nAgents = getSubOption(uint16(20), 'uint16', options, 'nAgents');
dWidth = uint16(ceil(sqrt(double(nAgents))));
dHeight = nAgents / dWidth;
gridSize = getSubOption([dWidth dHeight], 'uint16', options, 'gridSize');
doWrap = getSubOption('', char, options, 'doWrap');

%%
w = gridSize(1);
h = gridSize(2);
nAgents = w*h;

% Horizontal connections
s = bsxfun(@plus, (1:w:nAgents)', 0:(w-2));
e = [s(:) s(:)+1];

if strfind(doWrap, 'x')
    e = [e; (w:w:nAgents)' (1:w:nAgents)'];
end

% Vertical connections
s = (1:(w*(h-1)))';
e = [e; s s+w];

if strfind(doWrap, 'y')
    e = [e; 1+((w*(h-1)):(nAgents-1))' (1:w)'];
end

% Sort them nicely
e = sortrows(e);

end

