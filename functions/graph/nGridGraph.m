function e = nGridGraph(options)
%NGRIDGRAPH Summary of this function goes here
%   Detailed explanation goes here
% 

nAgents = double(getSubOption(uint16(20), 'uint16', options, 'nAgents'));
nDims = double(getSubOption(uint16(2), 'uint16', options, 'nDims'));

dSize = repmat(uint16(floor(nAgents^(1/nDims))),1,nDims);
for i = 1:(nDims-1)
    dSize(i+1) = floor((nAgents / prod(dSize(1:i))) ^ (1 / (nDims - i)));
end

gridSize = getSubOption(dSize, 'uint16', options, 'gridSize')
doWrap = getSubOption('', char, options, 'doWrap');

%%

e = uint16([]);
for d = 1:numel(gridSize);
    l = gridSize(d);
    n = prod(gridSize(1:d-1));
    
    % Duplicate previous edges in new dim
    [A,B] = meshgrid(e,n*(0:(l-1)));
    e = reshape(A(:)+B(:),[],2); 
    
    % Connect between dimensions
    from = 1:(n*(l-1));
    to = from + n;
    e = [e; from' to'];
end
e = sortrows(e);

%%

% w = gridSize(1);
% h = gridSize(2);
% nAgents = w*h;
% 
% % Horizontal connections
% s = bsxfun(@plus, (1:w:nAgents)', 0:(w-2));
% e = [s(:) s(:)+1];
% 
% if strfind(doWrap, 'x')
%     e = [e; (w:w:nAgents)' (1:w:nAgents)'];
% end
% 
% % Vertical connections
% s = (1:(w*(h-1)))';
% e = [e; s s+w];
% 
% if strfind(doWrap, 'y')
%     e = [e; 1+((w*(h-1)):(nAgents-1))' (1:w)'];
% end
% 
% % Sort them nicely
% e = sortrows(e);

end

