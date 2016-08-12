%% PREPARERESULTS
% *Summary of this function goes here*
%
% Detailed explanation goes here
%
%% Copyright
% * *2016 - TNO*
% * *Author*: Coen van Leeuwen
% * *Since*: January 16, 2016
% 
%% See also:
%

%% Function Definition
function matresults = prepareResults(cellresults, range)

algos = sort(fieldnames(cellresults));

for i = 1:numel(algos)
    if nargin < 2 || isempty(range)
        % range = 1:max(structfun(@(x) max(x.iterations), cellresults));
        % range = 1:max([cellresults.(algos{i}).iterations{:}]);
        range = 1:max(cellresults.(algos{i}).iterations);
    end    
    
    for e = 1:numel(cellresults.(algos{i}).iterations)
        matresults.(algos{i}).iterations(e) = cellresults.(algos{i}).iterations(e);
        matresults.(algos{i}).costs(:,e) = concatResults(cellresults.(algos{i}).costs{e}, range);
        matresults.(algos{i}).evals(:,e) = concatResults(cellresults.(algos{i}).evals{e}, range);
        matresults.(algos{i}).msgs(:,e) = concatResults(cellresults.(algos{i}).msgs{e}, range);
        matresults.(algos{i}).times(:,e) = concatResults(cellresults.(algos{i}).times{e}, range);
    end
end

end

function mat = concatResults(data, range)
    if numel(data) == 1
        mat = data;
    elseif numel(data) >= numel(range)
        mat = data(range)';
    else
        mat = ones(numel(range), 1) * data(end);
        mat(1:numel(data)) = data;
    end
end