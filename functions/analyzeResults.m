%% ANALYZERESULTS
% *Summary of this function goes here*
%
% Detailed explanation goes here
%
%% Copyright
% * *2016 - TNO*
% * *Author*: Coen van Leeuwen
% * *Since*: January 15, 2016
% 
%% See also:
%

%% Function Definition
function out = analyzeResults(results)

algos = sort(fieldnames(results));

% Compute for all algorithms
for i = 1:numel(algos)
    
    % Initialize one result for each experiment
    y = results.(algos{i}).costs;
    iters = nan(1,numel(y));
    costs = nan(1,numel(y));
    msgs = nan(1,numel(y));
    evals = nan(1,numel(y));
    
    % For all experiments
    for e = 1:numel(y);
        if all(isnan(y{e})); continue; end
        
        n = find(y{e} < (1.01*min(y{e})), 1, 'first');
        if isempty(n); n = numel(y); end
       
        iters(e) = n;
        
        costs(e) = results.(algos{i}).costs{e}(n);
        msgs(e) = results.(algos{i}).msgs{e}(n);
        evals(e) = results.(algos{i}).evals{e}(n);
    end  

    out.(algos{i}).iterations = ceil(nanmean(iters));
    out.(algos{i}).costs = nanmean(costs);
    out.(algos{i}).msgs = nanmean(msgs);
    out.(algos{i}).evals = nanmean(evals);
%     
%     fprintf('Results for %s:\n', algos{i});
%     fprintf('Iterations: %d\n', ceil(nanmean(iters)));
%     fprintf('Costs: %0.1f\n', nanmean(costs));
%     fprintf('Messages: %0.1f\n', nanmean(msgs));
%     fprintf('Evaluations: %0.1f\n\n', nanmean(evals));
    
end