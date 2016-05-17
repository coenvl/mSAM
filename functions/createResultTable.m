%% CREATERESULTTABLE
% *Summary of this function goes here*
%
% Detailed explanation goes here
%
%% Copyright
% * *2016 - TNO*
% * *Author*: Coen van Leeuwen
% * *Since*: January 30, 2016
% 
%% See also:
%

%% Function Definition
function varargout = createResultTable(exp)

if numel(exp) == 1
    exp = struct('results', exp);
end

for i = 1:numel(exp)
    A(i) = analyzeResults(exp(i).results);
end

header = sprintf('\\toprule\n\\hline % -10s & % -7s & % -7s & % -7s & % -7s & % -7s\\\\ \n\\midrule \n', 'Algorithm', 'I', 'S', 'M', 'E', 'T');

algos = fieldnames(A)';
algoStr = cell(size(algos));
for i = 1:numel(algos)
    algoResults = [A.(algos{i})];
    iterStr = sprintf(repmat('&% -7d ', 1, numel(algoResults)), round([algoResults.iterations]));
    costStr = sprintf(repmat('&% -7d ', 1, numel(algoResults)), round([algoResults.costs]));
    msgsStr = sprintf(repmat('&% -7d ', 1, numel(algoResults)), round([algoResults.msgs] / 1000));
    evalStr = sprintf(repmat('&% -7d ', 1, numel(algoResults)), round([algoResults.evals] / 1000));
    timeStr = sprintf(repmat('&% -6.1f ', 1, numel(algoResults)), [algoResults.times]);
    algoStr{i} = sprintf('\\hline % -10s %s %s %s %s %s', algos{i}, iterStr, costStr, msgsStr, evalStr, timeStr);
end

str = [header sprintf('%s \\\\\n', algoStr{:}) '\bottomrule'];

if nargout > 0
    varargout{1} = str;
else
    fprintf('%s\n', str);
end