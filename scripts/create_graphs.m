%% Get Default options
options = getGraphOptions();
options.plot.colors = cubehelix(7, .5, -1.5, 3, 1);
options.export.do = false;
options.export.format = 'eps';
options.export.folder = 'C:\Data\Documents\PHD\ijcai2016\images';
options.plot.errorbar = false;
% options.plot.hi_error_fun = @(x) nanmean(x,2) + nanstd(x, [], 2);
% options.plot.lo_error_fun = @(x) nanmean(x,2) - nanstd(x, [], 2);
options.plot.emphasize = 'CoCoA';
% options.plot.styles = {':', '-', '--', '-.', '-', '--'};
options.label.Y = 'Solution cost';

%% Graph coloring experiment
exp(1) = load('C:\Develop\matlab\CoCoA\data\ijcai2016\exp_LocalInequalityConstraintCostFunction_delaunayGraph_i100_d3_n200_t20160127T165702_results.mat');
options.export.name = 'graph_coloring';

% Convert from cells to matrix
exp(1).settings.nMaxIterations = 70;
options.plot.range = 1:exp(1).settings.nMaxIterations;
resultsMat = prepareResults(exp(1).results, options.plot.range);
resultsMat = rmfield(resultsMat, 'MaxSumADVP');

% Create figure
options.figure.number = 187;
createResultGraph(resultsMat, exp(1).settings, 'costs', options);
options.figure.number = 188;
createResultGraph(resultsMat, exp(1).settings, 'msgs', options);
options.figure.number = 189;
createResultGraph(resultsMat, exp(1).settings, 'evals', options);

% Get results from best 1%
% analyzeResults(exp(1).results);

%% Game theory experiment
exp(2) = load('C:\Develop\matlab\CoCoA\data\ijcai2016\exp_LocalGameTheoreticCostFunction_delaunayGraph_i100_d3_n200_t20160126T151322_results.mat');
options.export.name = 'game_theory';

% Convert from cells to matrix
exp(2).settings.nMaxIterations = 100;
options.plot.range = 1:exp(2).settings.nMaxIterations;
resultsMat = prepareResults(exp(2).results, options.plot.range);

% Create figure
options.figure.number = 190;
createResultGraph(resultsMat, exp(2).settings, 'costs', options);
options.figure.number = 191;
createResultGraph(resultsMat, exp(2).settings, 'msgs', options);
options.figure.number = 192;
createResultGraph(resultsMat, exp(2).settings, 'evals', options);
% Get results from best 1%
% analyzeResults(exp(2).results);

%% Semi random experiments
exp(3) = load('C:\Develop\matlab\CoCoA\data\ijcai2016\exp_SemiRandomCostFunction_scalefreeGraph_i100_d10_n200_t20160125T075348_results.mat');
options.export.name = 'semi_random';

% Convert from cells to matrix
exp(3).settings.nMaxIterations = 800;
options.plot.range = 1:exp(3).settings.nMaxIterations;
resultsMat = prepareResults(exp(3).results, options.plot.range);

% Create figure
options.figure.number = 193;
createResultGraph(resultsMat, exp(3).settings, 'costs', options);
options.figure.number = 194;
createResultGraph(resultsMat, exp(3).settings, 'msgs', options);
options.figure.number = 195;
options.label.Y = 'Function evaluations';
createResultGraph(resultsMat, exp(3).settings, 'evals', options);
% Get results from best 1%
% analyzeResults(exp(3).results);

miniexp = load('C:\Develop\matlab\CoCoA\data\ijcai2016\exp_SemiRandomCostFunction_scalefreeGraph_i100_d10_n200_t20160129T161432_partial_results.mat');
% exp(3)
exp(3).results.MaxSumADVP = miniexp.results.MaxSumADVP;
%% Create result table

str = createResultTable(exp);
% clipboard('copy', str)
disp(str);

