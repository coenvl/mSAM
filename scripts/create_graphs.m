%% Get Default options
options = getGraphOptions();
options.plot.colors = cubehelix(8, .5, -1.5, 3, 1);
options.export.do = true;
options.export.format = 'eps';
options.export.folder = 'C:\Data\Documents\PHD\ecai2016\images';
options.plot.errorbar = false;
% options.plot.hi_error_fun = @(x) nanmean(x,2) + nanstd(x, [], 2);
% options.plot.lo_error_fun = @(x) nanmean(x,2) - nanstd(x, [], 2);
options.plot.emphasize = {'CoCoA_UF', 'CoCoA'};
% options.plot.styles = {':', '-', '--', '-.', '-', '--'};
options.label.Y = 'Solution cost';

%% Graph coloring experiment
% exp(1) = load('C:\Develop\matlab\CoCoA\data\ijcai2016\exp_LocalInequalityConstraintCostFunction_delaunayGraph_i100_d3_n200_t20160127T165702_results.mat');
exp(1) = load('C:\Develop\matlab\CoCoA\data\exp_InequalityConstraint_delaunayGraph_i100_d3_n200_t20160410T173328_results.mat');
exp(1).results = fixSleepyLaptop(exp(1).results);
options.export.name = 'graph_coloring';

% Convert from cells to matrix
% exp(1).settings.nMaxIterations = 70;
% options.plot.range = 1:exp(1).settings.nMaxIterations;
resultsMat = prepareResults(exp(1).results);
resultsMat.CoCoA_UF = resultsMat.CoCoA;
resultsMat.CoCoA = resultsMat.CoCoS;
resultsMat = rmfield(resultsMat, 'CoCoS');
resultsMat = rmfield(resultsMat, 'MaxSumADVP');

% Create figure
options.figure.number = 187;
options.axes.xmax = 5;
options.label.X = 'Running time (s)';
createResultGraph(resultsMat, 'times', 'costs', options);
% options.figure.number = 188;
% createResultGraph(resultsMat, 'times', 'msgs', options);
% options.figure.number = 189;
% createResultGraph(resultsMat, 'times', 'evals', options);

% Get results from best 1%
% analyzeResults(exp(1).results);

%% Game theory experiment
% exp(2) = load('C:\Develop\matlab\CoCoA\data\ijcai2016\exp_LocalGameTheoreticCostFunction_delaunayGraph_i100_d3_n200_t20160126T151322_results.mat');
exp(2) = load('C:\Develop\matlab\CoCoA\data\exp_CostMatrixConstraint_delaunayGraph_i100_d3_n200_t20160411T080710_results.mat');;
exp(2).results = fixSleepyLaptop(exp(2).results);
options.export.name = 'game_theory';

% Convert from cells to matrix
% exp(2).settings.nMaxIterations = 100;
% options.plot.range = 1:exp(2).settings.nMaxIterations;
resultsMat = prepareResults(exp(2).results);
resultsMat.CoCoA_UF = resultsMat.CoCoA;
resultsMat.CoCoA = resultsMat.CoCoS;
resultsMat = rmfield(resultsMat, 'CoCoS');

% Create figure
options.figure.number = 190;
options.axes.xmax = 25;
options.label.X = 'Running time (s)';
createResultGraph(resultsMat, 'times', 'costs', options);
% options.figure.number = 191;
% createResultGraph(resultsMat, 'times', 'msgs', options);
% options.figure.number = 192;
% createResultGraph(resultsMat, 'times', 'evals', options);
% Get results from best 1%
% analyzeResults(exp(2).results);

%% Semi random experiments
% exp(3) = load('C:\Develop\matlab\CoCoA\data\ijcai2016\exp_SemiRandomCostFunction_scalefreeGraph_i100_d10_n200_t20160125T075348_results.mat');
exp(3) = load('C:\Develop\matlab\CoCoA\data\exp_SemiRandomConstraint_scalefreeGraph_i100_d10_n200_t20160413T123851_results.mat');
exp(3).results = fixSleepyLaptop(exp(3).results);
options.export.name = 'semi_random';

% Convert from cells to matrix
% exp(3).settings.nMaxIterations = 800;
% options.plot.range = 1:exp(3).settings.nMaxIterations;
resultsMat = prepareResults(exp(3).results);
resultsMat.CoCoA_UF = resultsMat.CoCoA;
resultsMat.CoCoA = resultsMat.CoCoS;
resultsMat = rmfield(resultsMat, 'CoCoS');

% Create figure
options.figure.number = 193;
options.axes.xmax = 180;
options.label.X = 'Running time (s)';
createResultGraph(resultsMat, 'times', 'costs', options);
% options.figure.number = 194;
% createResultGraph(resultsMat, 'times', 'msgs', options);
% options.figure.number = 195;
% options.label.Y = 'Function evaluations';
% createResultGraph(resultsMat, 'times', 'evals', options);
% Get results from best 1%
% analyzeResults(exp(3).results);

%%
exp(4) = load('C:\Develop\matlab\CoCoA\data\exp_CostMatrixConstraint_nGridGraph_i100_d5_n200_t20160415T153119_results.mat');
exp(4).results = fixSleepyLaptop(exp(4).results);
options.export.name = 'ngrid_random';

resultsMat = prepareResults(exp(4).results);
resultsMat.CoCoA_UF = resultsMat.CoCoA;
resultsMat.CoCoA = resultsMat.CoCoS;
resultsMat = rmfield(resultsMat, 'CoCoS');

options.figure.number = 194;
options.axes.xmax = 10;
options.axes.xscale = 'linear';
options.axes.ymin = 3500;
options.label.X = 'Running time (s)';
createResultGraph(resultsMat, 'times', 'costs', options);

%% Create result table

str = createResultTable(exp);
% clipboard('copy', str)
disp(str);

