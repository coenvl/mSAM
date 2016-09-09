%% Get Default options
options = getGraphOptions();
options.plot.colors = cubehelix(8, .5, -1.5, 3, 1);
options.export.do = true;
options.export.format = 'eps';
options.export.folder = 'C:\Data\Documents\PHD\papers\aaai17\images';
options.plot.errorbar = false;
% options.plot.hi_error_fun = @(x) nanmean(x,2) + nanstd(x, [], 2);
% options.plot.lo_error_fun = @(x) nanmean(x,2) - nanstd(x, [], 2);
options.plot.y_fun = @(x) nanmean(x,2);
options.plot.x_fun = @(x) nanmean(x,2);
options.plot.emphasize = {'CoCoA_UF', 'CoCoA'};
% options.plot.styles = {':', '-', '--', '-.', '-', '--'};
options.label.Y = 'Solution cost';
options.label.X = 'Running time (s)';

%% Graph coloring experiment
exp(1) = load('data\aaai17\results_graphColoring_delaunayGraph_i100_d3_n200_t20160819T172541.mat');
% exp(1).results = fixSleepyLaptop(exp(1).results);
options.export.name = 'graph_coloring';

% Convert from cells to matrix
resultsMat = prepareResults(exp(1).results);

% Create figure
options.figure.number = 187;
options.axes.xmax = 3;
options.axes.ymin = 50;
createResultGraph(resultsMat, 'times', 'costs', options);
createResultTable(exp(1).results);
%% Semirandom experiment
exp(2) = load('data\aaai17\results_semirandom_scalefreeGraph_i100_d10_n200_t20160820T064440.mat');
exp(2).results = fixSleepyLaptop(exp(2).results);
options.export.name = 'semirandom';

% Convert from cells to matrix
resultsMat = prepareResults(exp(2).results);

% Create figure
options.figure.number = 190;
options.axes.xmax = 75;
options.axes.ymin = [];
createResultGraph(resultsMat, 'times', 'costs', options);
createResultTable(exp(2).results);
%% Meeting scheduling
exp(3) = load('data\aaai17\results_scheduling_i100_d20_n50_t20160826T095605.mat');
exp(3).results = fixSleepyLaptop(exp(3).results);
options.export.name = 'meetingScheduling';

% Convert from cells to matrix
resultsMat = prepareResults(exp(3).results);

% Create figure
options.figure.number = 193;
options.axes.ymax = 5e9;
options.axes.xmax = 3;
options.axes.yscale = 'log';
createResultGraph(resultsMat, 'times', 'costs', options);
createResultTable(exp(3).results);
%%
exp(4) = load('data\aaai17\results_sensornet_scalefreeGraph_i100_d11_n50_t20160829T142608');
exp(4).results = fixSleepyLaptop(exp(4).results);
options.export.name = 'sensornet';

% Convert from cells to matrix
resultsMat = prepareResults(exp(4).results);

options.figure.number = 194;
options.axes.ymax = [];
options.axes.xmin = .07;
options.axes.xmax = 10;
options.axes.yscale = 'log';
options.axes.xscale = 'log';
options.label.X = 'Running time (s)';
createResultGraph(resultsMat, 'times', 'costs', options);
createResultTable(exp(4).results);
%% Create result table

% str = createResultTable(exp);
% % clipboard('copy', str)
% disp(str);

