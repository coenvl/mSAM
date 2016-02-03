%#ok<*FNDSB>

%% Get Default options
options = getGraphOptions();
options.axes.yscale = 'linear'; % True for most situations
options.export.do = true;

%% Channel Allocation experiment
load data\exp_channel_alloc_distanceGraph3_11_20150717T093411_i10_results.mat
options.export.name = 'channel_alloc';

% Costs
options.figure.number = 187;
options.label.Y = 'Solution cost';
createResultGraph(results, settings, 'costs', options);

% Messages
options.figure.number = 188;
options.label.Y = 'Messages transmitted';
createResultGraph(results, settings, 'msgs', options);

% Evaluations
options.figure.number = 189;
options.label.Y = 'Costfun. evaluations';
createResultGraph(results, settings, 'evals', options);

%% Random Graph experiment
load data\exp_density_randomGraph_4_20150717T102628_i10_results.mat
options.export.name = 'density';

% Costs
options.figure.number = 190;
options.label.Y = 'Solution cost';
createResultGraph(results, settings, 'costs', options);

% Messages
options.figure.number = 191;
options.label.Y = 'Messages transmitted';
createResultGraph(results, settings, 'msgs', options);

% Evaluations
options.figure.number = 192;
options.label.Y = 'Costfun. evaluations';
createResultGraph(results, settings, 'evals', options);

%% Game Theory experiment
load data\exp_LocalGameTheoreticCostFunction_delaunayGraph_i10_c3_t20150717T142129_results.mat
options.export.name = 'gametheory';

% Costs
options.figure.number = 193;
options.label.Y = 'Solution cost';
createResultGraph(results, settings, 'costs', options);

% Messages
options.figure.number = 194;
options.label.Y = 'Messages transmitted';
createResultGraph(results, settings, 'msgs', options);

% Evaluations
options.figure.number = 195;
options.label.Y = 'Costfun. evaluations';
createResultGraph(results, settings, 'evals', options);

%% Graph Coloring experiment
load data\exp_LocalInequalityConstraintCostFunction_delaunayGraph_i10_c3_t20150717T141847_results.mat
options.export.name = 'coloring';

% Costs
options.figure.number = 196;
options.label.Y = 'Solution cost';
createResultGraph(results, settings, 'costs', options);

% Messages
options.figure.number = 197;
options.axes.yscale = 'log';
options.label.Y = 'Messages transmitted';
createResultGraph(results, settings, 'msgs', options);

% Evaluations
options.figure.number = 198;
options.label.Y = 'Costfun. evaluations';
createResultGraph(results, settings, 'evals', options);

% Reset yscale setting
options.axes.yscale = 'linear';

%% Task Scheduling experiment
load data\exp_task_scheduling_randomGraph_5_20150716T131749_i10_results.mat
options.export.name = 'task_scheduling';

% Costs
options.figure.number = 199;
options.label.Y = 'Solution cost';
createResultGraph(results, settings, 'costs', options);

% Messages
options.figure.number = 200;
options.label.Y = 'Messages transmitted';
createResultGraph(results, settings, 'msgs', options);

% Evaluations
options.figure.number = 201;
options.label.Y = 'Costfun. evaluations';
createResultGraph(results, settings, 'evals', options);