%#ok<*SAGROW>
superclear
warning('off', 'MATLAB:legend:PlotEmpty');
warning('off', 'MATLAB:legend:IgnoringExtraEntries');

%% Overall experiment settings
settings.numExps = 1; % i.e. number of problems generated
settings.nStableIterations = 100;
settings.nMaxIterations = 0;
settings.nagents = 4;
settings.ncolors = 200;
settings.visualizeProgress = true;
settings.graphType = 'customWPT';
settings.series = 'wpt';

%% Create the experiment options
options.ncolors = uint16(settings.ncolors);
options.graph.nAgents = uint16(settings.nagents);
options.nStableIterations = uint16(settings.nStableIterations);
options.nMaxIterations = uint16(settings.nMaxIterations);
options.initSolverType = '';
options.iterSolverType = 'nl.coenvl.sam.solvers.MaxSumVariableSolver';

% Transmitter positions
agentPos = {[1 2.3], [4 8], [6.5 4.5], [13 3]};
for i = 1:numel(agentPos)
    options.agentProperties(i).position = agentPos{i};
end

% Receiver positions
options.constraint.arguments = {[1.5 6], [5 2], [8 8], [10.3 5.2], [12.8 8]};

% Sensor positions
options.sensorConstraint.arguments = {[3.5 5], [9.2 2], [10.8 8]};

%% Do the experiment
%edges = {{1, [1 3], [2 3], [3 4], 4}, {[1 2 3], [3 4], 3}};
edges = {{1, [1 3], [2 3], [3 4], 4}, {[1 2], [3 4], 3}};

exp = WPTExperiment(edges, options);
    
exp.run();
fprintf('Finished in t = %0.1f seconds\n\n', exp.results.time(end));

return

%% Save results
saveResults

%% Create graph

graphoptions = getGraphOptions();
graphoptions.figure.number = 188;
graphoptions.axes.yscale = 'log';
graphoptions.axes.xscale = 'log';
% graphoptions.axes.xmax = 5;
graphoptions.export.do = false;
graphoptions.label.Y = 'Solution Cost';
graphoptions.label.X = 'Time (s)';
graphoptions.plot.errorbar = false;
% graphoptions.plot.emphasize = {'CoCoA'};
% graphoptions.legend.location = 'NorthEast';
% graphoptions.legend.orientation = 'Horizontal';
% graphoptions.plot.x_fun = @(x) 1:max(x);
% graphoptions.plot.range = 1:1600;
resultsMat = prepareResults(results); %, graphoptions.plot.range);
createResultGraph(resultsMat, 'times', 'costs', graphoptions);
