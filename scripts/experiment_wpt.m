%#ok<*SAGROW>
superclear
warning('off', 'MATLAB:legend:PlotEmpty');
warning('off', 'MATLAB:legend:IgnoringExtraEntries');

%% Overall experiment settings
settings.numExps = 1; % i.e. number of problems generated
settings.nStableIterations = 30;
settings.nMaxIterations = 30;
settings.nagents = 150;
settings.nreceivers = 140;
settings.nsensors = 130;
settings.ncolors = 20;
settings.visualizeProgress = true;
settings.graphType = 'customWPT';
settings.series = 'wpt';

%% Create the experiment options
options.ncolors = uint16(settings.ncolors);
options.graph.nAgents = uint16(settings.nagents);
options.nStableIterations = uint16(settings.nStableIterations);
options.nMaxIterations = uint16(settings.nMaxIterations);

if false
    options.initSolverType = '';
    options.iterSolverType = 'nl.coenvl.sam.solvers.MaxSumVariableSolver';
else
    options.initSolverType = 'nl.coenvl.sam.solvers.CoCoAWPTSolver';
    options.iterSolverType = '';%nl.coenvl.sam.solvers.DSASolver';
end

%% Build the scenario
[agentPos, receiverPos, sensorPos, edges, transmitter_to_receiver, transmitter_to_sensor] = generateWPTScenario(settings.nagents, settings.nreceivers, settings.nsensors);

% Do not go for quad-constraints they are TOO much
% while (~isempty(findstr(options.iterSolverType, 'MaxSum')) && ...
%         max(cellfun(@(x) max(cellfun(@(y) numel(y), x)), edges)) > 3)
% while (max(cellfun(@(x) max(cellfun(@(y) numel(y), x)), edges)) > 3)
%    [agentPos, receiverPos, sensorPos, edges, transmitter_to_receiver, transmitter_to_sensor] = generateWPTScenario(settings.nagents, settings.nreceivers, settings.nsensors);
% end   

while (~higherOrderGraphIsConnected(edges, settings.nagents))
   [agentPos, receiverPos, sensorPos, edges, transmitter_to_receiver, transmitter_to_sensor] = generateWPTScenario(settings.nagents, settings.nreceivers, settings.nsensors);
end   

drawnow;
% waitforbuttonpress
%% Do the experiment
% Transmitter positions
for i = 1:numel(agentPos)
    options.agentProperties(i).position = agentPos{i};
end

% Receiver positions
options.constraint.arguments = receiverPos;

% Sensor positions
options.sensorConstraint.arguments = sensorPos;

exp = WPTExperiment(edges, options);

exp.reset();
exp.run();
fprintf('Finished in t = %0.1f seconds\n\n', exp.results.time(end));

results = exp.results;
% filename = sprintf('results_wpt_d%d_n%d_t%s.mat', settings.ncolors, settings.nagents, datestr(now,30))
% save(fullfile('data', filename), 'settings', 'options', 'results');

x_cocoa = cellfun(@(x) double(x.getValue()), exp.variable)
fval_cocoa = exp.getCost()
plot(exp.results.cost)
%% Now call Sinan's code
alpha = 1;
MAX_POWER = 10;
MIN_POWER = 0;
EMR_Threshold = 0.018;

% initial_powers = [10 10 10 10];

% centralized LP solutions
[x_lp, fval_lp] = LP_solution(transmitter_to_receiver,transmitter_to_sensor,EMR_Threshold,MAX_POWER,MIN_POWER)
% EMR_LP = transpose(path_loss_factor(transmitter_to_sensor))*x_lp;

return
%% validate

% f'.x : total transmitted power
f = -sum(path_loss_factor(transmitter_to_receiver),2);
f'*x_cocoa

% A.x <= b = EMR threshold 
A = transpose(path_loss_factor(transmitter_to_sensor));
A*x_cocoa
%% compute distance matrix
% A = cell2mat(agentPos');
% B = cell2mat(options.constraint.arguments');
% C = cell2mat(options.sensorConstraint.arguments');
% 
% t2r = sqrt(pdist2(A,B));
% t2s = sqrt(pdist2(A,C));
% 
% % receiver..
% t2r_mask = false(size(t2r));
% for a = 1:numel(edges{1})
%     for r = edges{1}{a}
%         t2r_mask(r,a) = true;
%     end
% end
% t2r(~t2r_mask) = inf
% 
% % sensor..
% t2s_mask = false(size(t2s));
% for a = 1:numel(edges{2})
%     for s = edges{2}{a}
%         t2s_mask(s,a) = true;
%     end
% end
% 
% t2s(~t2s_mask) = inf;

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
