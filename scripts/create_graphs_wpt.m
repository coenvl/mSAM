%% Get Default options
options = getGraphOptions();
% options.plot.colors = cubehelix(6, .5, -1.5, 3, 1);
options.export.do = true;
options.export.format = 'eps';
options.export.folder = fullfile(pwd, 'figures');
options.plot.errorbar = false;

%% Experiment with multiple algorithms
exp = load('data\results_wpt_customWPT_i200_d20_n100_t20170409T025034.mat');
exp.results = fixSleepyLaptop(exp(1).results);

% Convert from cells to matrix
resultsMat = prepareResults(exp.results);
% resultsMat = rmfield(resultsMat, 'Max_Sum');
resultsMat = renameStructField(resultsMat, 'CoCoA_UF', 'CoCoA_WPT');

% Create figure
options.figure.number = 187;
options.axes.xmax = 1.5;
options.export.name = 'algorithms';
options.label.Y = 'Solution Cost';
options.label.X = 'Time (s)';

createResultGraph(resultsMat, 'times', 'costs', options);
createResultTable(rmfield(exp.results, 'LPSolver'));

% Some additional results...
costs_max_sum = cellfun(@(x) x(end), exp.results.Max_Sum.costs);
costs_mcsmgm = cellfun(@(x) x(end), exp.results.MCSMGM.costs);
costs_cocoa = cellfun(@(x) x(end), exp.results.CoCoA_UF.costs);

fprintf('Max Sum violated constraints %0.03f%% of the time\n', 100*mean(costs_max_sum > 0));
fprintf('If it does not, the mean cost is %0.03f\n', mean(costs_max_sum(costs_max_sum < 0)));
fprintf('MCS-MGM only finds a good solution in %0.03f%% of the time\n', 100*mean(costs_mcsmgm > 0));
fprintf('And there CoCoA is %0.03f\n', mean(costs_cocoa(costs_max_sum < 0)));
%% Semirandom experiment
exp = load('data\results_wpt_sizes_i100_d20_n12_t20170414T220031.mat');
options.export.name = 'sizes';
exp.results = renameStructField(exp.results, 'CoCoA_UF', 'CoCoA_WPT');

% Create figure
options.figure.number = 188;
options.label.Y = 'Transmitted Power (W)';
options.label.X = '# Transmitters';
options.plot.styles = {':','-'};
options.axes.xmax = [];
createResultGraph(exp.results, 'size', 'costs', options);

%% Noisy experiment
exp = load('data\results_wpt_noisy_i200_d20_n9_t20170418T125941.mat');
options.export.name = 'noisy';

exp.results = rmfield(exp.results, 'LPSolverNoisy');
exp.results = renameStructField(exp.results, 'CoCoA_UF', 'CoCoA_WPT');

% Create figure
options.figure.number = 189;
options.label.Y = 'Transmitted Power (W)';
options.label.X = 'Model Noise (\sigma)';

createResultGraph(exp.results, 'noise', 'costs', options);

%% Plot the dynamic run experiment
exp = load('data\dynamic_withpower.mat')

fig = figure(190);
set(fig, 'Units', options.figure.units, ...
    'Position', [3 3 options.figure.width options.figure.height], ...
    'name', 'Transmitted power for dynamic simulation');

ax = cla;
hold(ax, 'on');

a = strcmp(exp.events, 'Adding agent');
r = strcmp(exp.events, 'Removing agent');

colors = cubehelix(2, .5, -1.5, 3, 1);

plot(exp.cost, ':', 'linewidth', options.plot.linewidth, 'color', colors(1,:));
plot(exp.lpcost, '-', 'linewidth', options.plot.linewidth, 'color', colors(2,:));

scatter(find(a), cost(a), 'g^', 'filled', 'sizedata', 200)
scatter(find(r), cost(r), 'rv', 'filled', 'sizedata', 200)

hl = legend(ax, 'TESSA', 'LP Solver');
set (hl, 'Location', options.legend.location, ...
    'fontsize', options.legend.fontsize, 'fontname', options.legend.font, 'linewidth', ...
    options.legend.linewidth, 'Box', options.legend.box, 'Interpreter', 'none');

set(ax, 'fontsize', options.axes.fontsize, 'fontname', options.axes.font, ...
    'linewidth', options.axes.linewidth, ...
    'YMinorGrid', options.axes.minorgrid, 'YMinorTick', options.axes.minortick, ...
    'XMinorGrid', options.axes.minorgrid, 'XMinorTick', options.axes.minortick, ...
    'Box', options.axes.box, 'YGrid', options.axes.grid, 'XGrid', options.axes.grid);

xlabel(ax, 'Time (s)', 'fontsize', options.label.fontsize, 'fontname', options.label.font);
ylabel(ax, 'Transmitted Power (W)', 'fontsize', options.label.fontsize, 'fontname', options.label.font);

xlim([0 360])
ylim([-0.25001 -0.09999])
if options.export.do
    filename = fullfile(options.export.folder, sprintf('dynamic.%s', options.export.format));
    export_fig(fig, filename, options.export.arguments{:}); 
end

%% Also plot all received powers on receivers
fig = figure(191);
set(fig, 'Units', options.figure.units, ...
    'Position', [3 3 options.figure.width options.figure.height], ...
    'name', 'Transmitted power for dynamic simulation');

ax = cla;
hold(ax, 'on');

colors = cubehelix(3, .5, -1.5, 3, 1);
colors(3,:) = [.7 .7 1];

power = exp.wpt(:,1:exp.settings.nreceivers);
      
plot(max(power,[],2), ':', 'linewidth', options.plot.linewidth, 'color', colors(1,:));
plot(mean(power,2), '-', 'linewidth', options.plot.linewidth, 'color', colors(2,:));
plot(min(power,[],2), ':', 'linewidth', options.plot.linewidth, 'color', colors(3,:));

set(ax, 'fontsize', options.axes.fontsize, 'fontname', options.axes.font, ...
    'linewidth', options.axes.linewidth, ...
    'YMinorGrid', options.axes.minorgrid, 'YMinorTick', options.axes.minortick, ...
    'XMinorGrid', options.axes.minorgrid, 'XMinorTick', options.axes.minortick, ...
    'Box', options.axes.box, 'YGrid', options.axes.grid, 'XGrid', options.axes.grid);

xlim([0 360])
xlabel(ax, 'Time (s)', 'fontsize', options.label.fontsize, 'fontname', options.label.font);
ylabel(ax, 'Received Power (W)', 'fontsize', options.label.fontsize, 'fontname', options.label.font);

if options.export.do
    filename = fullfile(options.export.folder, sprintf('dynamic_recv.%s', options.export.format));
    export_fig(fig, filename, options.export.arguments{:}); 
end
%%

fig = figure(192);
clf(fig);
set(fig, 'Units', options.figure.units, ...
    'Position', [3 3 options.figure.width options.figure.height], ...
    'name', 'Transmitted power for dynamic simulation');

ax = cla;
hold(ax, 'on');

colors = cubehelix(3, .5, -1.5, 3, 1);
colors(3,:) = [.7 .7 1];
power = exp.wpt(:,exp.settings.nreceivers+1:end);

plot(repmat(0.018,1,360), 'r-', 'linewidth', options.plot.linewidth);
plot(max(power,[],2), ':', 'linewidth', options.plot.linewidth, 'color', colors(1,:));
plot(mean(power,2), '-', 'linewidth', options.plot.linewidth, 'color', colors(2,:));
plot(min(power,[],2), ':', 'linewidth', options.plot.linewidth, 'color', colors(3,:));

set(ax, 'fontsize', options.axes.fontsize, 'fontname', options.axes.font, ...
    'linewidth', options.axes.linewidth, ...
    'YMinorGrid', options.axes.minorgrid, 'YMinorTick', options.axes.minortick, ...
    'XMinorGrid', options.axes.minorgrid, 'XMinorTick', options.axes.minortick, ...
    'Box', options.axes.box, 'YGrid', options.axes.grid, 'XGrid', options.axes.grid);

xlim([0 360])
ylim([0.01 0.02])
xlabel(ax, 'Time (s)', 'fontsize', options.label.fontsize, 'fontname', options.label.font);
ylabel(ax, 'EMR (W)', 'fontsize', options.label.fontsize, 'fontname', options.label.font);

if options.export.do
    filename = fullfile(options.export.folder, sprintf('dynamic_sens.%s', options.export.format));
    export_fig(fig, filename, options.export.arguments{:}); 
end

%% Show an example of the graph
[agentPos, receiverPos, sensorPos, edges, transmitter_to_receiver, transmitter_to_sensor] = generateWPTScenario(70,60,50);
while (~higherOrderGraphIsConnected(edges, 70))
   [agentPos, receiverPos, sensorPos, edges, transmitter_to_receiver, transmitter_to_sensor] = generateWPTScenario(70,60,50);
end  

hh = findobj;
t = arrayfun(@class, hh, 'UniformOutput', false);

% Prettify markers
hm = hh(strcmp(t, 'matlab.graphics.chart.primitive.Scatter'));
set(hm, 'SizeData', 100)

% Prettify legend
hl = hh(strcmp(t, 'matlab.graphics.illustration.Legend'));
set(hl, 'FontName', options.legend.font, 'FontSize', options.legend.fontsize, ...
    'Linewidth', options.legend.linewidth);

% Prettify axes
ax = hh(strcmp(t, 'matlab.graphics.axis.Axes'));
set(ax, 'linewidth', options.axes.linewidth, 'XTick',[], 'YTick', [], ...
    'YMinorGrid', 'off', 'YMinorTick', 'off', 'XMinorGrid', 'off', ...
    'XMinorTick', 'off', 'Box', 'on', ...
    'YGrid', 'off', 'XGrid', 'off');

% Make a square figure
fig = hh(strcmp(t, 'matlab.ui.Figure'));
set(fig, 'units', options.figure.units, 'Position', [3 3 20 20]);

% Zoom in a little bit
zf = 0.8;
xl = xlim;
yl = ylim;
xlim([1-zf zf] .* xl(2))
ylim([1-zf zf] .* yl(2))

filename = fullfile(options.export.folder, sprintf('graph.%s', options.export.format));
% export_fig(fig, filename, options.export.arguments{:}); 