%% Get Default options
graphOptions = getGraphOptions();
graphOptions.export.do = true;
graphOptions.export.dotables = true;
graphOptions.export.format = 'pdf';
graphOptions.export.folder = 'figures\aaai18\';
graphOptions.export.tables = 'figures\aaai18\tables';
graphOptions.plot.fixedStyles = getFixedAlgoStyles();
graphOptions.label.Y = 'Solution cost';
graphOptions.label.X = 'Running time (s)';

styles = {'-', '--', ':', '-.'};
colors = [0 0 0; .8 0 0; .4 .4 1; .9 .7 0];

%% Pseudo random coloring experiment
load('C:\Users\leeuwencjv\Develop\Projects\SAM\mSAM\data\hybrid\results_GraphColoringExperiment_SemiRandomConstraint_delaunayGraph_i200_d10_n200_t20170610T024736.mat')
graphOptions.export.name = 'asymmetric-speedup';

% Convert from cells to matrix
results = fixSleepyLaptop(results);
resultsMat = prepareResults(results);
createResultTable(results, graphOptions);

% initSolver.SSLA = initSolver.CoCoA_WPT;
initSolver = rmfield(initSolver, {'CoCoA', 'CoCoA_UF'});
iterSolver = rmfield(iterSolver, {'DSA', 'MGM2'});
legendEntries = {'Random', 'Greedy', 'SSLA'};

for iter = fieldnames(iterSolver)'
    fig = figure('name', sprintf('Results for %s', iter{:}));
    ax = gca(fig);
    hold(ax, 'on');
    
    s = 1;
    for init = fieldnames(initSolver)'
        solvername = sprintf('%s - %s', init{:}, iter{:});
        solverfield = matlab.lang.makeValidName(solvername);
        plot(mean(resultsMat.(solverfield).times, 2), ...
            mean(resultsMat.(solverfield).costs, 2), 'LineWidth', 3, ...
            'LineStyle', styles{s}, 'Color', colors(s,:));
        s = s + 1;
    end
    
    hl = legend(legendEntries);
    label_x = xlabel(ax, 'Time (s)');
    label_y = ylabel(ax, 'Solution cost');
 
    xlim([0 5]);
    ylim([1e4 3e4]);
    
    prettyExportFig(sprintf('%s_%s', graphOptions.export.name, iter{:}),...
        fig, ax, label_x, label_y, hl, graphOptions)
end

%% Graph coloring coloring experiment
load('C:\Users\leeuwencjv\Develop\Projects\SAM\mSAM\data\hybrid\results_GraphColoringExperiment_InequalityConstraint_delaunayGraph_i200_d3_n200_t20170804T224610.mat')
graphOptions.export.name = 'symmetric-speedup';

% Convert from cells to matrix
results = fixSleepyLaptop(results);
resultsMat = prepareResults(results);
createResultTable(results, graphOptions);

% initSolver.SSLA = initSolver.CoCoA_WPT;
legendEntries = {'Random', 'Greedy', 'SSLA'};

for iter = fieldnames(iterSolver)'
    fig = figure('name', sprintf('Results for %s', iter{:}));
    ax = gca(fig);
    hold(ax, 'on');
    
    s = 1;
    for init = fieldnames(initSolver)'
        solvername = sprintf('%s - %s', init{:}, iter{:});
        solverfield = matlab.lang.makeValidName(solvername);
        plot(mean(resultsMat.(solverfield).times, 2), ...
            mean(resultsMat.(solverfield).costs, 2), 'LineWidth', 3, ...
            'LineStyle', styles{s}, 'Color', colors(s,:));
        s = s + 1;
    end
    
    hl = legend(legendEntries);
    label_x = xlabel(ax, 'Time (s)');
    label_y = ylabel(ax, 'Solution cost');
 
    xlim([0 3]);
    ylim([40 200]);
    
    prettyExportFig(sprintf('%s_%s', graphOptions.export.name, iter{:}),...
        fig, ax, label_x, label_y, hl, graphOptions)
end

%%
load('C:\Users\leeuwencjv\Develop\Projects\SAM\mSAM\data\hybrid\results_repeatInit_InequalityConstraint_delaunayGraph_i200_d3_n100_t20170804T120652.mat')
graphOptions.export.name = 'repeated-exp';

% Convert from cells to matrix
results = fixSleepyLaptop(results);
resultsMat = prepareResults(results); %, graphoptions.plot.range);

for iter = fieldnames(iterSolver)'
    fig = figure('name', sprintf('Results for %s', iter{:}));
    clf(fig);
    ax = gca(fig);
    cla(ax);
    hold(ax, 'on');
    
    s = 1;
    for init = fieldnames(initSolver)'
        solvername = sprintf('%s - %s', init{:}, iter{:});
        solverfield = matlab.lang.makeValidName(solvername);
        plot(mean(resultsMat.(solverfield).times, 2), ...
            mean(resultsMat.(solverfield).costs, 2), 'LineWidth', 3, ...
            'LineStyle', styles{s}, 'Color', colors(s,:));
        s = s + 1;
    end

    s = 1;
    for init = fieldnames(initSolver)'
        solvername = sprintf('%s - %s', init{:}, iter{:});
        solverfield = matlab.lang.makeValidName(solvername);
        plot(mean(resultsMat.(solverfield).times, 2), ...
            min(resultsMat.(solverfield).costs, [], 2), 'LineWidth', 1, ...
            'LineStyle', styles{s}, 'Color', colors(s,:));
        plot(mean(resultsMat.(solverfield).times, 2), ...
            max(resultsMat.(solverfield).costs, [], 2), 'LineWidth', 1, ...
            'LineStyle', styles{s}, 'Color', colors(s,:));
        s = s + 1;
    end
    
    hl = legend('Random', 'SSLA');   
    label_x = xlabel(ax, 'Time (s)');
    label_y = ylabel(ax, 'Solution cost');
    
    xlim([0 2]);
    ylim([20 80]);
    
    prettyExportFig(sprintf('%s_%s', graphOptions.export.name, iter{:}),...
        fig, ax, label_x, label_y, hl, graphOptions);
end

%% Get the correlation
fig = figure('name', 'Correlation between init and end');
clf(fig);
ax = gca(fig);
randomStart = resultsMat.Random_DSA.costs(1,:);
randomEnd = resultsMat.Random_DSA.costs(end,:);
scatter(randomStart, randomEnd, 'filled', 'SizeData', 70);
label_x = xlabel('Final cost');
label_y = ylabel('Initial cost');
prettyExportFig('start-end-correlation',...
        fig, ax, label_x, label_y, [], graphOptions);
fprintf('Pearson correlation coefficient is %f', corr2(randomStart, randomEnd))

%% Export the bridge topology

edges = [1 2; 1 3; 1 4; 1 5; 1 6; 1 7; 1 8; 1 9; 1 10; 2 3; 3 4; 4 5; 5 6; 6 7; 7 8; 8 9; 9 10;
    11 12; 11 13; 11 14; 11 15; 11 16; 11 17; 11 18; 11 19; 11 20; 12 13; 13 14; 14 15; 15 16; 16 17; 17 18; 18 19; 19 20;
    1 11];

exp = GraphColoringExperiment(edges, struct);
exp.init();
exp.print(fullfile(graphOptions.export.folder, sprintf('bridge_topology.%s', graphOptions.export.format)));

%% Show the unfortunate pickings
load('C:\Users\leeuwencjv\Develop\Projects\SAM\mSAM\data\hybrid\results_GraphColoringExperiment_InequalityConstraint_manualGraph_i100_d3_n12_t20170707T173100.mat')

graphOptions.export.name = 'ManualGraph';
resultsMat = prepareResults(results); %, graphoptions.plot.range);

close all;

for iter = fieldnames(iterSolver)'
    fig = figure('name', sprintf('Results for %s', iter{:}));
    ax = gca(fig);
    hold(ax, 'on');
    
    s = 1;
    for init = fieldnames(initSolver)'
        solvername = sprintf('%s - %s', init{:}, iter{:});
        solverfield = matlab.lang.makeValidName(solvername);
        plot(mean(resultsMat.(solverfield).times, 2), ...
            mean(resultsMat.(solverfield).costs, 2), 'LineWidth', 3, ...
            'LineStyle', styles{s}, 'Color', colors(s,:));
        s = s + 1;
    end
    
%     hl = legend(legendEntries);
    hl = legend(fieldnames(initSolver));
    label_x = xlabel(ax, 'Time (s)');
    label_y = ylabel(ax, 'Solution cost');
 
    xlim([0 0.2]);
    ylim([0 8]);
    
    prettyExportFig(sprintf('%s_%s', graphOptions.export.name, iter{:}),...
        fig, ax, label_x, label_y, hl, graphOptions)
end