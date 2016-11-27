%#ok<*SAGROW>
superclear
warning('off', 'MATLAB:legend:PlotEmpty');
warning('off', 'MATLAB:legend:IgnoringExtraEntries');

%% Overall experiment settings
settings.numExps = 100; % i.e. number of problems generated
settings.nStableIterations = 100;
settings.nMaxIterations = 0;
settings.nagents = 50;
settings.ncolors = 11;
settings.visualizeProgress = true;
settings.graphType = @scalefreeGraph;
settings.series = 'aaai17';

%% Create the experiment options
options.ncolors = uint16(settings.ncolors);

modelCosts = load('modelCosts.mat');
options.constraint.type = 'nl.coenvl.sam.constraints.CostMatrixConstraint';
options.constraint.arguments = {repmat(modelCosts.Qcost, numel(modelCosts.Qcost), 1)};

if isequal(settings.graphType, @scalefreeGraph)
    options.graphType = @scalefreeGraph;
    options.graph.maxLinks = uint16(4);
    options.graph.initialsize = uint16(10);
elseif isequal(settings.graphType, @randomGraph)
    options.graphType = @randomGraph;
    options.graph.density = 0.05;
elseif isequal(settings.graphType, @delaunayGraph)
    options.graphType = @delaunayGraph;
    options.graph.sampleMethod = 'poisson';
elseif isequal(settings.graphType, @nGridGraph)
    options.graphType = @nGridGraph;
    options.graph.nDims = uint16(3);
    options.graph.doWrap = '';
end

options.graph.nAgents = uint16(settings.nagents);
options.nStableIterations = uint16(settings.nStableIterations);
options.nMaxIterations = uint16(settings.nMaxIterations);

%% Solvers
solvers = getExperimentSolvers(settings.series);

%% Do the experiment
for e = 1:settings.numExps
    edges = feval(options.graphType, options.graph);

    exp = SensorNetworkExperiment(edges, options);
    
    for a = 1:numel(solvers)
        solvername = solvers(a).name;
        solverfield = matlab.lang.makeValidName(solvername);
        exp.initSolverType = solvers(a).initSolverType;
        exp.iterSolverType = solvers(a).iterSolverType; 
        
        try
            fprintf('Performing experiment with %s (%d/%d)\n', solvername, e, settings.numExps);

            exp.run();
            fprintf('Finished in t = %0.1f seconds\n\n', exp.results.time(end));
            
            results.(solverfield).costs{e} = exp.results.cost; 
            results.(solverfield).evals{e} = exp.results.evals;
            results.(solverfield).msgs{e} = exp.results.msgs;
            results.(solverfield).times{e} = exp.results.time;
            results.(solverfield).iterations(e) = exp.results.numIters;
            
            if settings.visualizeProgress
                visualizeProgress(exp, solverfield);
            end
        catch err
            warning('Timeout or error occured:');
            rethrow(err);
        end
    end
end

%% Save results
filename = sprintf('results_sensornet_%s_i%d_d%d_n%d_t%s.mat', func2str(settings.graphType), settings.numExps, settings.ncolors, settings.nagents, datestr(now,30))
save(fullfile('data', settings.series, filename), 'settings', 'options', 'solvers', 'results');

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
