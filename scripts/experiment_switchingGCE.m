%#ok<*SAGROW>
superclear
warning('off', 'MATLAB:legend:PlotEmpty');
warning('off', 'MATLAB:legend:IgnoringExtraEntries');

%% Overall experiment settings
settings.numExps = 100; % i.e. number of problems generated
settings.nMaxIterations = 0;
settings.nStableIterations = 100;
settings.nagents = 200;
settings.ncolors = 10;
settings.visualizeProgress = false;
settings.graphType = @randomGraph;
settings.series = 'hybrid+';

%% Create the experiment options
options.ncolors = uint16(settings.ncolors);
% options.constraint.type = 'nl.coenvl.sam.constraints.InequalityConstraint';
options.constraint.type = 'nl.coenvl.sam.constraints.SemiRandomConstraint';

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

% Take out max-sum, it is no use for coloring problems
k = arrayfun(@(x) isempty(strfind(x.name, 'Max-Sum')), solvers);
solvers = solvers(k);

% Take out the non-hybrids
k = arrayfun(@(x) ~isempty(x.iterSolverType), solvers);
solvers = solvers(k);

%% Do the experiment
for e = 1:settings.numExps
    edges = feval(options.graphType, options.graph);

    exp = SwitchingSolverGCE(edges, options);
    
    for a = 1:numel(solvers)
        solvername = solvers(a).name;
        solverfield = matlab.lang.makeValidName(solvername);
        exp.initSolverType = solvers(a).initSolverType;
        exp.iterSolverType = solvers(a).iterSolverType; 
        
%         try
            fprintf('Performing experiment with %s (%d/%d)\n', solvername, e, settings.numExps);
            
            exp.runI(100);
            exp.switchSolver();
            exp.run();
            fprintf('\nFinished in t = %0.1f seconds\n', exp.results.time(end));
            
            results.(solverfield).costs{e} = exp.results.cost; 
            results.(solverfield).evals{e} = exp.results.evals;
            results.(solverfield).msgs{e} = exp.results.msgs;
            results.(solverfield).times{e} = exp.results.time;
            results.(solverfield).iterations(e) = exp.results.numIters;
            
            if settings.visualizeProgress
                visualizeProgress(exp, solverfield);
            end
			drawnow;
			pause(0.1);
%         catch err
%             warning('Timeout or error occured:');
%             disp(err);
%         end
    end
end

%% Save results
saveResults

%% Create graph

graphoptions = getGraphOptions();
graphoptions.axes.yscale = 'linear'; % True for most situations
graphoptions.axes.ymin = .4;
graphoptions.axes.xmax = 10;
% graphoptions.export.do = false;
% graphoptions.export.name = expname;
graphoptions.label.Y = 'Solution Cost';
graphoptions.label.X = 'Time';
graphoptions.plot.errorbar = false;
% graphoptions.plot.emphasize = {'CoCoA'};
% graphoptions.legend.location = 'NorthEast';
% graphoptions.legend.orientation = 'Horizontal';
% graphoptions.plot.x_fun = @(x) 1:x;
graphoptions.plot.range = [];
resultsMat = prepareResults(results, graphoptions.plot.range);
createResultGraph(resultsMat, 'times', 'costs', graphoptions);

