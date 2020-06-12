%#ok<*SAGROW>
superclear
warning('off', 'MATLAB:legend:PlotEmpty');
warning('off', 'MATLAB:legend:IgnoringExtraEntries');

%% Overall experiment settings
settings.numExps = 200; % i.e. number of times the same problem is repeated
settings.nMaxIterations = 0;
settings.nStableIterations = 100;
settings.nagents = 100;
settings.ncolors = 3;
settings.visualizeProgress = true;
settings.graphType = @delaunayGraph;
settings.series = 'hybrid';
settings.name = 'repeatInit';
settings.description = 'In this experiment the same problem is solver over and over with different initializations.';

%% Create the experiment options
options.ncolors = uint16(settings.ncolors);
options.constraint.type = 'nl.coenvl.sam.constraints.InequalityConstraint';

if isequal(settings.graphType, @scalefreeGraph)
    options.graphType = @scalefreeGraph;
    options.graph.maxLinks = uint16(4);
    options.graph.initialsize = uint16(10);
elseif isequal(settings.graphType, @randomGraph)
    options.graphType = @randomGraph;
    options.graph.density = 0.1;
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
initSolver.Random = 'nl.coenvl.sam.solvers.RandomSolver';
% initSolver.Greedy = 'nl.coenvl.sam.solvers.GreedySolver';
% initSolver.CoCoA = 'nl.coenvl.sam.solvers.CoCoSolver';
% initSolver.CoCoA_UF = 'nl.coenvl.sam.solvers.CoCoASolver';
initSolver.CoCoA_WPT = 'nl.coenvl.sam.solvers.CoCoAWPTSolver';

clear iterSolver;
% iterSolver.NULL = '';
% iterSolver.DSA = 'nl.coenvl.sam.solvers.DSASolver';
% iterSolver.MGM2 = 'nl.coenvl.sam.solvers.MGM2Solver';
iterSolver.ACLS = 'nl.coenvl.sam.solvers.ACLSSolver';
% iterSolver.ACLSUB = 'nl.coenvl.sam.solvers.ACLSUBSolver';
% iterSolver.MCSMGM = 'nl.coenvl.sam.solvers.MCSMGMSolver';

solvers = struct([]);
for init = fieldnames(initSolver)'
    for iter = fieldnames(iterSolver)'
        solvers(end + 1).name = sprintf('%s - %s', init{:}, iter{:});
        solvers(end).initSolverType = initSolver.(init{:});
        solvers(end).iterSolverType = iterSolver.(iter{:});
    end
end

%% Do the experiment
edges = feval(options.graphType, options.graph);

while ~graphIsConnected(edges)
    edges = feval(options.graphType, options.graph);
end

exp = GraphColoringExperiment(edges, options);
for e = 1:settings.numExps
    for a = 1:numel(solvers)
        solvername = solvers(a).name;
        solverfield = matlab.lang.makeValidName(solvername);
        exp.initSolverType = solvers(a).initSolverType;
        exp.iterSolverType = solvers(a).iterSolverType;

        %         try
        fprintf('Performing experiment with %s (%d/%d)\n', solvername, e, settings.numExps);
        exp.run();
        fprintf('\nFinished in t = %0.1f seconds\n', exp.results.time(end));
        
        results.(solverfield).costs{e} = exp.results.cost;
        results.(solverfield).evals{e} = exp.results.evals;
        results.(solverfield).msgs{e} = exp.results.msgs;
        results.(solverfield).times{e} = exp.results.time;
        results.(solverfield).iterations(e) = exp.results.numIters;
        results.(solverfield).density(e) = exp.graph.density;

        if settings.visualizeProgress
            visualizeProgress(exp, solverfield);
        end
        drawnow;
        pause(0.1);
    end
end

%% Save results
saveResults

%% Create graph
resultsMat = prepareResults(results); %, graphoptions.plot.range);
close all;
incoroporateUnsetComparison = true;
for iter = fieldnames(iterSolver)'
    figure();
    cla;
    hold on;
    title(iter);
    for init = fieldnames(initSolver)'
        solvername = sprintf('%s - %s', init{:}, iter{:});
        solverfield = matlab.lang.makeValidName(solvername);
        
        plot(mean(resultsMat.(solverfield).times, 2), mean(resultsMat.(solverfield).costs, 2), 'LineWidth', 3);
        %         density = mean(results.(solverfield).density);
        %         uniquevalexplored = mean([results.(solverfield).uniquevalexplored{:}]);
        %         numvalexplored = mean([results.(solverfield).allvalexplored{:}]);
        % %         allcombos = (density/2) * settings.nagents * settings.nagents * (settings.ncolors + incoroporateUnsetComparison) * (settings.ncolors + incoroporateUnsetComparison);
        % %         supercombos = (settings.ncolors + incoroporateUnsetComparison) ^ settings.nagents;
        %         fprintf('%s average values explored: %1.2f in %1.2f tries\n', ... %, (of %1.2f, so coverage %1.2f %%, precision %1.2f %%)\n', ...
        %             solvername, uniquevalexplored, numvalexplored); %, supercombos, 100 * uniquevalexplored / supercombos, 100 * uniquevalexplored / numvalexplored);
    end
    fprintf('\n');
    h = legend(fieldnames(initSolver));
    set(h,'interpreter', 'none');
    
    for init = fieldnames(initSolver)'
        solvername = sprintf('%s - %s', init{:}, iter{:});
        solverfield = matlab.lang.makeValidName(solvername);
        plot(mean(resultsMat.(solverfield).times, 2), min(resultsMat.(solverfield).costs, [], 2), 'LineWidth', 1);
        %         plot(mean(resultsMat.(solverfield).costs, 2), 'LineWidth', 3);
        plot(mean(resultsMat.(solverfield).times, 2), max(resultsMat.(solverfield).costs, [], 2), 'LineWidth', 1);
    end
end

% fprintf(


