%#ok<*SAGROW>
superclear
warning('off', 'MATLAB:legend:PlotEmpty');
warning('off', 'MATLAB:legend:IgnoringExtraEntries');

%% Overall experiment settings
settings.numExps = 50; % i.e. number of problems generated
settings.nMaxIterations = 0;
settings.nStableIterations = 100;
settings.nagents = 200;
settings.ncolors = 3;
settings.visualizeProgress = true;
settings.graphType = @randomGraph;
settings.series = 'hybrid';
settings.densities = [0.01 0.03 0.05 0.075 0.1 0.15 0.2 0.25 0.3];

%% Create the experiment options
options.ncolors = uint16(settings.ncolors);
options.constraint.type = 'nl.coenvl.sam.constraints.InequalityConstraint';

options.constraint.arguments = {};
options.agentProperties = struct;
options.maxTime = 300;
options.waitTime = .01;
options.useRootedSolvers = false;
options.debug = false;
options.ssetrack = false;

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
initSolver.Greedy = 'nl.coenvl.sam.solvers.GreedySolver';
% initSolver.CoCoA = 'nl.coenvl.sam.solvers.CoCoSolver';
% initSolver.CoCoA_UF = 'nl.coenvl.sam.solvers.CoCoASolver';
initSolver.CoCoA_WPT = 'nl.coenvl.sam.solvers.CoCoAWPTSolver';

clear iterSolver;
% iterSolver.NULL = '';
iterSolver.DSA = 'nl.coenvl.sam.solvers.DSASolver';
iterSolver.MGM2 = 'nl.coenvl.sam.solvers.MGM2Solver';
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
for e = 1:settings.numExps
    for d = 1:numel(settings.densities)
        options.graph.density = settings.densities(d);
        edges = feval(options.graphType, options.graph);

        while ~graphIsConnected(edges)
            edges = feval(options.graphType, options.graph);
        end

        exp = GraphColoringExperiment(edges, options);

        for a = 1:numel(solvers)
            solvername = solvers(a).name;
            solverfield = matlab.lang.makeValidName(solvername);
            exp.initSolverType = solvers(a).initSolverType;
            exp.iterSolverType = solvers(a).iterSolverType;

            if strfind(solvers(a).iterSolverType, 'MCSMGM')
                exp.nStableIterations = settings.nStableIterations * 5;
            else
                exp.nStableIterations = settings.nStableIterations;
            end

            %         try
            fprintf('Performing experiment with %s, d=%f (%d/%d)\n', solvername, options.graph.density, e, settings.numExps);
            exp.run();
            fprintf('\nFinished in t = %0.1f seconds\n', exp.results.time(end));

            results(d).(solverfield).costs{e} = exp.results.cost;
            results(d).(solverfield).evals{e} = exp.results.evals;
            results(d).(solverfield).msgs{e} = exp.results.msgs;
            results(d).(solverfield).times{e} = exp.results.time;
            results(d).(solverfield).iterations(e) = exp.results.numIters;
            results(d).(solverfield).density(e) = exp.graph.density;
            
            if settings.visualizeProgress
                visualizeProgress(exp, solverfield);
            end
            drawnow;
            pause(0.1);
            % return
            % catch err
            % warning('Timeout or error occured:');
            % disp(err);
            % end
            %         exp.reset();
        end
    end
end

%% Save results
saveResults

%% Create graph
close all;

density_times = struct;
density_costs = struct;

iterNames = fieldnames(iterSolver);
for d = 1:numel(settings.densities)
    resultsMat = prepareResults(results(d)); %, graphoptions.plot.range);
    for i = 1:numel(iterNames)
%         iter = {'DSA'};
        iter = iterNames{i};
        figure(186+i);
        subplot(3,3,d);
        cla;
        hold on;
        title(sprintf('%s - density %f', iter, settings.densities(d)));
        for init = fieldnames(initSolver)'
            solvername = sprintf('%s - %s', init{:}, iter);
            solverfield = matlab.lang.makeValidName(solvername);

            times = mean(resultsMat.(solverfield).times, 2);
            costs = mean(resultsMat.(solverfield).costs, 2);
            density_times.(solverfield)(d) = times(end);
            density_costs.(solverfield)(d) = costs(end);
            
            plot(times, costs, 'LineWidth', 3);
        end
%         h = legend(fieldnames(initSolver));
%         set(h,'interpreter', 'none'); 
    end
end

%% Timing improvements
figure(190);
cla
hold on
plot(settings.densities, density_times.Greedy_MGM2 ./ density_times.Random_MGM2);
plot(settings.densities, density_times.CoCoA_WPT_MGM2 ./ density_times.Random_MGM2);
xlabel('Graph density');
ylabel('Normalized time difference');
legend('Greedy', 'SSLA')

%% Cost improvements
figure(191);
cla
hold on
plot(settings.densities, density_costs.Random_MGM2 ./ density_costs.Random_MGM2);
plot(settings.densities, density_costs.Greedy_MGM2 ./ density_costs.Random_MGM2);
plot(settings.densities, density_costs.CoCoA_WPT_MGM2 ./ density_costs.Random_MGM2);
xlabel('Graph density');
ylabel('Normalized cost difference');
legend('Baseline', 'Greedy', 'SSLA')


%% one time fix
clear results2
for e = 1:4
    for d = 1:numel(settings.densities)
        for a = 1:numel(solvers)
            solvername = solvers(a).name;
            solverfield = matlab.lang.makeValidName(solvername);
            
            results2(d).(solverfield).costs{e} = results.(solverfield).costs{e, d};
            results2(d).(solverfield).evals{e} = results.(solverfield).evals{e, d};
            results2(d).(solverfield).msgs{e} = results.(solverfield).msgs{e, d};
            results2(d).(solverfield).times{e} = results.(solverfield).times{e, d};
            results2(d).(solverfield).iterations(e) = results.(solverfield).iterations(e, d);
            results2(d).(solverfield).density(e) = results.(solverfield).density(e, d);
        end
    end
end
