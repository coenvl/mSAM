%#ok<*SAGROW>
superclear

settings.nagents = [5 10 15 20 25]; % [50 100 150 200 250]; %
settings.numExps = 10;

options.ncolors = uint16(3);
% options.costFunction = 'nl.coenvl.sam.costfunctions.LocalInequalityConstraintCostFunction';
options.costFunction = 'nl.coenvl.sam.costfunctions.LocalGameTheoreticCostFunction';
% options.costFunction = 'nl.coenvl.sam.costfunctions.RandomCostFunction';

% options.graphType = @scalefreeGraph;
% options.graphType = @randomGraph;
options.graphType = @delaunayGraph;
options.graph.sampleMethod = 'poisson';
% options.graph.sampleMethod = 'random';
options.keepCostGraph = false;

options.nStableIterations = uint16(25);
options.maxTime = 120;
options.waitTime = 1;

solvers.DSA = 'nl.coenvl.sam.solvers.DSASolver';
solvers.CoCoA = 'nl.coenvl.sam.solvers.UniqueFirstCooperativeSolver';
solvers.Greedy = 'nl.coenvl.sam.solvers.GreedyLocalSolver';
% solvers.MGM = 'nl.coenvl.sam.solvers.MGMSolver';
solvers.MGM2 = 'nl.coenvl.sam.solvers.MGM2Solver';
% solvers.SCA2 = 'nl.coenvl.sam.solvers.SCA2Solver';
% solvers.AFB = 'nl.coenvl.sam.solvers.FBSolver';
% solvers.CFL = 'nl.coenvl.sam.solvers.TickCFLSolver';
solvers.ACLS = 'nl.coenvl.sam.solvers.ACLSSolver';
solvers.MCSMGM = 'nl.coenvl.sam.solvers.MCSMGMSolver';

solvertypes = fieldnames(solvers);

C = strsplit(options.costFunction, '.');
expname = sprintf('exp_%s_%s_i%d_c%d_t%s', C{end}, func2str(options.graphType), settings.numExps, options.ncolors, datestr(now,30));

for n = 1:numel(settings.nagents)
    options.graph.nAgents = uint16(settings.nagents(n));
    for e = 1:settings.numExps
        edges = feval(options.graphType, options.graph);

        for a = 1:numel(solvertypes)
            solvername = solvertypes{a};
            options.solverType = solvers.(solvername);
            
%             try
                fprintf('Performing experiment with %s(%d)\n', solvername, settings.nagents(n));
                exp = doExperiment(edges, options);
%             catch err
%                 warning('Timeout or error occured:');
%                 disp(err);
%                 cost = nan; evals = nan; msgs = nan;
%             end
            results.(solvername).costs(n,e) = exp.cost; 
            results.(solvername).evals(n,e) = exp.evals;
            results.(solvername).msgs(n,e) = exp.msgs;
        end
    end
end

%% Create graph

options = getGraphOptions();
options.axes.yscale = 'linear'; % True for most situations
options.export.do = false;
options.label.Y = 'Solution cost';
options.plot.hi_error_fun = @(x)x + 1;
options.plot.low_error_fun = @(x)x - 1;
createResultGraph(results, settings, 'costs', options);

%% Save results

save(fullfile('data', sprintf('%s_results.mat', expname)), 'settings', 'solvers', 'results');

% create_graphs;

