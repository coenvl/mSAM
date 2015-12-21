%#ok<*SAGROW>
superclear

settings.nagents = 100;
settings.density = [.1 .2 .4 .5 .6 .7 .8 .9];
settings.numExps = 10;

solvers.DSA = 'nl.coenvl.sam.solvers.DSASolver';
solvers.CoCoA = 'nl.coenvl.sam.solvers.UniqueFirstCooperativeSolver';
solvers.Greedy = 'nl.coenvl.sam.solvers.GreedyLocalSolver';
solvers.MGM2 = 'nl.coenvl.sam.solvers.MGM2Solver';
% solvers.AFB = 'nl.coenvl.sam.solvers.FBSolver';
% solvers.CFL = 'nl.coenvl.sam.solvers.TickCFLSolver';

options.ncolors = uint16(4);
options.costFunction = 'nl.coenvl.sam.costfunctions.RandomCostFunction';
options.graphType = @randomGraph;

options.nStableIterations = uint16(25);
options.maxTime = 120;
options.waitTime = 15;
options.keepCostGraph = false;

solvertypes = fieldnames(solvers);

expname = sprintf('exp_density_%s_%d_%s_i%d', func2str(options.graphType), options.ncolors, datestr(now,30), settings.numExps);

for n = 1:numel(settings.density)
    options.graph.density = settings.density(n);
    options.graph.nAgents = uint16(settings.nagents);
    
    for e = 1:settings.numExps
        
        edges = feval(options.graphType, options.graph);
        fprintf('Graph generated with density %f (%s)\n', graphDensity(edges), func2str(options.graphType));
    
        for a = 1:numel(solvertypes)
            solvername = solvertypes{a};
            options.solverType = solvers.(solvername);
  
%             try
                fprintf('Performing experiment with %s(%d)\n', solvername, options.graph.nAgents);
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

%% Save results

save(fullfile('data', sprintf('%s_results.mat', expname)), 'settings', 'solvers', 'results');

% create_density_graphs;
