%#ok<*SAGROW>
superclear

settings.nagents = 100;
settings.density = [.05 .1 .15 .2];
settings.ncolors = 3;
settings.numExps = 1; %10;

settings.graphType = @randomGraph;

solvers.DSA = 'nl.coenvl.sam.solvers.DSASolver';
solvers.CoCoA = 'nl.coenvl.sam.solvers.UniqueFirstCooperativeSolver';
solvers.Greedy = 'nl.coenvl.sam.solvers.GreedyLocalSolver';
% solvers.AFB = 'nl.coenvl.sam.solvers.FBSolver';
% solvers.CFL = 'nl.coenvl.sam.solvers.TickCFLSolver';

solvertypes = fieldnames(solvers);

expname = sprintf('exp_density_%s_%d_%s_i%d', func2str(settings.graphType), settings.ncolors, datestr(now,30), settings.numExps);

for n = 1:numel(settings.density)
    for e = 1:settings.numExps
        edges = feval(settings.graphType, settings.nagents, settings.density(n));

        for a = 1:numel(solvertypes)
            solvername = solvertypes{a};
            solverclass = solvers.(solvername);
            
            try
                fprintf('Performing experiment with %s(%f)\n', solvername, settings.density(n));
                [cost, evals, msgs] = doExperiment(settings.nagents, settings.ncolors, solverclass, edges);
            catch err
                warning('Timeout or error occured:');
                disp(err);
                cost = nan; evals = nan; msgs = nan;
            end
            results.(solvername).costs(n,e) = cost; 
            results.(solvername).evals(n,e) = evals;
            results.(solvername).msgs(n,e) = msgs;
        end
    end
end

%% Save results

save(sprintf('data\%s_results.mat', expname), 'settings', 'solvers', 'results');

% create_graphs;
