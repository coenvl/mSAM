superclear

settings.nagents = [50 100 150 200 250];
settings.numExps = 10;

solvers.DSA = 'nl.coenvl.sam.solvers.DSASolver';
solvers.CoCoA = 'nl.coenvl.sam.solvers.UniqueFirstCooperativeSolver';
solvers.Greedy = 'nl.coenvl.sam.solvers.GreedyLocalSolver';
solvers.MGM2 = 'nl.coenvl.sam.solvers.MGM2Solver';
% solvers.AFB = 'nl.coenvl.sam.solvers.FBSolver';
% solvers.CFL = 'nl.coenvl.sam.solvers.TickCFLSolver';

options.ncolors = uint16(5);
options.costFunction = 'nl.coenvl.sam.costfunctions.TaskSchedulingCostFunction';
options.graphType = @randomGraph;
%options.graphType = @distanceGraph3;

% options.graph.sampleMethod = 'poisson';
options.graph.sampleMethod = 'random';
% options.graph.density = 0.2;

options.keepCostGraph = false;
options.nStableIterations = uint16(25);
options.maxTime = 120;
options.waitTime = 15;

solvertypes = fieldnames(solvers);

expname = sprintf('exp_task_scheduling_%s_%d_%s_i%d', func2str(options.graphType), options.ncolors, datestr(now,30), settings.numExps);

for n = 1:numel(settings.nagents)
    options.graph.nAgents = uint16(settings.nagents(n));
    for e = 1:settings.numExps
        edges = feval(options.graphType, options.graph);
        
        fprintf('Graph generated with density %f (%s)\n', graphDensity(edges), func2str(options.graphType));
        
        options.agentProperties.data = randi(20, 1, graphSize(edges));
        options.agentProperties.ops = randi(20, 1, graphSize(edges));
        
        for a = 1:numel(solvertypes)
            solvername = solvertypes{a};
            options.solverType = solvers.(solvername);
            
            %         try
            fprintf('Performing experiment with %s(%d)\n', solvername, settings.nagents(n));
            exp = doExperiment(edges, options);
            %         catch err
            %             warning('Timeout or error occured:');
            %             disp(err);
            %             cost = nan; evals = nan; msgs = nan;
            %         end
            
            results.(solvername).costs(n,e) = exp.cost;
            results.(solvername).evals(n,e) = exp.evals;
            results.(solvername).msgs(n,e) = exp.msgs;
        end
    end
end

save(fullfile('data', sprintf('%s_results.mat', expname)), 'settings', 'solvers', 'results');

return;
%%

solverNames = fieldnames(solvers);
ytype = 'linear'; %'log'; % log or linear

costfig = figure(187); costax = gca; hold on; title('Costs'); ylabel('cost');
msgfig = figure(188); msgax = gca; hold on; title('Messages'); ylabel('msgs');
evalfig = figure(189); evalax = gca; hold on; title('Evals'); ylabel('evals');

set(costax, 'YScale', ytype);
set(msgax, 'YScale', ytype);
set(evalax, 'YScale', ytype);

colors = cubehelix(numel(solverNames) + 2, .5, -1.5, 3, 1);

n = 1;
for s = 1:numel(solverNames)
    solver = solverNames{s};
    bar(costax, n, mean(results.(solver).costs), 'FaceColor', colors(s+1,:));
    bar(msgax, n, mean(results.(solver).msgs), 'FaceColor', colors(s+1,:));
    bar(evalax, n, mean(results.(solver).evals), 'FaceColor', colors(s+1,:));
    n = n + 1;
end

set(costax, 'XTick', 1:numel(solverNames), 'XTickLabel', solverNames);
set(msgax, 'XTick', 1:numel(solverNames), 'XTickLabel', solverNames);
set(evalax, 'XTick', 1:numel(solverNames), 'XTickLabel', solverNames);