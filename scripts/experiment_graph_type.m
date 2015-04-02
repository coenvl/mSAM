%#ok<*SAGROW>
superclear

settings.nagents = 150;
settings.ncolors = 3;
settings.numExps = 1;

settings.graphType.Delaunay = @delaunayGraph;
settings.graphType.ScaleFree = @scalefreeGraph;
settings.graphType.Random = @randomGraph;
settings.graphType.CElegans = 'celegansneural';
settings.graphType.Surfnet = 'Surfnet';
settings.graphType.Electr208 = 'ElectrC_s208';
settings.graphType.Electr420 = 'ElectrC_s420';
settings.graphType.Electr838 = 'ElectrC_s838';
% settings.graphType.Power = 'Power';

solvers.Greedy = 'nl.coenvl.sam.solvers.GreedyLocalSolver';
% solvers.AFB = 'nl.coenvl.sam.solvers.FBSolver';
solvers.CFL = 'nl.coenvl.sam.solvers.TickCFLSolver';
solvers.DSA = 'nl.coenvl.sam.solvers.DSASolver';
solvers.CoCoA = 'nl.coenvl.sam.solvers.UniqueFirstCooperativeSolver';

%%
graphTypes = fieldnames(settings.graphType);
solvertypes = fieldnames(solvers);

expname = sprintf('exp_realgraphs_%d_%s_i%d', settings.ncolors, datestr(now,30), settings.numExps);

for n = 1:numel(graphTypes)
    for e = 1:settings.numExps
        
        graphName = graphTypes{n};
        graphFun = settings.graphType.(graphName);
        if ischar(graphFun)
            edges = realGraph(graphFun);
        else
            edges = feval(graphFun, settings.nagents);
            graphFun = func2str(graphFun);
        end
        
        numNodes = numel(unique(edges(:)));
        graphs.(graphName).size(e) = numNodes;
        graphs.(graphName).density(e) = (2*size(edges,1))/(numNodes*(numNodes-1));

        for a = 1:numel(solvertypes)
            solvername = solvertypes{a};
            solverclass = solvers.(solvername);
            
            try
                fprintf('Performing experiment with %s (%s)\n', solvername, graphName);
                [cost, evals, msgs] = doExperiment(numNodes, settings.ncolors, solverclass, edges);
            catch err
                warning('Timeout or error occured:');
                disp(err);
                cost = nan; evals = nan; msgs = nan;
            end
            results.(graphName).(solvername).costs(e) = cost; 
            results.(graphName).(solvername).evals(e) = evals;
            results.(graphName).(solvername).msgs(e) = msgs;
        end
    end
end

%% Save results

save(sprintf('%s_results.mat', expname), 'settings', 'solvers', 'results', 'graphs');

%% Create latex table;

graphNames = fieldnames(results);
solverNames = fieldnames(solvers);

nCols = 2 + (2 * numel(solverNames));
fprintf('\\begin{tabular}{|l|%s|}\n', repmat('|r',1,nCols));
fprintf('\\hline %14s & %s & %10s %s %s \\\\ \n', ...
    'Graph', 'Size', 'Density', ...
    sprintf('& %8s E ', solverNames{:}), ...
    sprintf('& %8s C ', solverNames{:}));
fprintf('\\hline\n'); 
for i = 1:numel(graphNames)
    graph = graphs.(graphNames{i});
    
    solverCostStr = '';
    solverEvalStr = '';
    for s = 1:numel(solverNames)
        cost = mean(results.(graphNames{i}).(solverNames{s}).costs);
        evals = mean(results.(graphNames{i}).(solverNames{s}).evals);
        solverCostStr = [solverCostStr sprintf('& % 10.1f ', cost)];
        solverEvalStr = [solverEvalStr sprintf('& % 10.1f ', evals)];
    end
    
    fprintf('\\hline %14s & %4d & % 10.04f %s %s \\\\ \n', ...
        graphNames{i}, mean(graph.size), mean(graph.density), ...
        solverEvalStr, solverCostStr);
end
fprintf('\\hline\n');
fprintf('\\end{tabular}\n');

%% Create graph



