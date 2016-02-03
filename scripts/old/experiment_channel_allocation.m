superclear

settings.nagents = [50 100 150 200 250];
settings.numExps = 10;

solvers.DSA = 'nl.coenvl.sam.solvers.DSASolver';
solvers.CoCoA = 'nl.coenvl.sam.solvers.UniqueFirstCooperativeSolver';
solvers.Greedy = 'nl.coenvl.sam.solvers.GreedyLocalSolver';
solvers.MGM2 = 'nl.coenvl.sam.solvers.MGM2Solver';
% solvers.AFB = 'nl.coenvl.sam.solvers.FBSolver';
% solvers.CFL = 'nl.coenvl.sam.solvers.TickCFLSolver';

options.ncolors = uint16(11);
options.costFunction = 'nl.coenvl.sam.costfunctions.ChannelAllocationCostFunction';
options.graphType = @distanceGraph3;

% options.graph.sampleMethod = 'poisson';
options.graph.sampleMethod = 'random';
options.graph.maxDist = 30;
options.graph.scale = 50;

options.nStableIterations = uint16(25);
options.maxTime = 120;
options.waitTime = 15;
options.keepCostGraph = false;

solvertypes = fieldnames(solvers);

expname = sprintf('exp_channel_alloc_%s_%d_%s_i%d', func2str(options.graphType), options.ncolors, datestr(now,30), settings.numExps);

for n = 1:numel(settings.nagents)
    options.graph.nAgents = uint16(settings.nagents(n));
    for e = 1:settings.numExps
        [edges, pos] = feval(options.graphType, options.graph);
        
        fprintf('Graph generated with density %f (%s)\n', graphDensity(edges), func2str(options.graphType));
        
        options.agentProperties.xpos = pos(:,1);
        options.agentProperties.ypos = pos(:,2);
        if size(pos,2) > 2
            options.agentProperties.zpos = pos(:,3);
        end
        
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

%%

% solverNames = fieldnames(solvers);
% ytype = 'log'; % log or linear
% 
% costfig = figure(187); costax = gca; hold on; title('Costs'); ylabel('cost');
% msgfig = figure(188); msgax = gca; hold on; title('Messages'); ylabel('msgs');
% evalfig = figure(189); evalax = gca; hold on; title('Evals'); ylabel('evals');
% 
% set(costax, 'YScale', ytype);
% set(msgax, 'YScale', ytype);
% set(evalax, 'YScale', ytype);
% 
% colors = cubehelix(numel(solverNames) + 2, .5, -1.5, 3, 1);
% 
% n = 1;
% for s = 1:numel(solverNames)
%     solver = solverNames{s};
%     bar(costax, n, mean(results.(solver).costs), 'FaceColor', colors(s+1,:));
%     bar(msgax, n, mean(results.(solver).msgs), 'FaceColor', colors(s+1,:));
%     bar(evalax, n, mean(results.(solver).evals), 'FaceColor', colors(s+1,:));
%     n = n + 1;
% end
% 
% set(costax, 'XTick', 1:numel(solverNames), 'XTickLabel', solverNames);
% set(msgax, 'XTick', 1:numel(solverNames), 'XTickLabel', solverNames);
% set(evalax, 'XTick', 1:numel(solverNames), 'XTickLabel', solverNames);

return
%% Have a look at the interference functions

%           Select      Parabola     Limit
f = @(x,d) (x<3*d) .* ((x-3*d).^2 / (9*d.^2));

x = 0:0.1:90;
clf;
plot(x, f(x,5), 'r-');
hold on;
plot(x, f(x,10), 'b-');
plot(x, f(x,30), 'k-');

legend({'\Delta c = 2', '\Delta c = 1', '\Delta c = 0'})

title('Interference')
xlabel('Distance to neighbor')
ylabel('Interference ratio (cost)')