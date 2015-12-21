% superclear
clear all
clear java

rng(1, 'twister');

colornames = {'red', 'green', 'cyan', 'blue', 'yellow', 'magenta'};
options.ncolors = uint16(3);
% options.costFunction = 'nl.coenvl.sam.costfunctions.LocalInequalityConstraintCostFunction';
options.costFunction = 'nl.coenvl.sam.costfunctions.LocalGameTheoreticCostFunction';
% options.costFunction = 'nl.coenvl.sam.costfunctions.RandomCostFunction';

% options.solverType = 'nl.coenvl.sam.solvers.DSASolver';
% options.solverType = 'nl.coenvl.sam.solvers.UniqueFirstCooperativeSolver';
% options.solverType = 'nl.coenvl.sam.solvers.GreedyCooperativeSolver';
% options.solverType = 'nl.coenvl.sam.solvers.GreedyLocalSolver';
% options.solverType = 'nl.coenvl.sam.solvers.TickCFLSolver';
% options.solverType = 'nl.coenvl.sam.solvers.FBSolver';
% options.solverType = 'nl.coenvl.sam.solvers.MGMSolver';
% options.solverType = 'nl.coenvl.sam.solvers.SCA2Solver';
% options.solverType = 'nl.coenvl.sam.solvers.MGM2Solver';
% options.solverType = 'nl.coenvl.sam.solvers.ACLSSolver';
options.solverType = 'nl.coenvl.sam.solvers.MCSMGMSolver';

options.graph.nAgents = uint16(50);
options.graphType = @delaunayGraph;
options.graph.sampleMethod = 'poisson';
% options.graphType = @scalefreeGraph;
% options.graphType = @randomGraph;
% options.graph.density = .2;

options.nStableIterations = uint16(20);
% options.nMaxIterations = uint16(50);
options.keepCostGraph = true;

% Do the experiment
edges = feval(options.graphType, options.graph);
experimentResult = doExperiment(edges, options);

% plot(experimentResult.allcost)
% line([0 numel(experimentResult.allcost)], [experimentResult.cost experimentResult.cost]);
% shg

%% Show results
fprintf('\nExperiment results:\n');
fprintf('\tFound cost: %d\n', experimentResult.cost);
fprintf('\tFunction evals: %d\n', experimentResult.evals);
fprintf('\tNumber of msgs: %d\n', experimentResult.msgs);
fprintf('\tSolver type: %s\n', class(experimentResult.vars.solver(1)));
fprintf('\tCost function: %s\n', class(experimentResult.vars.costfun(1)));
fprintf('\tGraph type: %s\n', func2str(options.graphType));
fprintf('\t\t- size: %d\n', experimentResult.graph.nAgents);
fprintf('\t\t- density: %1.5f\n', experimentResult.graph.density);

figure(187)
if isfield(experimentResult, 'allcost')
    plot(experimentResult.allcost)
end

%% Create pretty graph (not required)

fid = fopen('test.gv', 'w');
fprintf(fid, 'graph {\n');
fprintf(fid, 'graph [K=0.03, overlap=false];\n');
fprintf(fid, 'node [style=filled];\n');

for i = 1:experimentResult.graph.nAgents
    var = experimentResult.vars.variable(i);
    if var.isSet()
        fprintf(fid, ' %d [fillcolor=%s, width=.2, height=.1];\n', i, colornames{double(var.getValue)});
    else
        fprintf(fid, ' %d [fillcolor=white, width=.2, height=.1];\n', i);
    end
end

for i = 1:size(experimentResult.graph.edges,1)
    fprintf(fid, ' %d -- %d;\n', experimentResult.graph.edges(i,1), experimentResult.graph.edges(i,2));
end
fprintf(fid, '}\n');
fclose(fid);

searchpath = 'c:\Progra~1';
files = dir(fullfile(searchpath, 'Graphviz*'));
if (numel(files) == 0)
    % Also search (x86)
    searchpath = 'c:\Progra~2';
    files = dir(fullfile(searchpath, 'Graphviz*'));
end

if (numel(files) == 0)
    error('Unable to complete printing graph, could not find Graphviz installation');
end
graphvizpath = fullfile(searchpath, files(1).name, 'bin');

layout = fullfile(graphvizpath, 'sfdp.exe');
% layout = fullfile(graphvizpath, 'twopi.exe');

gvfile = fullfile(pwd, 'test.gv');
pngfile = fullfile(pwd, 'test_colored.png');
[a,b] = system(sprintf('%s %s -T png -Gdpi=100 -o %s', layout, gvfile, pngfile));

img = imread(pngfile);
figure(188);
imshow(img);

fprintf('%s - Done!\n', datestr(now));