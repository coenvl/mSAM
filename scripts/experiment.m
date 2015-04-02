% superclear
clear all
clear java

rng(1, 'twister');

colornames = {'red', 'green', 'blue', 'yellow', 'cyan', 'magenta'};
options.ncolors = uint16(3);

% options.solverType = 'nl.coenvl.sam.solvers.DSASolver';
% options.solverType = 'nl.coenvl.sam.solvers.UniqueFirstCooperativeSolver';
% options.solverType = 'nl.coenvl.sam.solvers.GreedyCooperativeSolver';
% options.solverType = 'nl.coenvl.sam.solvers.GreedyLocalSolver';
% options.solverType = 'nl.coenvl.sam.solvers.TickCFLSolver';
% options.solverType = 'nl.coenvl.sam.solvers.FBSolver';
% options.solverType = 'nl.coenvl.sam.solvers.MGMSolver';
options.solverType = 'nl.coenvl.sam.solvers.SCA2Solver';

options.graph.nAgents = uint16(40);
options.graphType = @delaunayGraph;
% options.graphType = @scalefreeGraph;
% options.graphType = @randomGraph;
% options.graph.density = .1;

options.nIterations = uint16(40);

experimentResult = doExperiment(options);

fprintf('\nExperiment results:\n');
fprintf('\tFound cost: %d\n', experimentResult.cost);
fprintf('\tFunction evals: %d\n', experimentResult.evals);
fprintf('\tNumber of msgs: %d\n', experimentResult.msgs);
fprintf('\tSolver type: %s\n', options.solverType);
fprintf('\tGraph type: %s\n', func2str(options.graphType));
fprintf('\t\t- size: %d\n', experimentResult.graph.nAgents);
fprintf('\t\t- density: %1.5f\n', experimentResult.graph.density);

return
%% The rest makes a nice graph, which is not required

fid = fopen('test.gv', 'w');
fprintf(fid, 'graph {\n');
fprintf(fid, 'graph [K=0.03, overlap=false];\n');
fprintf(fid, 'node [style=filled];\n');

if exist('variable', 'var')
    for i = 1:nagents
        var = variable(i);
        if var.isSet()
            fprintf(fid, ' %d [fillcolor=%s, width=.2, height=.1];\n', i, colornames{double(var.getValue)});
        else
            fprintf(fid, ' %d [fillcolor=white];\n', i);
        end
    end
end

for i = 1:size(edges,1)
    fprintf(fid, ' %d -- %d;\n', edges(i,1), edges(i,2));
end
fprintf(fid, '}\n');
fclose(fid);


graphvizpath = 'c:\Progra~2/Graphviz2.34/bin';

% layout = fullfile(graphvizpath, 'sfdp.exe');
layout = fullfile(graphvizpath, 'twopi.exe');

gvfile = fullfile(pwd, 'test.gv');
pngfile = fullfile(pwd, 'test.png');
[a,b] = system(sprintf('%s %s -T png -o %s', layout, gvfile, pngfile));

fprintf('%s - Done!\n', datestr(now));