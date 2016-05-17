%#ok<*SAGROW>
superclear
warning('off', 'MATLAB:legend:PlotEmpty');
warning('off', 'MATLAB:legend:IgnoringExtraEntries');

%% Overall experiment settings
settings.numExps = 10; % i.e. number of problems generated
settings.nMaxIterations = [];
settings.nStableIterations = 200;
settings.nagents = 200;
settings.visualizeProgress = true;
settings.makeRandomConstraintCosts = false;

%% Create the experiment options
options.ncolors = uint16(10);
% options.constraint.type = 'nl.coenvl.sam.constraints.InequalityConstraint';
% options.constraint.arguments = {1};
% options.constraint.type = 'nl.coenvl.sam.constraints.CostMatrixConstraint';
% options.constraint.arguments = {[[1 0 3];[3 1 0];[0 3 1]], [[1 0 3];[3 1 0];[0 3 1]]};
options.constraint.type = 'nl.coenvl.sam.constraints.SemiRandomConstraint';
% options.constraint.type = 'nl.coenvl.sam.constraints.RandomConstraint';

options.graphType = @scalefreeGraph;
options.graph.maxLinks = uint16(4);
options.graph.initialsize = uint16(10);

% options.graphType = @randomGraph;
% options.graph.density = 0.1;

% options.graphType = @delaunayGraph;
% options.graph.sampleMethod = 'poisson';

% options.graphType = @nGridGraph;
% options.graph.nDims = uint16(3);
% options.graph.doWrap = '';

options.graph.nAgents = uint16(settings.nagents);

options.nStableIterations = uint16(settings.nStableIterations);
options.nMaxIterations = uint16(settings.nMaxIterations);
options.maxTime = 120;
options.waitTime = 1;
options.keepCostGraph = true;

solvers.ACLS = 'nl.coenvl.sam.solvers.ACLSSolver';
% solvers.ACLSUB = 'nl.coenvl.sam.solvers.ACLSUBSolver';
% solvers.ACLSProb = 'nl.coenvl.sam.solvers.ACLSProbSolver';
% solvers.AFB = 'nl.coenvl.sam.solvers.FBSolver';
% solvers.CFL = 'nl.coenvl.sam.solvers.TickCFLSolver';
solvers.CoCoA = 'nl.coenvl.sam.solvers.CoCoASolver';
solvers.CoCoS = 'nl.coenvl.sam.solvers.CoCoSolver';
solvers.ReCoCoS = 'nl.coenvl.sam.solvers.ReCoCoSolver';
% solvers.ReCoCoS2 = 'nl.coenvl.sam.solvers.ReCoCoSolverWorksGreat';
solvers.DSA = 'nl.coenvl.sam.solvers.DSASolver';
% solvers.Greedy = 'nl.coenvl.sam.solvers.GreedySolver';
% solvers.MaxSum = 'nl.coenvl.sam.solvers.MaxSumVariableSolver';
% solvers.MaxSumAD = 'nl.coenvl.sam.solvers.MaxSumADVariableSolver';
solvers.MaxSumADVP = 'nl.coenvl.sam.solvers.MaxSumADVPVariableSolver';
solvers.MCSMGM = 'nl.coenvl.sam.solvers.MCSMGMSolver';
% solvers.MGM = 'nl.coenvl.sam.solvers.MGMSolver';
solvers.MGM2 = 'nl.coenvl.sam.solvers.MGM2Solver';
% solvers.Random = 'nl.coenvl.sam.solvers.RandomSolver';
% solvers.SCA2 = 'nl.coenvl.sam.solvers.SCA2Solver';

%%
solvertypes = fieldnames(solvers);

C = strsplit(options.constraint.type, '.');
expname = sprintf('exp_%s_%s_i%d_d%d_n%d_t%s', C{end}, func2str(options.graphType), settings.numExps, options.ncolors, settings.nagents, datestr(now,30));

% Do the experiment
clear handles;
for e = 1:settings.numExps
    edges = feval(options.graphType, options.graph);
    
    if isfield(settings, 'makeRandomConstraintCosts') && settings.makeRandomConstraintCosts
        constraintCosts = randi(10, options.ncolors, options.ncolors, numel(edges));
        options.constraint.arguments = arrayfun(@(x) constraintCosts(:,:,x), 1:numel(edges), 'UniformOutput', false);
    else
        options.constraint.arguments = {};
    end
    
    for a = 1:numel(solvertypes)
        solvername = solvertypes{a};
        options.solverType = solvers.(solvername);

%         try
            fprintf('Performing experiment with %s (%d/%d)\n', solvername, e, settings.numExps);
            exp = doExperiment(edges, options);
            fprintf('Finished in t = %0.1f seconds\n', exp.time);
%         catch err
%             warning('Timeout or error occured:');
%             disp(err);
%             
%             exp.time = nan;
%             exp.allcost = nan;
%             exp.allevals = nan;
%             exp.allmsgs = nan;
%             exp.iterations = nan;
%             exp.alltimes = nan;
%         end
            
        results.(solvername).costs{e} = exp.allcost; 
        results.(solvername).evals{e} = exp.allevals;
        results.(solvername).msgs{e} = exp.allmsgs;
        results.(solvername).times{e} = exp.alltimes;
        results.(solvername).iterations(e) = exp.iterations;
        
        if settings.visualizeProgress
            visualizeProgress(exp, solvername);
        end
    end
end

%% Save results

save(fullfile('data', sprintf('%s_results.mat', expname)), 'settings', 'options', 'solvers', 'results');

%% Create graph

graphoptions = getGraphOptions();
graphoptions.figure.number = 188;
graphoptions.axes.yscale = 'linear'; % True for most situations
graphoptions.axes.xscale = 'linear';
graphoptions.axes.ymin = [];
% graphoptions.axes.xmax = 100;
graphoptions.export.do = false;
% graphoptions.export.name = expname;
graphoptions.label.Y = 'Solution Cost';
% graphoptions.label.X = 'Time';
graphoptions.plot.errorbar = false;
graphoptions.plot.emphasize = []; %'CoCoA';
% graphoptions.legend.location = 'NorthEast';
% graphoptions.legend.orientation = 'Horizontal';
% graphoptions.plot.x_fun = @(x) 1:max(x);
% graphoptions.plot.range = 1:1600;
resultsMat = prepareResults(results); %, graphoptions.plot.range);
createResultGraph(resultsMat, 'times', 'costs', graphoptions);
createResultTable(results)

