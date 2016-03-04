%#ok<*SAGROW>
superclear

%% Overall experiment settings
settings.numExps = 1; % i.e. number of problems generated
settings.nMaxIterations = 500;
settings.nStableIterations = [];
settings.nagents = 50;

%% Create the experiment options
options.ncolors = uint16(3);
options.costFunction = 'nl.coenvl.sam.costfunctions.LocalInequalityConstraintCostFunction';
% options.costFunction = 'nl.coenvl.sam.costfunctions.LocalGameTheoreticCostFunction';
% options.costFunction = 'nl.coenvl.sam.costfunctions.SemiRandomCostFunction';
% options.costFunction = 'nl.coenvl.sam.costfunctions.RandomCostFunction';
% options.costFunction = 'nl.coenvl.sam.costfunctions.NoisyLocalInequalityConstraintCostFunction';

% options.graphType = @scalefreeGraph;
% options.graph.maxLinks = uint16(4);
% options.graph.initialsize = uint16(10);

% options.graphType = @randomGraph;
% options.graph.density = 0.05;

options.graphType = @delaunayGraph;
options.graph.sampleMethod = 'poisson';
% options.graph.sampleMethod = 'random';

options.graph.nAgents = uint16(settings.nagents);

options.nStableIterations = uint16(settings.nStableIterations);
options.nMaxIterations = uint16(settings.nMaxIterations);
options.maxTime = 120;
options.waitTime = 1;
options.keepCostGraph = true;

% solvers.ACLS = 'nl.coenvl.sam.solvers.ACLSSolver';
% solvers.ACLSUB = 'nl.coenvl.sam.solvers.ACLSUBSolver';
% solvers.ACLSProb = 'nl.coenvl.sam.solvers.ACLSProbSolver';
% solvers.AFB = 'nl.coenvl.sam.solvers.FBSolver';
% solvers.CFL = 'nl.coenvl.sam.solvers.TickCFLSolver';
% solvers.CoCoA = 'nl.coenvl.sam.solvers.UniqueFirstCooperativeSolver';
solvers.DSA = 'nl.coenvl.sam.solvers.DSASolver';
% solvers.Greedy = 'nl.coenvl.sam.solvers.GreedyLocalSolver';
% solvers.MaxSum = 'nl.coenvl.sam.solvers.MaxSumVariableSolver';
% solvers.MaxSumAD = 'nl.coenvl.sam.solvers.MaxSumADVariableSolver';
% solvers.MaxSumADVP = 'nl.coenvl.sam.solvers.MaxSumADVPVariableSolver';
% solvers.MCSMGM = 'nl.coenvl.sam.solvers.MCSMGMSolver';
% solvers.MGM = 'nl.coenvl.sam.solvers.MGMSolver';
% solvers.MGM2 = 'nl.coenvl.sam.solvers.MGM2Solver';
solvers.Random = 'nl.coenvl.sam.solvers.RandomSolver';
% solvers.SCA2 = 'nl.coenvl.sam.solvers.SCA2Solver';

solvertypes = fieldnames(solvers);

C = strsplit(options.costFunction, '.');
expname = sprintf('exp_%s_%s_i%d_d%d_n%d_t%s', C{end}, func2str(options.graphType), settings.numExps, options.ncolors, settings.nagents, datestr(now,30));

% Do the experiment
clear handles;
for e = 1:settings.numExps
    edges = feval(options.graphType, options.graph);

    for a = 1:numel(solvertypes)
        solvername = solvertypes{a};
        options.solverType = solvers.(solvername);

        try
            fprintf('Performing experiment with %s (%d/%d)\n', solvername, e, settings.numExps);
            exp = doDynamicExperiment(edges, options);
%             exp = struct('allcost', rand(1,100), 'iterations', 100, 'allevals', rand(1,100), 'allmsgs', rand(1,100));
        catch err
            warning('Timeout or error occured:');
            disp(err);
            
            exp.allcost = nan;
            exp.allevals = nan;
            exp.allmsgs = nan;
            exp.iterations = nan;
        end
            
        results.(solvername).costs{e} = exp.allcost; 
        results.(solvername).evals{e} = exp.allevals;
        results.(solvername).msgs{e} = exp.allmsgs;
        results.(solvername).iterations{e} = exp.iterations;
        
        if ~exist('handles', 'var') || ~isfield(handles, 'fig') || ~ishandle(handles.fig)
            handles.fig = figure(007);
            handles.ax = gca(handles.fig);
            hold(handles.ax, 'on');
            handles.legend = legend(handles.ax, solvertypes);
        end
        
        if ~isfield(handles, solvername) || ~ishandle(handles.(solvername))
            handles.(solvername) = plot(exp.allcost, 'parent', handles.ax);
            handles.legend = legend(handles.ax, solvertypes);
        else
            set(handles.(solvername), 'XData', 1:numel(exp.allcost), 'YData', exp.allcost);
        end
        drawnow;
    end
end

%% Save results

save(fullfile('data', sprintf('%s_results.mat', expname)), 'settings', 'solvers', 'results');

%% Create graph

graphoptions = getGraphOptions();
graphoptions.axes.yscale = 'linear'; % True for most situations
graphoptions.axes.ymin = [];
% graphoptions.export.do = false;
% graphoptions.export.name = expname;
graphoptions.label.Y = 'Solution cost';
graphoptions.plot.errorbar = false;
graphoptions.plot.emphasize = []; %'CoCoA';
% graphoptions.legend.location = 'NorthEast';
% graphoptions.legend.orientation = 'Horizontal';
graphoptions.plot.range = 1:settings.nMaxIterations;
resultsMat = prepareResults(results, graphoptions.plot.range);
createResultGraph(resultsMat, settings, 'costs', graphoptions);

