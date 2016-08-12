%#ok<*SAGROW>
superclear
warning('off', 'MATLAB:legend:PlotEmpty');
warning('off', 'MATLAB:legend:IgnoringExtraEntries');

%% Overall experiment settings
settings.numExps = 20; % i.e. number of problems generated
settings.nMaxIterations = [];
settings.nStableIterations = 200;
settings.nagents = 50;
settings.nmeetings = 15;
settings.ncolors = 10; % i.e. timeslots
settings.visualizeProgress = true;

%% Create the experiment options
options.ncolors = uint16(settings.ncolors);

options.graphType = @meetingSchedulingGraph;
options.graph.nAgents = uint16(settings.nagents);
options.graph.nMeetings = uint16(settings.nmeetings);
% options.graph.meetingSizeFun = @(x) randi(x,varargin{:});
% options.graph.meetingSizeFun = @(x) 2 + ((x-2) ./ ((x+1) - randi(x)));
options.graph.meetingSizeFun = @(x) 2 + randi(x/10);

options.nStableIterations = uint16(settings.nStableIterations);
options.nMaxIterations = uint16(settings.nMaxIterations);
options.maxTime = 120;
options.waitTime = .01;
options.keepCostGraph = true;

solvers = {};

solvers(end+1).name = 'Greedy';
solvers(end).initSolverType = 'nl.coenvl.sam.solvers.GreedySolver';
solvers(end).iterSolverType = '';

solvers(end+1).name = 'CoCoA_simple';
solvers(end).initSolverType = 'nl.coenvl.sam.solvers.CoCoSolver';
solvers(end).iterSolverType = '';

solvers(end+1).name = 'CoCoA';
solvers(end).initSolverType = 'nl.coenvl.sam.solvers.CoCoASolver';
solvers(end).iterSolverType = '';

solvers(end+1).name = 'ACLS';
solvers(end).initSolverType = '';
solvers(end).iterSolverType = 'nl.coenvl.sam.solvers.ACLSSolver';

solvers(end+1).name = 'CoCoA - ACLS';
solvers(end).initSolverType = 'nl.coenvl.sam.solvers.CoCoASolver';
solvers(end).iterSolverType = 'nl.coenvl.sam.solvers.ACLSSolver';
% 
solvers(end+1).name = 'CoCoA - ACLSUB';
solvers(end).initSolverType = 'nl.coenvl.sam.solvers.CoCoASolver';
solvers(end).iterSolverType = 'nl.coenvl.sam.solvers.ACLSUBSolver';

solvers(end+1).name = 'DSA';
solvers(end).initSolverType = '';
solvers(end).iterSolverType = 'nl.coenvl.sam.solvers.DSASolver';

% solvers(end+1).name = 'CoCoA - DSA';
% solvers(end).initSolverType = 'nl.coenvl.sam.solvers.CoCoASolver';
% solvers(end).iterSolverType = 'nl.coenvl.sam.solvers.DSASolver';

solvers(end+1).name = 'MCSMGM';
solvers(end).initSolverType = '';
solvers(end).iterSolverType = 'nl.coenvl.sam.solvers.MCSMGMSolver';

solvers(end+1).name = 'CoCoA - MCSMGM';
solvers(end).initSolverType = 'nl.coenvl.sam.solvers.CoCoASolver';
solvers(end).iterSolverType = 'nl.coenvl.sam.solvers.MCSMGMSolver';

solvers(end+1).name = 'MGM2';
solvers(end).initSolverType = '';
solvers(end).iterSolverType = 'nl.coenvl.sam.solvers.MGM2Solver';

solvers(end+1).name = 'CoCoA - MGM2';
solvers(end).initSolverType = 'nl.coenvl.sam.solvers.CoCoASolver';
solvers(end).iterSolverType = 'nl.coenvl.sam.solvers.MGM2Solver';

solvers(end+1).name = 'Max-Sum';
solvers(end).initSolverType = '';
solvers(end).iterSolverType = 'nl.coenvl.sam.solvers.MaxSumVariableSolver';

% solvers(end+1).name = 'Max-Sum_ADVP';
% solvers(end).initSolverType = '';
% solvers(end).iterSolverType = 'nl.coenvl.sam.solvers.MaxSumADVPVariableSolver';

% solvers(end+1).name = 'CoCoA - Max-Sum_ADVP';
% solvers(end).initSolverType = 'nl.coenvl.sam.solvers.CoCoASolver';
% solvers(end).iterSolverType = 'nl.coenvl.sam.solvers.MaxSumADVPVariableSolver';


%%
expname = sprintf('exp_meetingScheduling_i%d_d%d_n%d_t%s', settings.numExps, options.ncolors, settings.nagents, datestr(now,30));

% Do the experiment
clear handles;
for e = 1:settings.numExps
    edges = feval(options.graphType, options.graph);
    
    while (~graphIsConnected(edges))
        warning('EXPERIMENT:GRAPHNOTCONNECTED', 'Graph not connected, retrying...');
        edges = feval(options.graphType, options.graph);
    end
    
    for j = 1:settings.nagents
        % Preference for certain timeslots
        options.agentProperties(j).preference = rand(1,settings.ncolors);
    end
    
    for a = 1:numel(solvers)
        solvername = solvers(a).name;
        options.initSolverType = solvers(a).initSolverType;
        options.iterSolverType = solvers(a).iterSolverType;

        try
            fprintf('Performing experiment with %s (%d/%d)\n', solvername, e, settings.numExps);
            exp = doMeetingSchedulingExperiment(edges, options);
            fprintf('Finished in t = %0.1f seconds\n', exp.time);
        catch err
            warning('Timeout or error occured:');
            disp(err);
            
            exp.time = nan;
            exp.allcost = nan;
            exp.allevals = nan;
            exp.allmsgs = nan;
            exp.iterations = nan;
            exp.alltimes = nan;
        end
        
        solverfield = matlab.lang.makeValidName(solvername);
        results.(solverfield).solver = solvers(a);
        results.(solverfield).costs{e} = exp.allcost; 
        results.(solverfield).evals{e} = exp.allevals;
        results.(solverfield).msgs{e} = exp.allmsgs;
        results.(solverfield).times{e} = exp.alltimes;
        results.(solverfield).iterations(e) = exp.iterations;
        
        if settings.visualizeProgress
            visualizeProgress(exp, solverfield);
        end
    end
end

%% Save results

save(fullfile('data', sprintf('%s_results.mat', expname)), 'settings', 'options', 'solvers', 'results');

%% Create graph

graphoptions = getGraphOptions();
graphoptions.figure.number = 188;
graphoptions.figure.height = 12;
graphoptions.axes.yscale = 'log'; % True for most situations
graphoptions.axes.xscale = 'linear';
graphoptions.axes.ymin = [];
graphoptions.axes.xmax = 30;
graphoptions.export.do = false;
% graphoptions.export.format = 'eps';
% graphoptions.export.name = expname;
graphoptions.label.Y = 'Solution Cost';
% graphoptions.label.X = 'Time';
graphoptions.plot.errorbar = false;
graphoptions.plot.emphasize = {}; %'CoCoA';
graphoptions.legend.box = 'off';
% graphoptions.legend.orientation = 'Horizontal';
graphoptions.plot.y_fun = @(x) nanmean(x,2);
graphoptions.plot.x_fun = @(x) nanmean(x,2);
% graphoptions.plot.x_fun = @(x) 1:max(x);
% graphoptions.plot.range = 1:1600;
resultsMat = prepareResults(results); %, graphoptions.plot.range);
createResultGraph(resultsMat, 'times', 'costs', graphoptions);
createResultTable(results)

