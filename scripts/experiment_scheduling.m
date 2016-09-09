%#ok<*SAGROW>
superclear
warning('off', 'MATLAB:legend:PlotEmpty');
warning('off', 'MATLAB:legend:IgnoringExtraEntries');

%% Overall experiment settings
settings.numExps = 100; % i.e. number of problems generated
settings.nMaxIterations = 0;
settings.nStableIterations = 100;
settings.nagents = 50;
settings.nmeetings = 10;
settings.ncolors = 20; % i.e. timeslots
settings.visualizeProgress = true;
settings.series = 'aaai17';

%% Create the experiment options
options.ncolors = uint16(settings.ncolors);

options.graphType = @meetingSchedulingGraph;
options.graph.nMeetings = uint16(settings.nmeetings);
options.graph.meetingSizeFun = @(x) randi(10);

options.graph.nAgents = uint16(settings.nagents);
options.nStableIterations = uint16(settings.nStableIterations);
options.nMaxIterations = uint16(settings.nMaxIterations);

%% Solvers
solvers = {};

solvers(end+1).name = 'CoCoA';
solvers(end).initSolverType = 'nl.coenvl.sam.solvers.CoCoSolver';
solvers(end).iterSolverType = '';

solvers(end+1).name = 'CoCoA_UF';
solvers(end).initSolverType = 'nl.coenvl.sam.solvers.CoCoASolver';
solvers(end).iterSolverType = '';

solvers(end+1).name = 'ACLS';
solvers(end).initSolverType = '';
solvers(end).iterSolverType = 'nl.coenvl.sam.solvers.ACLSSolver';

% solvers(end+1).name = 'CoCoA - ACLS';
% solvers(end).initSolverType = 'nl.coenvl.sam.solvers.CoCoASolver';
% solvers(end).iterSolverType = 'nl.coenvl.sam.solvers.ACLSSolver';

% solvers(end+1).name = 'CoCoA - ACLSUB';
% solvers(end).initSolverType = 'nl.coenvl.sam.solvers.CoCoASolver';
% solvers(end).iterSolverType = 'nl.coenvl.sam.solvers.ACLSUBSolver';

solvers(end+1).name = 'DSA';
solvers(end).initSolverType = '';
solvers(end).iterSolverType = 'nl.coenvl.sam.solvers.DSASolver';

% solvers(end+1).name = 'CoCoA - DSA';
% solvers(end).initSolverType = 'nl.coenvl.sam.solvers.CoCoASolver';
% solvers(end).iterSolverType = 'nl.coenvl.sam.solvers.DSASolver';

solvers(end+1).name = 'MCSMGM';
solvers(end).initSolverType = '';
solvers(end).iterSolverType = 'nl.coenvl.sam.solvers.MCSMGMSolver';

% solvers(end+1).name = 'CoCoA - MCSMGM';
% solvers(end).initSolverType = 'nl.coenvl.sam.solvers.CoCoASolver';
% solvers(end).iterSolverType = 'nl.coenvl.sam.solvers.MCSMGMSolver';

solvers(end+1).name = 'MGM2';
solvers(end).initSolverType = '';
solvers(end).iterSolverType = 'nl.coenvl.sam.solvers.MGM2Solver';

% solvers(end+1).name = 'CoCoA - MGM2';
% solvers(end).initSolverType = 'nl.coenvl.sam.solvers.CoCoASolver';
% solvers(end).iterSolverType = 'nl.coenvl.sam.solvers.MGM2Solver';

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
for e = 1:settings.numExps
    edges = feval(options.graphType, options.graph);

    while (~graphIsConnected(edges))
        warning('EXPERIMENT:GRAPHNOTCONNECTED', 'Graph not connected, retrying...');
        edges = feval(options.graphType, options.graph);
    end
    
    for j = 1:settings.nagents
        % Add preference for certain timeslots
        options.agentProperties(j).preference = rand(1,settings.ncolors);
    end
    
    exp = SchedulingExperiment(edges, options); 
    
    for a = 1:numel(solvers)
        solvername = solvers(a).name;
        solverfield = matlab.lang.makeValidName(solvername);
        exp.initSolverType = solvers(a).initSolverType;
        exp.iterSolverType = solvers(a).iterSolverType; 

        try
            fprintf('Performing experiment with %s (%d/%d)\n', solvername, e, settings.numExps);

            exp.run();
            fprintf('Finished in t = %0.1f seconds\n', exp.results.time(end));
            
            results.(solverfield).costs{e} = exp.results.cost; 
            results.(solverfield).evals{e} = exp.results.evals;
            results.(solverfield).msgs{e} = exp.results.msgs;
            results.(solverfield).times{e} = exp.results.time;
            results.(solverfield).iterations(e) = exp.results.numIters;
            
            if settings.visualizeProgress
                visualizeProgress(exp, solverfield);
            end
        catch err
            warning('Timeout or error occured:');
            disp(err);
        end
    end
end

%% Save results
filename = sprintf('results_scheduling_i%d_d%d_n%d_t%s.mat', settings.numExps, settings.ncolors, settings.nagents, datestr(now,30))
save(fullfile('data', settings.series, filename), 'settings', 'options', 'solvers', 'results');

%% Create graph

graphoptions = getGraphOptions();
graphoptions.figure.number = 188;
graphoptions.figure.height = 12;
graphoptions.axes.yscale = 'log'; % True for most situations
graphoptions.axes.xscale = 'log';
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
graphoptions.plot.hi_error_fun = @(x) nanmean(x,2) + nanstd(x,[],2);
graphoptions.plot.lo_error_fun = @(x) nanmean(x,2) - nanstd(x,[],2);
resultsMat = prepareResults(results);
createResultGraph(resultsMat, 'times', 'costs', graphoptions);

