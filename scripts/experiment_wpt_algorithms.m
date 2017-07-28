%#ok<*SAGROW>
superclear
warning('off', 'MATLAB:legend:PlotEmpty');
warning('off', 'MATLAB:legend:IgnoringExtraEntries');

%% Overall experiment settings
settings.numExps = 200; % i.e. number of problems generated
settings.nMaxIterations = 0;
settings.nStableIterations = 40;
settings.nagents = 100;
settings.nreceivers = 75;
settings.nsensors = 50;
settings.ncolors = 20;
settings.graphType = 'customWPT';
settings.series = 'wpt';

%% Create the experiment options
options.ncolors = uint16(settings.ncolors);
options.graph.nAgents = uint16(settings.nagents);
options.nStableIterations = uint16(settings.nStableIterations);
options.nMaxIterations = uint16(settings.nMaxIterations);

solvers = getExperimentSolvers(settings.series);

%% Do the experiment
for e = 1:settings.numExps
    
    % Generate an experiment scenario
    [agentPos, receiverPos, sensorPos, edges, transmitter_to_receiver, transmitter_to_sensor] = generateWPTScenario(settings.nagents, settings.nreceivers, settings.nsensors);

    while (~higherOrderGraphIsConnected(edges, settings.nagents))
       [agentPos, receiverPos, sensorPos, edges, transmitter_to_receiver, transmitter_to_sensor] = generateWPTScenario(settings.nagents, settings.nreceivers, settings.nsensors);
    end   
    drawnow;
    
    % Set all positions
    for i = 1:numel(agentPos)
        options.agentProperties(i).position = agentPos{i};
    end
    options.constraint.arguments = receiverPos;
    options.sensorConstraint.arguments = sensorPos;

    exp = WPTExperiment(edges, options);
    
    for a = 1:numel(solvers)
        solvername = solvers(a).name;
        solverfield = matlab.lang.makeValidName(solvername);
        exp.initSolverType = solvers(a).initSolverType;
        exp.iterSolverType = solvers(a).iterSolverType; 
        
%         try
            fprintf('Performing experiment with %s (%d/%d)\n', solvername, e, settings.numExps);
            
            exp.run();
            fprintf('Finished in t = %0.1f seconds\n', exp.results.time(end));
            
            results.(solverfield).costs{e} = exp.results.cost; 
            results.(solverfield).evals{e} = exp.results.evals;
            results.(solverfield).msgs{e} = exp.results.msgs;
            results.(solverfield).times{e} = exp.results.time;
            results.(solverfield).iterations(e) = exp.results.numIters;
            visualizeProgress(exp, solvername);
%         catch err
%             warning('Timeout or error occured:');
%             disp(err);
%         end
    end
    
    % Calculate optimal value using Sinan's code
    EMR_Threshold = 0.018;
    tic;
    [x_lp, fval_lp] = LP_solution(transmitter_to_receiver,transmitter_to_sensor,EMR_Threshold,10,0);
    t_lp = toc;
    results.LPSolver.iterations(e) = 1;
    results.LPSolver.evals{e} = nan;
    results.LPSolver.msgs{e} = nan;
    results.LPSolver.costs{e} = fval_lp;
    results.LPSolver.times{e} = t_lp;
    
end

%% Save results
filename = sprintf('results_wpt_%s_i%d_d%d_n%d_t%s.mat', settings.graphType, settings.numExps, settings.ncolors, settings.nagents, datestr(now,30))
save(fullfile('data', filename), 'settings', 'options', 'solvers', 'results');

%% Create graph

graphoptions = getGraphOptions();
graphoptions.axes.yscale = 'linear'; % True for most situations
graphoptions.label.Y = 'Solution Cost';
graphoptions.label.X = 'Time';
graphoptions.plot.errorbar = false;
resultsMat = prepareResults(results);
createResultGraph(resultsMat, 'times', 'costs', graphoptions);

