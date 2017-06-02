%#ok<*SAGROW>
superclear
warning('off', 'MATLAB:legend:PlotEmpty');
warning('off', 'MATLAB:legend:IgnoringExtraEntries');

%% Overall experiment settings
settings.numExps = 100; % i.e. number of problems generated
settings.nMaxIterations = 0;
settings.nStableIterations = 40;
settings.sizes = [4 8 16 32 64 128 192 256 384 512 768 1024];
settings.ncolors = 20;
settings.graphType = 'customWPT';
settings.series = 'wpt';

%% Create the experiment options
options.ncolors = uint16(settings.ncolors);
options.nStableIterations = uint16(settings.nStableIterations);
options.nMaxIterations = uint16(settings.nMaxIterations);

solvers = getExperimentSolvers(settings.series);
solvers = solvers(strcmp({solvers.name}, 'CoCoA_UF'));

%% Do the experiment
for e = 1:settings.numExps
    
    for a = 1:numel(settings.sizes)
        nagents = settings.sizes(a);
        nreceivers = ceil(nagents * .8);
        nsensors = ceil(nagents * .6);
        
        % Generate an experiment scenario
        [agentPos, receiverPos, sensorPos, edges, transmitter_to_receiver, transmitter_to_sensor] = generateWPTScenario(nagents, nreceivers, nsensors);

        while (~higherOrderGraphIsConnected(edges, nagents))
           [agentPos, receiverPos, sensorPos, edges, transmitter_to_receiver, transmitter_to_sensor] = generateWPTScenario(nagents, nreceivers, nsensors);
        end   
        drawnow;

        % Set all positions
        for i = 1:numel(agentPos)
            options.agentProperties(i).position = agentPos{i};
        end
        options.graph.nAgents = nagents;
        options.constraint.arguments = receiverPos;
        options.sensorConstraint.arguments = sensorPos;

        exp = WPTExperiment(edges, options);
    
        solvername = solvers.name;
        solverfield = matlab.lang.makeValidName(solvername);
        exp.initSolverType = solvers.initSolverType;
        exp.iterSolverType = solvers.iterSolverType; 
        
%         try
            fprintf('Performing experiment with %d agents (%d/%d)\n', nagents, e, settings.numExps);
            
            exp.run();
            fprintf('Finished in t = %0.1f seconds\n', exp.results.time(end));
            
            results.(solverfield).size(a) = nagents;
            results.(solverfield).costs(a, e) = exp.results.cost; 
            results.(solverfield).evals(a, e) = exp.results.evals;
            results.(solverfield).msgs(a, e) = exp.results.msgs;
            results.(solverfield).times(a, e) = exp.results.time;
            results.(solverfield).iterations(a, e) = exp.results.numIters;
%         catch err
%             warning('Timeout or error occured:');
%             disp(err);
%         end
    
        exp.reset();

        % Calculate optimal value using Sinan's code
        EMR_Threshold = 0.018;
        tic;
        [x_lp, fval_lp] = LP_solution(transmitter_to_receiver,transmitter_to_sensor,EMR_Threshold,10,0);
        t_lp = toc;
        
        results.LPSolver.size(a) = nagents;
        results.LPSolver.costs(a, e) = fval_lp;
        results.LPSolver.times(a, e) = t_lp;
    end
    
end

%% Save results
filename = sprintf('results_wpt_sizes_i%d_d%d_n%d_t%s.mat', settings.numExps, settings.ncolors, numel(settings.sizes), datestr(now,30))
save(fullfile('data', filename), 'settings', 'options', 'solvers', 'results');

%% Create graph

graphoptions = getGraphOptions();
graphoptions.axes.yscale = 'linear'; % True for most situations
graphoptions.label.Y = 'Solution Cost';
graphoptions.label.X = 'size';
graphoptions.plot.errorbar = false;
graphoptions.plot.styles = {'-'};
createResultGraph(results, 'size', 'costs', graphoptions);

