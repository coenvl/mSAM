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
settings.noise = linspace(0,0.2,9); % The stddev of the noise factor
settings.ncolors = 20;
settings.graphType = 'customWPT';
settings.series = 'wpt';

%% Create the experiment options
options.ncolors = uint16(settings.ncolors);
options.graph.nAgents = uint16(settings.nagents);
options.nStableIterations = uint16(settings.nStableIterations);
options.nMaxIterations = uint16(settings.nMaxIterations);

solvers = getExperimentSolvers(); %settings.series);
% solvers = solvers(strcmp({solvers.name}, 'CoCoA - ACLSUB'));
solvers = solvers(strcmp({solvers.name}, 'CoCoA_UF'));

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

    for a = 1:numel(settings.noise)       
        options.measurementNoise = settings.noise(a);
        
        % Generate the experiment
        exp = WPTExperiment(edges, options);
        
        solvername = solvers.name;
        solverfield = matlab.lang.makeValidName(solvername);
        exp.initSolverType = solvers.initSolverType;
        exp.iterSolverType = solvers.iterSolverType; 
        
%         try
            fprintf('Performing experiment with noise %0.03f (%d/%d)\n', options.measurementNoise, e, settings.numExps);
            
            exp.run();
            fprintf('Finished in t = %0.1f seconds\n', exp.results.time(end));
            
            results.(solverfield).noise(a) = options.measurementNoise;
            results.(solverfield).costs(a, e) = exp.results.cost(end); 
            results.(solverfield).evals(a, e) = exp.results.evals(end);
            results.(solverfield).msgs(a, e) = exp.results.msgs(end);
            results.(solverfield).times(a, e) = exp.results.time(end);
            results.(solverfield).iterations(a, e) = exp.results.numIters;
%         catch err
%             warning('Timeout or error occured:');
%             disp(err);
%         end
    
        

        % Calculate optimal value using Sinan's code
        EMR_Threshold = 0.018;
        tic;
        [x_lp, fval_lp] = LP_solution(transmitter_to_receiver,transmitter_to_sensor,EMR_Threshold,10,0);
        t_lp = toc;

        results.LPSolver.noise(a) = options.measurementNoise;
        results.LPSolver.costs(a, e) = fval_lp;
        results.LPSolver.times(a, e) = t_lp;
        
        % Also store the "noisy" cost
        domain = linspace(0, 10, settings.ncolors + 1);
        for i = 1:numel(exp.variable)
            x_lp_ind(i) = domain(find(x_lp(i) >= domain, 1, 'last'));
            assert(x_lp_ind(i) <= x_lp(i));
            exp.variable{i}.setValue(x_lp_ind(i));
        end
        results.LPSolverNoisy.noise(a) = options.measurementNoise;
        results.LPSolverNoisy.costs(a, e) = exp.getCost();
        
%         exp.reset();
    end 
end

%% Save results
filename = sprintf('results_wpt_noisy_i%d_d%d_n%d_t%s.mat', settings.numExps, settings.ncolors, numel(settings.noise), datestr(now,30))
save(fullfile('data', filename), 'settings', 'options', 'solvers', 'results');

%% Create graph

graphoptions = getGraphOptions();
graphoptions.axes.yscale = 'linear'; % True for most situations
graphoptions.label.Y = 'Solution Cost';
graphoptions.label.X = 'noise';
graphoptions.plot.errorbar = false;
graphoptions.plot.styles = {'-'};
createResultGraph(results, 'noise', 'costs', graphoptions);

