%#ok<*SAGROW>
superclear
warning('off', 'MATLAB:legend:PlotEmpty');
warning('off', 'MATLAB:legend:IgnoringExtraEntries');

%% Overall experiment settings
settings.numExps = 1; % i.e. number of problems generated
settings.nStableIterations = 30;
settings.nMaxIterations = 30;
settings.nagents = 10;
settings.nreceivers = 8;
settings.nsensors = 6;
settings.ncolors = 20;
settings.visualizeProgress = true;
settings.graphType = 'customWPT';
settings.series = 'wpt';

%% Create the experiment options
options.ncolors = uint16(settings.ncolors);
options.graph.nAgents = uint16(settings.nagents);
options.nStableIterations = uint16(settings.nStableIterations);
options.nMaxIterations = uint16(settings.nMaxIterations);

if false
    options.initSolverType = '';
    options.iterSolverType = 'nl.coenvl.sam.solvers.MaxSumVariableSolver';
else
    options.initSolverType = 'nl.coenvl.sam.solvers.CoCoAWPTSolver';
    options.iterSolverType = '';%nl.coenvl.sam.solvers.DSASolver';
end

%% Build the scenario
while (~exist('tp', 'var') || ~higherOrderGraphIsConnected(edges, settings.nagents))
    % Ballpark something that might work
    MAX_ARITY = 4;
    
    tp = 100 * (rand(1,2) + poissonSample(settings.nagents));
    tp = tp - min(tp);
    for i = settings.nagents:-1:1
        agentPos{i} = tp(i,:);
    end

    rp = rand(1,2) + poissonSample(settings.nreceivers);
    rp = rp - min(rp);
    rp = rp .* (max(tp) ./ max(rp)); % Scale to same area as transmitters
    for i = settings.nreceivers:-1:1
        receiverPos{i} = rp(i,:);
    end

    sp = 100 * (rand(1,2) + poissonSample(settings.nsensors));
    sp = sp - min(sp);
    sp = sp .* (max(tp) ./ max(sp)); % Scale to same area as transmitters
    for i = settings.nsensors:-1:1
        sensorPos{i} = sp(i,:);
    end

    transmitter_to_receiver = pdist2(tp, rp, 'euclidean');
    transmitter_to_sensor = pdist2(tp, sp, 'euclidean');

    receiverEdges = getEdgesFromDistances(transmitter_to_receiver, MAX_ARITY);
    sensorEdges = getEdgesFromDistances(transmitter_to_sensor, MAX_ARITY);

    edges = {receiverEdges, sensorEdges};

    % receiver..
    t2r_mask = false(settings.nagents, settings.nreceivers);
    for a = 1:numel(receiverEdges)
        for r = receiverEdges{a}
            t2r_mask(r,a) = true;
        end
    end

    % sensor..
    t2s_mask = false(settings.nagents, settings.nsensors);
    for a = 1:numel(sensorEdges)
        for s = sensorEdges{a}
            t2s_mask(s,a) = true;
        end
    end

    transmitter_to_receiver(~t2r_mask) = inf;
    transmitter_to_sensor(~t2s_mask) = inf;

    fprintf('generated problem\n');
end   

%% Do the experiment
tic;
time = [];
cost = [];
events = {};
for t = 1:360
%     [agentPos, receiverPos, sensorPos, edges, transmitter_to_receiver, transmitter_to_sensor] = generateWPTScenario(settings.nagents, settings.nreceivers, settings.nsensors);
% 
%     while (~higherOrderGraphIsConnected(edges, settings.nagents))
%        [agentPos, receiverPos, sensorPos, edges, transmitter_to_receiver, transmitter_to_sensor] = generateWPTScenario(settings.nagents, settings.nreceivers, settings.nsensors);
%     end 
    
    for i = 1:numel(agentPos)
        options.agentProperties(i).position = agentPos{i};
    end

    % Receiver positions
    options.graph.nAgents = settings.nagents;
    options.constraint.arguments = receiverPos;
    options.sensorConstraint.arguments = sensorPos;

    exp = WPTExperiment(edges, options);

    exp.reset();
%     pause(.1);
    exp.run();
%     pause(.1);
    % fprintf('Finished in t = %0.1f seconds\n\n', exp.results.time(end));

    wpt(t,:) = cellfun(@(x) x.getPower(), exp.constraint);
    time(t) = toc;
    cost(t) = exp.getCost();
    
    [x_lp, fval_lp] = LP_solution(transmitter_to_receiver,transmitter_to_sensor,0.018,10,0);
    lpcost(t) = fval_lp;
    
    if (lpcost(t) > cost(t))
        keyboard;
    end
    
    plot([cost; lpcost]');
    drawnow;
%     pause(0.7);
    
    clear exp;
    options = rmfield(options, 'agentProperties');

    if (rand < .05)
        if (rand > .5 || settings.nagents <= MAX_ARITY)
            % Add one agent
            fprintf('ADDING AN AGENT!');
            events{t} = 'Adding agent';
            settings.nagents = settings.nagents + 1;
            tp(end + 1,:) = rand(1,2) .* max(tp);
            agentPos{end + 1} = tp(end,:);
        else
            % Remove one agent
            fprintf('REMOVING AN AGENT!');
            events{t} = 'Removing agent';
            k = ceil(rand * settings.nagents);
            tp(k,:) = [];
            settings.nagents = settings.nagents - 1;
            
            clear agentPos
            for i = settings.nagents:-1:1
                agentPos{i} = tp(i,:);
            end
        end
        
        transmitter_to_receiver = pdist2(tp, rp, 'euclidean');
        transmitter_to_sensor = pdist2(tp, sp, 'euclidean');
        
        receiverEdges = getEdgesFromDistances(transmitter_to_receiver, MAX_ARITY);
        sensorEdges = getEdgesFromDistances(transmitter_to_sensor, MAX_ARITY);
        
        edges = {receiverEdges, sensorEdges};

        t2r_mask = false(settings.nagents, settings.nreceivers);
        for a = 1:numel(receiverEdges)
            for r = receiverEdges{a}
                t2r_mask(r,a) = true;
            end
        end

        % sensor..
        t2s_mask = false(settings.nagents, settings.nsensors);
        for a = 1:numel(sensorEdges)
            for s = sensorEdges{a}
                t2s_mask(s,a) = true;
            end
        end

        transmitter_to_receiver(~t2r_mask) = inf;
        transmitter_to_sensor(~t2s_mask) = inf;
    end
end

%%
clf;
a = strcmp(events, 'Adding agent');
r = strcmp(events, 'Removing agent');
plot([cost; lpcost]');
hold on;
scatter(find(a), cost(a), 'g^')
scatter(find(r), cost(r), 'rv')

power = wpt(:,1:settings.nreceivers)
emr = wpt(:,settings.nreceivers+1:end)

% save('dynamic_1.mat')