function results = doExperiment(edges, options)
%#ok<*AGROW>

%% Parse the options
nColors = getSubOption(uint16(3), 'uint16', options, 'ncolors');
nStableIterations = getSubOption(uint16([]), 'uint16', options, 'nStableIterations');
nMaxIterations = getSubOption(uint16([]), 'uint16', options, 'nMaxIterations');
solverType = getSubOption('nl.coenvl.sam.solvers.UniqueFirstCooperativeSolver', ...
    'char', options, 'solverType');
costFunctionType = getSubOption('nl.coenvl.sam.costfunctions.LocalInequalityConstraintCostFunction', ...
    'char', options, 'costFunction');
maxtime = getSubOption(180, 'double', options, 'maxTime'); %maximum delay in seconds
waittime = getSubOption(1/2, 'double', options, 'waitTime'); %delay between checks
agentProps = getSubOption(struct, 'struct', options, 'agentProperties');
keepCostGraph = getSubOption(false, 'logical', options, 'keepCostGraph');

nagents = graphSize(edges);

%% Setup the agents and variables
nl.coenvl.sam.ExperimentControl.ResetExperiment();

if ~isempty(strfind(solverType, 'MaxSum'))
    nStableIterations = nStableIterations * 2.5;
end

fields = fieldnames(agentProps);
for i = 1:nagents
    varName = sprintf('variable%05d', i);
    agentName = sprintf('agent%05d', i);
    
    variable(i) = nl.coenvl.sam.variables.IntegerVariable(int32(1), int32(nColors), varName);
    if strcmp(solverType, 'nl.coenvl.sam.solvers.FBSolver')
        agent(i) = nl.coenvl.sam.agents.OrderedSolverAgent(agentName, variable(i));
        costfun(i) = nl.coenvl.sam.costfunctions.InequalityConstraintCostFunction(agent(i).getSequenceID());
        solver(i) = feval(solverType, agent(i), costfun(i));
    elseif ~isempty(strfind(solverType, 'MaxSum'))
        % Don't create the costfunction yet, that is for the constraints
        agent(i) = nl.coenvl.sam.agents.LocalSolverAgent(agentName, variable(i));
        solver(i) = feval(solverType, agent(i));
    else
        agent(i) = nl.coenvl.sam.agents.LocalSolverAgent(agentName, variable(i));
        costfun(i) = feval(costFunctionType, agent(i));
        solver(i) = feval(solverType, agent(i), costfun(i));
    end
    
    agent(i).setSolver(solver(i));
    
    for f = fields'
        prop = f{:};
        if numel(agentProps.(prop)) == 1
            agent(i).set(prop, agentProps.(prop));
        elseif numel(agentProps.(prop)) >= nagents
            agent(i).set(prop, agentProps.(prop)(i));
        else
            error('DOEXPERIMENT:INCORRECTPROPERTYCOUNT', ...
                'Incorrect number of properties, must be either 1 or number of agents (%d)', ...
                nagents);
        end
    end
    
    variable(i).clear();
    agent(i).reset();
end

%% Add children and parent if required
if isa(agent(1), 'nl.coenvl.sam.agents.OrderedAgent')
    % We assume a static final ordering where each agent has just one child
    % Therefore the algorithm is never really asynchronous
    for i = 1:(nagents-1)
        agent(i).addChild(agent(i+1));
    end
    
    for i = 2:nagents
        agent(i).setParent(agent(i-1));
    end
end

%% Add the constraints

if ~isempty(strfind(solverType, 'MaxSum'))
    for i = 1:size(edges,1)
        % Create constraint agent
        agentName = sprintf('constraint%05d', i);
        functionAgent(i) = nl.coenvl.sam.agents.LocalSolverAgent(agentName, null(1));
        costfun(i) = feval(costFunctionType, functionAgent(i));
        functionSolverType = char(solver(edges(i,1)).getCounterPart().getCanonicalName());
        functionsolver(i) = feval(functionSolverType, functionAgent(i), costfun(i));
        functionAgent(i).setSolver(functionsolver(i));
        functionAgent(i).reset();
        
        % Connect constraints to variables
        functionAgent(i).addToNeighborhood(agent(edges(i,1)));
        functionAgent(i).addToNeighborhood(agent(edges(i,2)));
        
        % And vice versa
        agent(edges(i,1)).addToNeighborhood(functionAgent(i));
        agent(edges(i,2)).addToNeighborhood(functionAgent(i));
        
        functionAgent(i).init();
    end
else
    for i = 1:nagents
        k = find(edges(:,1) == i);
        if isempty(k)
            continue;
        end

        a = agent(i);
        %s = solver(i);
        for v = edges(k,2)'
            if strcmp(solverType, 'nl.coenvl.sam.solvers.FBSolver')
                costfun(i).addConstraintIndex(agent(v).getSequenceID());
                costfun(v).addConstraintIndex(agent(i).getSequenceID());
                %s.addConstraint(agent(v));
                %solver(v).addConstraint(a);
            else
                % We are not MaxSum so create a "regular" graph
                a.addToNeighborhood(agent(v));
                agent(v).addToNeighborhood(a);
            end
        end
    end
end

%% Init all agents
for i = nagents:-1:1
    agent(i).init();
    pause(.01);
end

%% Start the experiment
t_experiment_start = tic; % start the clock here

startidx = randi(nagents);
a = solver(startidx);
if isa(a, 'nl.coenvl.sam.solvers.GreedyCooperativeSolver')
    msg = nl.coenvl.sam.messages.HashMessage('GreedyCooperativeSolver:PickAVar');
    a.push(msg);
elseif isa(a, 'nl.coenvl.sam.solvers.UniqueFirstCooperativeSolver')
    msg = nl.coenvl.sam.messages.HashMessage('UniqueFirstCooperativeSolver:PickAVar');
    a.push(msg);
elseif isa(a, 'nl.coenvl.sam.solvers.GreedyLocalSolver')
    msg = nl.coenvl.sam.messages.HashMessage('GreedyLocalSolver:AssignVariable');
    a.push(msg);
end

%% Do the iterations
numIters = 0;
if isa(solver(1), 'nl.coenvl.sam.solvers.IterativeSolver')
    %bestSolution = getCost(costfun, variable, agent);
    bestSolution = inf;
    
    if keepCostGraph; 
        costList = []; %bestSolution;
        evalList = nl.coenvl.sam.ExperimentControl.getNumberEvals();
        msgList = nl.coenvl.sam.agents.AbstractSolverAgent.getTotalSentMessages();
        timeList = toc(t_experiment_start);
    end
    
    % Iterate for AT LEAST nStableIterations
    countDown = nStableIterations;
    fprintf('Iteration:');
    while ~doStop(numIters, nMaxIterations, countDown, nStableIterations)       
        countDown = countDown - 1;
        numIters = numIters + 1;
        if (mod(numIters,25) == 0)
            fprintf(' %d', numIters);
        end
        for j = 1:nagents
            solver(j).tick();
        end

        if exist('functionsolver', 'var')
            for j = 1:numel(functionsolver)
                functionsolver(j).tick();
            end
        end
            
        cost = getCost(costFunctionType, variable, agent, edges);
        if keepCostGraph; 
            costList(numIters) = cost;
            evalList(numIters) = nl.coenvl.sam.ExperimentControl.getNumberEvals();
            msgList(numIters) = nl.coenvl.sam.agents.AbstractSolverAgent.getTotalSentMessages();
            timeList(numIters) = toc(t_experiment_start);
        end
        
        % If a better solution is found, reset countDown
        if cost < bestSolution
            countDown = nStableIterations;
            bestSolution = cost;
        end
    end
    fprintf(' ... done\n');
end
%% Wat for the algorithms to converge

% This loop does not really work for algorithms that run iteratively
for t = 1:(maxtime / waittime)
    pause(waittime);
    isset = arrayfun(@(x) x.isSet(), variable);

    if all(isset), break; end
end

%% Gather results to return
results.time = toc(t_experiment_start);
results.vars.agent = agent;
results.vars.variable = variable;
results.vars.solver = solver;
results.vars.costfun = costfun;

if exist('bestSolution', 'var')
    results.cost = bestSolution;
else
    results.cost = getCost(costFunctionType, variable, agent, edges);
end

if keepCostGraph && exist('costList', 'var')
    results.allcost = costList; 
else
    results.allcost = results.cost;
end

if keepCostGraph && exist('msgList', 'var')
    results.allmsgs = msgList; 
else
    results.allmsgs = nl.coenvl.sam.agents.AbstractSolverAgent.getTotalSentMessages();
end

if keepCostGraph && exist('evalList', 'var')
    results.allevals = evalList; 
else
    results.allevals = nl.coenvl.sam.ExperimentControl.getNumberEvals();
end

if keepCostGraph && exist('timeList', 'var')
    results.alltimes = timeList; 
else
    results.alltimes = results.time;
end

results.iterations = numIters;
results.evals = nl.coenvl.sam.ExperimentControl.getNumberEvals();
results.msgs = nl.coenvl.sam.agents.AbstractSolverAgent.getTotalSentMessages();

results.graph.density = graphDensity(edges);
results.graph.edges = edges;
results.graph.nAgents = nagents;

% clean up java objects
arrayfun(@(x) x.reset, agent);
nl.coenvl.sam.ExperimentControl.ResetExperiment();

end

function cost = getCost(costFunctionType, variable, agent, edges)

costfun = feval(costFunctionType, null(1));

cost = 0;
for i = 1:size(edges,1)
    a = edges(i,1);
    b = edges(i,2);
    
    pc = nl.coenvl.sam.problemcontexts.LocalProblemContext(agent(a));
       
%     if (variable(a).isSet())
        v = variable(a).getValue();
        if isa(v, 'double')
            pc.setValue(agent(a), java.lang.Integer(v));
        else
            pc.setValue(agent(a), v);
        end
%     end
    
%     if (variable(b).isSet())
        v = variable(b).getValue();
        if isa(v, 'double')
            pc.setValue(agent(b), java.lang.Integer(v));
        else
            pc.setValue(agent(b), v);
        end
%     end

    cost = cost + costfun(1).evaluateFull(pc);
end

end

function cost = getCost_old(costfun, variable, agent)

%% Get the results
if isa(costfun(1), 'nl.coenvl.sam.costfunctions.InequalityConstraintCostFunction')
    pc = nl.coenvl.sam.problemcontexts.IndexedProblemContext(-1);
    for i = 1:numel(variable)
        if (variable(i).isSet())
            v = variable(i).getValue();
            if isa(v, 'double')
                pc.setValue(i-1, java.lang.Integer(v));
            else
                pc.setValue(i-1, v);
            end
        end
    end
else
    pc = nl.coenvl.sam.problemcontexts.LocalProblemContext(agent(1));
    for i = 1:numel(variable)
        if (variable(i).isSet())
            v = variable(i).getValue();
            if isa(v, 'double')
                pc.setValue(agent(i), java.lang.Integer(v));
            else
                pc.setValue(agent(i), v);
            end
        end
    end
end

cost = 0;
for i = 1:numel(costfun)
%     cost = cost + costfun(i).currentValue();
    cost = cost + costfun(i).evaluate(pc);
end

%cost = cost / 2; % Since symmetric cost functions

end

% Stop as soon as one of the stopping criteria was met
function bool = doStop(numIters, nMaxIterations, countDown, nStableIterations)

bool = false;
if ~isempty(nMaxIterations) && (numIters >= nMaxIterations)
    bool = true;
end

if ~isempty(nStableIterations) && (countDown <= 0)
    bool = true;
end

return

end

