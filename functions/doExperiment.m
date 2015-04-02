function results = doExperiment(options)
%#ok<*AGROW>

%% Parse the options
nColors = getSubOption(uint16(3), 'uint16', options, 'ncolors');
nIterations = getSubOption(uint16(10), 'uint16', options, 'nIterations');
solverType = getSubOption('nl.coenvl.sam.solvers.UniqueFirstCooperativeSolver', ...
    'char', options, 'solverType');
costFunctionType = getSubOption('nl.coenvl.sam.costfunctions.LocalInequalityConstraintCostFunction', ...
    'char', options, 'costFunction');
maxtime = getSubOption(180, 'double', options, 'maxTime'); %maximum delay in seconds
waittime = getSubOption(1/2, 'double', options, 'waitTime'); %delay between checks
graphType = getSubOption(@randomGraph, 'function_handle', options, 'graphType');

%% Create problem
edges = feval(graphType, options.graph);
nagents = numel(unique(edges));

%% Setup the agents and variables
nl.coenvl.sam.ExperimentControl.ResetExperiment();

for i = 1:nagents
    varName = sprintf('variable%05d', i);
    agentName = sprintf('agent%05d', i);
    
    variable(i) = nl.coenvl.sam.variables.IntegerVariable(int32(1), int32(nColors), varName);
    if strcmp(solverType, 'nl.coenvl.sam.solvers.FBSolver')
        agent(i) = nl.coenvl.sam.agents.OrderedSolverAgent(agentName, variable(i));
        costfun(i) = feval(costFunctionType, agent(i).getSequenceID());
    else
        agent(i) = nl.coenvl.sam.agents.LocalSolverAgent(agentName, variable(i));
        costfun(i) = feval(costFunctionType, agent(i));
    end

    solver(i) = feval(solverType, agent(i), costfun(i));
    
    agent(i).setSolver(solver(i));
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
for i = 1:nagents
    k = find(edges(:,1) == i);
    if isempty(k)
        continue;
    end
    
    a = agent(i);
    c = costfun(i);
    for v = edges(k,2)'
        if strcmp(solverType, 'nl.coenvl.sam.solvers.FBSolver')
            c.addConstraintIndex(agent(v).getSequenceID());
            costfun(v).addConstraintIndex(a.getSequenceID());
        else
            a.addToNeighborhood(agent(v));
            agent(v).addToNeighborhood(a);
        end
    end
end

%% Init all agents
for i = nagents:-1:1
    agent(i).init();
    pause(.01);
end

%% Start the experiment

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
if isa(solver(1), 'nl.coenvl.sam.solvers.IterativeSolver')
    for i = 1:nIterations
        for j = 1:nagents
            solver(j).tick();
        end
        pause(.1);
    end
end

%% Wat for the algorithms to converge

% keyboard
for t = 1:(maxtime / waittime)
% while true
    pause(waittime);
    % This loop does not really work for algorithms that run iteratively
    if variable(randi(nagents)).isSet
%         fprintf('Experiment done...\n');
        break
    end
end

%% Get the results
if isa(costfun(1), 'nl.coenvl.sam.costfunctions.InequalityConstraintCostFunction')
    pc = nl.coenvl.sam.problemcontexts.IndexedProblemContext(-1);
    for i = 1:numel(variable)
        v = variable(i).getValue();
        if isa(v, 'double')
            pc.setValue(i-1, java.lang.Integer(v));
        else
            pc.setValue(i-1, v);
        end
    end
else
    pc = nl.coenvl.sam.problemcontexts.LocalProblemContext(agent(1));
    for i = 1:numel(variable)
        v = variable(i).getValue();
        if isa(v, 'double')
            pc.setValue(agent(i), java.lang.Integer(v));
        else
            pc.setValue(agent(i), v);
        end
    end
end

cost = 0;
for i = 1:numel(costfun)
%     cost = cost + costfun(i).currentValue();
    cost = cost + costfun(i).evaluate(pc);
end
cost = cost / 2 % Since symmetric cost functions
%%
% cost = 0;
% for i = 1:size(e,1)
%     varA = variable(e(i,1));
%     varB = variable(e(i,2));
%     if varA.getValue() == varB.getValue()
%         cost = cost + 1;
%     end
% end

results.vars.agent = agent;
results.vars.variable = variable;
results.vars.solver = solver;
results.vars.costfun = costfun;

results.cost = cost;
results.evals = nl.coenvl.sam.ExperimentControl.getNumberEvals();
results.msgs = nl.coenvl.sam.agents.AbstractSolverAgent.getTotalSentMessages;

results.graph.density = (2*size(edges,1))/(nagents*(nagents-1));
results.graph.edges = edges;
results.graph.nAgents = nagents;

end

