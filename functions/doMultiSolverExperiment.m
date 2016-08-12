function results = doMultiSolverExperiment(edges, options)
%#ok<*AGROW>

%% Parse the options
try 
    nColors = getSubOption(uint16(3), 'uint16', options, 'ncolors');
    nStableIterations = getSubOption(uint16([]), 'uint16', options, 'nStableIterations');
    nMaxIterations = getSubOption(uint16([]), 'uint16', options, 'nMaxIterations');
    initSolverType = getSubOption('', 'char', options, 'initSolverType');
    iterSolverType = getSubOption('nl.coenvl.sam.solvers.CoCoASolver', ...
        'char', options, 'iterSolverType');
    constraintType = getSubOption('nl.coenvl.sam.constraints.InequalityConstraint', ...
        'char', options, 'constraint', 'type');
    constraintArgs = getSubOption({}, 'cell', options, 'constraint', 'arguments');
    maxtime = getSubOption(180, 'double', options, 'maxTime'); %maximum delay in seconds
    waittime = getSubOption(1/2, 'double', options, 'waitTime'); %delay between checks
    agentProps = getSubOption(struct, 'struct', options, 'agentProperties');
    keepCostGraph = getSubOption(false, 'logical', options, 'keepCostGraph');
catch err
    err.throwAsCaller();
end

nagents = graphSize(edges);

if strfind(iterSolverType, 'MCSMGM')
    nStableIterations = nStableIterations .* 2.5;
end

if isempty(initSolverType) && isempty(iterSolverType)
    error('Either set initSolverType, iterSolverType or both')
end

%% Setup the agents and variables
nl.coenvl.sam.ExperimentControl.ResetExperiment();

fields = fieldnames(agentProps);
for i = 1:nagents
    varName = sprintf('variable%05d', i);
    agentName = sprintf('agent%05d', i);
    
    variable(i) = nl.coenvl.sam.variables.IntegerVariable(int32(1), int32(nColors), varName);
%     if strfind(solverType, 'FBSolver')
%         agent(i) = nl.coenvl.sam.agents.LinkedAgent(variable(i), agentName);
%     else

%     if strfind(iterSolverType, 'MaxSum')
       agent(i) = nl.coenvl.sam.agents.VariableAgent(variable(i), agentName);
%     else
%        agent(i) = nl.coenvl.sam.agents.MultiSolverAgent(variable(i), agentName);
%     end
    
    if ~isempty(initSolverType)
        agent(i).setInitSolver(feval(initSolverType, agent(i)));
    end
    if ~isempty(iterSolverType)
        solver = feval(iterSolverType, agent(i));
        agent(i).setIterativeSolver(solver);
    end
    
    for f = fields'
        prop = f{:};
        if numel(agentPropts) == nagents
            agent(i).set(prop, agentProps(i).(prop));
        elseif numel(agentPropts) == 1 && numel(agentProps.(prop)) == 1
            agent(i).set(prop, agentProps.(prop));
        elseif numel(agentPropts) == 1 && numel(agentProps.(prop)) >= nagents
            agent(i).set(prop, agentProps.(prop)(i));
        else
            error('DOEXPERIMENT:INCORRECTPROPERTYCOUNT', ...
                'Incorrect number of properties, must be either 1 or number of agents (%d)', ...
                nagents);
        end
    end
    
    variable(i).clear();
end

%% Add the constraints

for i = 1:size(edges,1)
    a = edges(i,1);
    b = edges(i,2);
    
    if numel(constraintArgs) <= 1
        % Create a constraint always with same argument
        constraint(i) = feval(constraintType, variable(a), variable(b), constraintArgs{:});
    elseif numel(constraintArgs) == size(edges,1)
        % Create a constraint for every argument
        constraint(i) = feval(constraintType, variable(a), variable(b), constraintArgs{i});
    elseif mod(numel(constraintArgs), size(edges,1)) == 0
        error('Deprecated style of constraint argumentations!');
%         n = numel(constraintArgs) / size(edges,1);
%         k = (1+(i-1)*n):(i*n);
%         constraiint(i) = feval(constraintType, variable(a), variable(b), constraintArgs{k});
    else
        error('DOEXPERIMENT:INCORRECTARGUMENTCOUNT', ...
            'Incorrect number of constraint arguments, must be 0, 1 or number of edges (%d)', ...
            size(edges,1));
    end
    
    agent(a).addConstraint(constraint(i));
    agent(b).addConstraint(constraint(i));

    if ~isempty(strfind(iterSolverType, 'MaxSum'))
        % Create constraint agent
        agentName = sprintf('constraint%05d', i);
        constraintAgent(i) = nl.coenvl.sam.agents.ConstraintAgent(agentName, constraint(i), variable(a), variable(b));
        functionSolverType = char(solver.getCounterPart().getCanonicalName());
        constraintAgent(i).setSolver(feval(functionSolverType, constraintAgent(i)));
        
        % Set constraint agent address as targets
        agent(a).addFunctionAddress(constraintAgent(i).getID());
        agent(b).addFunctionAddress(constraintAgent(i).getID());
    end

end

%% Start the experiment
t_experiment_start = tic; % start the clock

% Special init for non-iterative solvers
a = agent(randi(nagents));
a.set(nl.coenvl.sam.solvers.CoCoSolver.ROOTNAME_PROPERTY, true);

arrayfun(@(x) x.init(), agent);
if exist('constraintAgent', 'var')
    arrayfun(@(x) x.init(), constraintAgent);
end

if ~isempty(initSolverType)
    pauseUntilVariablesAreSet(variable, maxtime, waittime)
end

%% Do the iterations
numIters = 0;

if ~isempty(iterSolverType)
    %bestSolution = getCost(costfun, variable, agent);
    bestSolution = inf;

    if keepCostGraph;
        costList = []; %bestSolution;
        evalList = nl.coenvl.sam.ExperimentControl.getNumberEvals();
        msgList = nl.coenvl.sam.MailMan.getTotalSentMessages();
        timeList = toc(t_experiment_start);
    end

    % Iterate for AT LEAST nStableIterations
    countDown = nStableIterations;
    fprintf('Iteration: ');
    while ~doStop(numIters, nMaxIterations, countDown, nStableIterations)
        countDown = countDown - 1;
        numIters = numIters + 1;
        if mod(numIters, 25) == 0
             fprintf(' %d', numIters);
        end

        arrayfun(@(x) x.tick, agent);

        if exist('constraintAgent', 'var')
            arrayfun(@(x) x.tick(), constraintAgent);
        end

        cost = getCost(constraint);
        if keepCostGraph;
            costList(numIters) = cost;
            evalList(numIters) = nl.coenvl.sam.ExperimentControl.getNumberEvals();
            msgList(numIters) = nl.coenvl.sam.MailMan.getTotalSentMessages();
            timeList(numIters) = toc(t_experiment_start);
        end

        % If a better solution is found, reset countDown
        if cost < bestSolution
            countDown = nStableIterations;
            bestSolution = cost;
        end
    end
    fprintf(' done\n')
end

% The solver is not iterative, but may take a while to complete
pauseUntilVariablesAreSet(variable, maxtime, waittime) 

%% Gather results to return
results.time = toc(t_experiment_start);
results.vars.agent = agent;
results.vars.variable = variable;
% results.vars.solver = solver;
results.vars.constraint = constraint;
if exist('constraintAgent', 'var')
    results.vars.constraintAgent = constraintAgent;
end

if exist('bestSolution', 'var')
    results.cost = bestSolution;
else
    results.cost = getCost(constraint);
end

if keepCostGraph && exist('costList', 'var')
    results.allcost = costList;
else
    results.allcost = results.cost;
end

if keepCostGraph && exist('msgList', 'var')
    results.allmsgs = msgList;
else
    results.allmsgs = nl.coenvl.sam.MailMan.getTotalSentMessages();
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

results.iterations = max(1,numIters);
results.evals = nl.coenvl.sam.ExperimentControl.getNumberEvals();
results.msgs = nl.coenvl.sam.MailMan.getTotalSentMessages();
%results.allmsgs = nl.coenvl.sam.MailMan.getSentMessages();

results.graph.density = graphDensity(edges);
results.graph.edges = edges;
results.graph.nAgents = nagents;

% clean up java objects
arrayfun(@(x) x.reset, agent);
nl.coenvl.sam.ExperimentControl.ResetExperiment();

end

function cost = getCost(constraint)
%% Get solution costs

cost = sum(arrayfun(@(x) x.getExternalCost(), constraint));

end

function pauseUntilVariablesAreSet(variable, maxtime, waittime)
% Wait until all variables are assigned a value

for t = 1:(maxtime / waittime)
    isset = arrayfun(@(x) x.isSet(), variable);

    if all(isset)
        return
    else
        pause(waittime);
    end
end

error('Solver did not terminate in %d seconds, increase maximum wait time if it is expected to terminate.', maxtime);

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

end

