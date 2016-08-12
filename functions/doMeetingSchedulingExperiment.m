function results = doMeetingSchedulingExperiment(edges, options)
%#ok<*AGROW>

%% Parse the options
try 
    nColors = getSubOption(uint16(3), 'uint16', options, 'ncolors');
    nStableIterations = getSubOption(uint16([]), 'uint16', options, 'nStableIterations');
    nMaxIterations = getSubOption(uint16([]), 'uint16', options, 'nMaxIterations');
    initSolverType = getSubOption('', 'char', options, 'initSolverType');
    iterSolverType = getSubOption('nl.coenvl.sam.solvers.CoCoASolver', ...
        'char', options, 'iterSolverType');
    maxtime = getSubOption(180, 'double', options, 'maxTime'); %maximum delay in seconds
    waittime = getSubOption(1/2, 'double', options, 'waitTime'); %delay between checks
    agentProps = getSubOption(struct, 'struct', options, 'agentProperties');
    keepCostGraph = getSubOption(false, 'logical', options, 'keepCostGraph');
catch err
    err.throwAsCaller();
end

pref_field = 'preference';
nagents = sum(cellfun(@(x) graphSize(x), edges));
nmeetings = numel(edges);

% if strfind(iterSolverType, 'MaxSum')
%     nStableIterations = nStableIterations .* 2.5;
% end

if isempty(initSolverType) && isempty(iterSolverType)
    error('Either set initSolverType, iterSolverType or both')
end

%% Setup the agents and variables
nl.coenvl.sam.ExperimentControl.ResetExperiment();

idx = 0;
for m = 1:nmeetings
    graph = edges{m};
    for i = unique(graph(:))'
        idx = idx + 1;
        varName = sprintf('variable%05d', idx);
        agentName = sprintf('agent%05d', idx);

        variable{m,i} = nl.coenvl.sam.variables.IntegerVariable(int32(1), int32(nColors), varName);
        newagent = nl.coenvl.sam.agents.VariableAgent(variable{m,i}, agentName);
        agent{m,i} = newagent;
        
        if ~isempty(initSolverType)
            newagent.setInitSolver(feval(initSolverType, newagent));
        end
        if ~isempty(iterSolverType)
            solver = feval(iterSolverType, newagent);
            newagent.setIterativeSolver(solver);
        end

        variable{m,i}.clear();
        
        % Set preference property
        newagent.set(pref_field, agentProps(i).(pref_field));
    end
end

%% Add the constraints

peqc = 'nl.coenvl.sam.constraints.PreferentialEqualityConstraint';
neqc = 'nl.coenvl.sam.constraints.InequalityConstraint';

constraintCost = 1e9;

idx = 0;
% First add constraints within meetings
for m = 1:nmeetings
    graph = edges{m};
    for i = 1:size(graph,1)
        a = graph(i,1);
        b = graph(i,2);
        idx = idx + 1;
        
        constraint{idx} = feval(peqc, variable{m,a}, variable{m,b}, agent{m,a}.get(pref_field), agent{m,b}.get(pref_field), constraintCost);

        agent{m,a}.addConstraint(constraint{idx});
        agent{m,b}.addConstraint(constraint{idx});

        if ~isempty(strfind(iterSolverType, 'MaxSum'))
            % Create constraint agent
            agentName = sprintf('constraint%05d', i);
            constraintAgent(idx) = nl.coenvl.sam.agents.ConstraintAgent(agentName, constraint{idx}, variable{m,a}, variable{m,b});
            functionSolverType = char(solver.getCounterPart().getCanonicalName());
            constraintAgent(idx).setSolver(feval(functionSolverType, constraintAgent(idx)));

            % Set constraint agent address as targets
            agent{m,a}.addFunctionAddress(constraintAgent(idx).getID());
            agent{m,b}.addFunctionAddress(constraintAgent(idx).getID());
        end
    end
end

% Next add constraints between (same) agents
for a = 1:size(agent,2)
    for m1 = 1:nmeetings
        if ~isempty(agent{m1,a})
            for m2 = (m1 + 1):nmeetings
                if ~isempty(agent{m2,a})
                    % add constraints between agent {m1,a} and {m2,a}
                    idx = idx + 1;
        
                    constraint{idx} = feval(neqc, variable{m1,a}, variable{m2,a}, constraintCost);

                    agent{m1,a}.addConstraint(constraint{idx});
                    agent{m2,a}.addConstraint(constraint{idx});

                    if ~isempty(strfind(iterSolverType, 'MaxSum'))
                        % Create constraint agent
                        agentName = sprintf('constraint%05d', i);
                        constraintAgent(idx) = nl.coenvl.sam.agents.ConstraintAgent(agentName, constraint{idx}, variable{m1,a}, variable{m2,a});
                        functionSolverType = char(solver.getCounterPart().getCanonicalName());
                        constraintAgent(idx).setSolver(feval(functionSolverType, constraintAgent(idx)));

                        % Set constraint agent address as targets
                        agent{m1,a}.addFunctionAddress(constraintAgent(idx).getID());
                        agent{m2,a}.addFunctionAddress(constraintAgent(idx).getID());
                    end
                end
            end
        end
    end 
end


allagents = [agent{:}];
allvariables = [variable{:}];

%% Start the experiment
t_experiment_start = tic; % start the clock

% Special init for non-iterative solvers
a = allagents(randi(numel(allagents)));
a.set(nl.coenvl.sam.solvers.CoCoSolver.ROOTNAME_PROPERTY, true);

arrayfun(@(x) x.init(), allagents);
if exist('constraintAgent', 'var')
    arrayfun(@(x) x.init(), constraintAgent);
end

if ~isempty(initSolverType)
    pauseUntilVariablesAreSet(allvariables, maxtime, waittime)
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

        arrayfun(@(x) x.tick, allagents);

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
pauseUntilVariablesAreSet(allvariables, maxtime, waittime) 

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

% results.graph.density = graphDensity(edges);
results.graph.edges = edges;
results.graph.nAgents = nagents;

% clean up java objects
arrayfun(@(x) x.reset, allagents);
nl.coenvl.sam.ExperimentControl.ResetExperiment();

end

function cost = getCost(constraint)
%% Get solution costs

cost = sum(cellfun(@(x) x.getExternalCost(), constraint));

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

