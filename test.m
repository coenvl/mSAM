distance = .5;
pos = rand(15,3);

%%

D = pdist2(pos, pos);
D(D == 0) = Inf;
[I,J] = find(D < distance);

e = [I J];


%%

options.nAgents = uint16(100);
options.maxDist = 90;
options.scale = 200;

[edges, pos] = distanceGraph3(options);

clf;
scatter3(pos(:,1), pos(:,2), pos(:,3))
axis equal;
hold on;
for i = 1:size(edges,1); 
    line(pos(edges(i,:),1),pos(edges(i,:),2),pos(edges(i,:),3)); 
end

%% Test threads

variable = nl.coenvl.sam.variables.IntegerVariable(0,10);
agent = nl.coenvl.sam.agents.LocalSolverAgent('ThreadTestAgent', variable);
costfun = nl.coenvl.sam.costfunctions.LocalInequalityConstraintCostFunction(agent);
solver = nl.coenvl.sam.solvers.UniqueFirstCooperativeSolver(agent, costfun);

agent.setSolver(solver);
agent.init();

agent.reset();
nl.coenvl.sam.ExperimentControl.ResetExperiment()
