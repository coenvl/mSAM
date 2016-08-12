%% GETSOLVERCOUNTERPART
% *Summary of this function goes here*
%
% Detailed explanation goes here
%
%% Copyright
% * *2016 - TNO*
% * *Author*: Coen van Leeuwen
% * *Since*: August 12, 2016
%
%% See also:
%

%% Function Definition
function type = getSolverCounterPart(solverType)

dummyVariable = nl.coenvl.sam.variables.IntegerVariable(int32(1), int32(1));
dummyAgent = nl.coenvl.sam.agents.VariableAgent(dummyVariable, 'dummy');
dummySolver = feval(solverType, dummyAgent);

assert(isa(dummySolver, 'nl.coenvl.sam.solvers.MaxSumVariableSolver'), ...
    'EXPERIMENT:initConstraintAgents:INVALIDSOLVERTYPE', ...
    'Unexpected solver type, constraint agents only apply to MaxSum');

type = char(dummySolver.getCounterPart().getCanonicalName());

end