function total = getNumberOfEvals()

total = nl.coenvl.sam.costfunctions.LocalInequalityConstraintCostFunction.nComparisons +...
    nl.coenvl.sam.costfunctions.InequalityConstraintCostFunction.nComparisons;