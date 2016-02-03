%load data\exp_LocalGameTheoreticCostFunction_delaunayGraph_i10_c3_t20150717T142129_results.mat
load data\exp_task_scheduling_randomGraph_5_20150716T131749_i10_results.mat
ratio = @(A,B) 100 * mean((mean(A,2) ./ mean(B,2)));
%%
mgm_cost_ratio = 100-ratio(results.CoCoA.costs, results.MGM2.costs)
mgm_msgs_ratio = ratio(results.CoCoA.msgs, results.MGM2.msgs)
mgm_evals_ratio = ratio(results.CoCoA.msgs, results.MGM2.evals)

dsa_cost_ratio = 100-ratio(results.CoCoA.costs, results.DSA.costs)
dsa_msgs_ratio = ratio(results.CoCoA.msgs, results.DSA.msgs)
dsa_evals_ratio = ratio(results.CoCoA.msgs, results.DSA.evals)