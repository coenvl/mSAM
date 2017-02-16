%% GETEXPERIMENTSOLVERS
% *Summary of this function goes here*
%
% Detailed explanation goes here
%
%% Copyright
% * *2016 - TNO*
% * *Author*: Coen van Leeuwen
% * *Since*: September 15, 2016
% 
%% See also:
%

%% Function Definition
function solvers = getExperimentSolvers(series)

%%
if nargin > 0 && ~isempty(series)
    switch series
        case 'aaai17'
            idx = [1 2 3 7 9 11 13];
        case 'ijcai17'
            idx = 1:14;
        case 'wpt'
            idx = 13;
        case 'hybrid+'
            idx = [3 5 7 9 11 16:35];
        otherwise
            error('Unknown series %s', series);
    end
end

solvers(1).name = 'CoCoA';
solvers(1).initSolverType = 'nl.coenvl.sam.solvers.CoCoSolver';
solvers(1).iterSolverType = '';

solvers(2).name = 'CoCoA_UF';
solvers(2).initSolverType = 'nl.coenvl.sam.solvers.CoCoASolver';
solvers(2).iterSolverType = '';

solvers(3).name = 'ACLS';
solvers(3).initSolverType = '';
solvers(3).iterSolverType = 'nl.coenvl.sam.solvers.ACLSSolver';

solvers(4).name = 'CoCoA - ACLS';
solvers(4).initSolverType = 'nl.coenvl.sam.solvers.CoCoASolver';
solvers(4).iterSolverType = 'nl.coenvl.sam.solvers.ACLSSolver';

solvers(5).name = 'ACLSUB';
solvers(5).initSolverType = '';
solvers(5).iterSolverType = 'nl.coenvl.sam.solvers.ACLSUBSolver';

solvers(6).name = 'CoCoA - ACLSUB';
solvers(6).initSolverType = 'nl.coenvl.sam.solvers.CoCoASolver';
solvers(6).iterSolverType = 'nl.coenvl.sam.solvers.ACLSUBSolver';

solvers(7).name = 'DSA';
solvers(7).initSolverType = '';
solvers(7).iterSolverType = 'nl.coenvl.sam.solvers.DSASolver';

solvers(8).name = 'CoCoA - DSA';
solvers(8).initSolverType = 'nl.coenvl.sam.solvers.CoCoASolver';
solvers(8).iterSolverType = 'nl.coenvl.sam.solvers.DSASolver';

solvers(9).name = 'MCSMGM';
solvers(9).initSolverType = '';
solvers(9).iterSolverType = 'nl.coenvl.sam.solvers.MCSMGMSolver';

solvers(10).name = 'CoCoA - MCSMGM';
solvers(10).initSolverType = 'nl.coenvl.sam.solvers.CoCoASolver';
solvers(10).iterSolverType = 'nl.coenvl.sam.solvers.MCSMGMSolver';

solvers(11).name = 'MGM2';
solvers(11).initSolverType = '';
solvers(11).iterSolverType = 'nl.coenvl.sam.solvers.MGM2Solver';

solvers(12).name = 'CoCoA - MGM2';
solvers(12).initSolverType = 'nl.coenvl.sam.solvers.CoCoASolver';
solvers(12).iterSolverType = 'nl.coenvl.sam.solvers.MGM2Solver';

solvers(13).name = 'Max-Sum';
solvers(13).initSolverType = '';
solvers(13).iterSolverType = 'nl.coenvl.sam.solvers.MaxSumVariableSolver';

solvers(14).name = 'Max-Sum_ADVP';
solvers(14).initSolverType = '';
solvers(14).iterSolverType = 'nl.coenvl.sam.solvers.MaxSumADVPVariableSolver';

solvers(15).name = 'CoCoA - Max-Sum_ADVP';
solvers(15).initSolverType = 'nl.coenvl.sam.solvers.CoCoASolver';
solvers(15).iterSolverType = 'nl.coenvl.sam.solvers.MaxSumADVPVariableSolver';

solvers(16).name = 'ACLS_ACLSUB';
solvers(16).initSolverType = 'nl.coenvl.sam.solvers.ACLSSolver';
solvers(16).iterSolverType = 'nl.coenvl.sam.solvers.ACLSUBSolver';

solvers(17).name = 'ACLS_DSA';
solvers(17).initSolverType = 'nl.coenvl.sam.solvers.ACLSSolver';
solvers(17).iterSolverType = 'nl.coenvl.sam.solvers.DSASolver';

solvers(18).name = 'ACLS_MGM2';
solvers(18).initSolverType = 'nl.coenvl.sam.solvers.ACLSSolver';
solvers(18).iterSolverType = 'nl.coenvl.sam.solvers.MGM2Solver';

solvers(19).name = 'ACLS_MCSMGM';
solvers(19).initSolverType = 'nl.coenvl.sam.solvers.ACLSSolver';
solvers(19).iterSolverType = 'nl.coenvl.sam.solvers.MCSMGMSolver';

solvers(20).name = 'ACLSUB_ACLS';
solvers(20).initSolverType = 'nl.coenvl.sam.solvers.ACLSUBSolver';
solvers(20).iterSolverType = 'nl.coenvl.sam.solvers.ACLSSolver';

solvers(21).name = 'ACLSUB_DSA';
solvers(21).initSolverType = 'nl.coenvl.sam.solvers.ACLSUBSolver';
solvers(21).iterSolverType = 'nl.coenvl.sam.solvers.DSASolver';

solvers(22).name = 'ACLSUB_MGM2';
solvers(22).initSolverType = 'nl.coenvl.sam.solvers.ACLSUBSolver';
solvers(22).iterSolverType = 'nl.coenvl.sam.solvers.MGM2Solver';

solvers(23).name = 'ACLSUB_MCSMGM';
solvers(23).initSolverType = 'nl.coenvl.sam.solvers.ACLSUBSolver';
solvers(23).iterSolverType = 'nl.coenvl.sam.solvers.MCSMGMSolver';

solvers(24).name = 'DSA_ACLS';
solvers(24).initSolverType = 'nl.coenvl.sam.solvers.DSASolver';
solvers(24).iterSolverType = 'nl.coenvl.sam.solvers.ACLSSolver';

solvers(25).name = 'DSA_ACLSUB';
solvers(25).initSolverType = 'nl.coenvl.sam.solvers.DSASolver';
solvers(25).iterSolverType = 'nl.coenvl.sam.solvers.ACLSUBSolver';

solvers(26).name = 'DSA_MGM2';
solvers(26).initSolverType = 'nl.coenvl.sam.solvers.DSASolver';
solvers(26).iterSolverType = 'nl.coenvl.sam.solvers.MGM2Solver';

solvers(27).name = 'DSA_MCSMGM';
solvers(27).initSolverType = 'nl.coenvl.sam.solvers.DSASolver';
solvers(27).iterSolverType = 'nl.coenvl.sam.solvers.MCSMGMSolver';

solvers(28).name = 'MGM2_ACLS';
solvers(28).initSolverType = 'nl.coenvl.sam.solvers.MGM2Solver';
solvers(28).iterSolverType = 'nl.coenvl.sam.solvers.ACLSSolver';

solvers(29).name = 'MGM2_ACLSUB';
solvers(29).initSolverType = 'nl.coenvl.sam.solvers.MGM2Solver';
solvers(29).iterSolverType = 'nl.coenvl.sam.solvers.ACLSUBSolver';

solvers(30).name = 'MGM2_DSA';
solvers(30).initSolverType = 'nl.coenvl.sam.solvers.MGM2Solver';
solvers(30).iterSolverType = 'nl.coenvl.sam.solvers.DSASolver';

solvers(31).name = 'MGM2_MCSMGM';
solvers(31).initSolverType = 'nl.coenvl.sam.solvers.MGM2Solver';
solvers(31).iterSolverType = 'nl.coenvl.sam.solvers.MCSMGMSolver';

solvers(32).name = 'MCSMGM_ACLS';
solvers(32).initSolverType = 'nl.coenvl.sam.solvers.MCSMGMSolver';
solvers(32).iterSolverType = 'nl.coenvl.sam.solvers.ACLSSolver';

solvers(33).name = 'MCSMGM_ACLSUB';
solvers(33).initSolverType = 'nl.coenvl.sam.solvers.MCSMGMSolver';
solvers(33).iterSolverType = 'nl.coenvl.sam.solvers.ACLSUBSolver';

solvers(34).name = 'MCSMGM_DSA';
solvers(34).initSolverType = 'nl.coenvl.sam.solvers.MCSMGMSolver';
solvers(34).iterSolverType = 'nl.coenvl.sam.solvers.DSASolver';

solvers(35).name = 'MCSMGM_MGM2';
solvers(35).initSolverType = 'nl.coenvl.sam.solvers.MCSMGMSolver';
solvers(35).iterSolverType = 'nl.coenvl.sam.solvers.MGM2Solver';

%% Select based on series
if exist('idx', 'var')
    solvers = solvers(idx);
end

%% Remove empty entries
k = arrayfun(@(x) ~isempty(x.name), solvers);
solvers = solvers(k);