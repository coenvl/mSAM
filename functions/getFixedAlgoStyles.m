%% GETFIXEDALGOSTYLES
% *Summary of this function goes here*
%
% Detailed explanation goes here
%
%% Copyright
% * *2017 - TNO*
% * *Author*: leeuwencjv
% * *Since*: January 06, 2017
% 
%% See also:
%

%% Function Definition
function fixedStyles = getFixedAlgoStyles()

colors = cubehelix(8, .5, -1.5, 3, 1);
fixedStyles.CoCoA =              {'Marker', 'o',    'LineStyle', 'none', 'Color', colors(1,:)};
fixedStyles.ACLS =               {'Marker', 'none', 'LineStyle', '-',    'Color', colors(2,:)};
fixedStyles.CoCoA_ACLS =         {'Marker', 'none', 'LineStyle', ':',    'Color', colors(2,:)};
fixedStyles.ACLSUB =             {'Marker', 'none', 'LineStyle', '--',   'Color', colors(3,:)};
fixedStyles.CoCoA_ACLSUB =       {'Marker', 'none', 'LineStyle', '-.',   'Color', colors(3,:)};
fixedStyles.DSA =                {'Marker', 'none', 'LineStyle', '-',    'Color', colors(4,:)};
fixedStyles.CoCoA_DSA =          {'Marker', 'none', 'LineStyle', ':',    'Color', colors(4,:)};
fixedStyles.MGM2 =               {'Marker', 'none', 'LineStyle', '--',   'Color', colors(6,:)};
fixedStyles.CoCoA_MGM2 =         {'Marker', 'none', 'LineStyle', '-.',   'Color', colors(6,:)};
fixedStyles.MCSMGM =             {'Marker', 'none', 'LineStyle', '-',    'Color', colors(5,:)};
fixedStyles.CoCoA_MCSMGM =       {'Marker', 'none', 'LineStyle', ':',    'Color', colors(5,:)};
fixedStyles.Max_Sum =            {'Marker', 'none', 'LineStyle', '--',   'Color', colors(7,:)};
fixedStyles.CoCoA_Max_Sum =      {'Marker', 'none', 'LineStyle', '-.',   'Color', colors(7,:)};
fixedStyles.CoCoA_Max_Sum_ADVP = {'Marker', 'none', 'LineStyle', ':',    'Color', colors(8,:)};
