%--------------------------------------------------------
% returns power levels according to the simple path loss
%--------------------------------------------------------
function power_factor = path_loss_factor(distances)
    alpha = 100;
    beta = 100;
    power_factor = alpha./(distances+beta).^2;
end