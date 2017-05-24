%--------------------------------------------------------
% Centralized linear programming solution
%--------------------------------------------------------
function [x,fval] = LP_solution(transmitter_to_receiver_distances,transmitter_to_sensor_distances,...
    EMR_Threshold,MAX_POWER, MIN_POWER)

    [Num_Transmitters,~] = size(transmitter_to_receiver_distances);
    [~,Num_Sensors] = size(transmitter_to_sensor_distances);
    
    % f'.x : total transmitted power
    f = -sum(path_loss_factor(transmitter_to_receiver_distances),2);
    
    % A.x <= b = EMR threshold 
    A = transpose(path_loss_factor(transmitter_to_sensor_distances));
    b = repmat(EMR_Threshold,Num_Sensors,1);
    
    % lower and upper bounds for the power search space
    lb = repmat(MIN_POWER,Num_Transmitters,1);
    ub = repmat(MAX_POWER,Num_Transmitters,1);
      
    [x,fval] = linprog(f,A,b,[],[],lb,ub);
end
