%% GENERATEWPTSCENARIO
% *Summary of this function goes here*
%
% Detailed explanation goes here
%
%% Copyright
% * *2017 - TNO*
% * *Author*: leeuwencjv
% * *Since*: March 24, 2017
% 
%% See also:
%

%% Function Definition
function [agentPos, receiverPos, sensorPos, edges, transmitter_to_receiver, transmitter_to_sensor] = ...
    generateWPTScenario( nagents, nreceivers, nsensors )

% Ballpark something that might work
MAX_ARITY = ceil(log10(nagents)+1);

tp = 100 * (rand(1,2) + poissonSample(nagents));
tp = tp - min(tp);
for i = nagents:-1:1
    agentPos{i} = tp(i,:);
end

rp = rand(1,2) + poissonSample(nreceivers);
rp = rp - min(rp);
rp = rp .* (max(tp) ./ max(rp)); % Scale to same area as transmitters
for i = nreceivers:-1:1
    receiverPos{i} = rp(i,:);
end

sp = 100 * (rand(1,2) + poissonSample(nsensors));
sp = sp - min(sp);
sp = sp .* (max(tp) ./ max(sp)); % Scale to same area as transmitters
for i = nsensors:-1:1
    sensorPos{i} = sp(i,:);
end

transmitter_to_receiver = sqrt(pdist2(tp, rp));
transmitter_to_sensor = sqrt(pdist2(tp, sp));

receiverEdges = getEdgesFromDistances(transmitter_to_receiver, MAX_ARITY);
sensorEdges = getEdgesFromDistances(transmitter_to_sensor, MAX_ARITY);

% % Every receiver and every sensor is going to be an edge
% t2r_mask = transmitter_to_receiver < 20 * min(transmitter_to_receiver,[],1);
% t2s_mask = transmitter_to_sensor < 20 * min(transmitter_to_sensor,[],1);
% 
% [a,b] = find(t2r_mask);
% receiverEdges = cell(1, nreceivers);
% for i = 1:numel(a)
%     r = b(i);
%     if isempty(receiverEdges{r})
%         receiverEdges{r} = a(i);
%     elseif numel(receiverEdges{r}) < MAX_ARITY
%         receiverEdges{r} = [receiverEdges{r} a(i)];
%     end
% end
% 
% [c,d] = find(t2s_mask);
% sensorEdges = cell(1, nsensors);
% for i = 1:numel(c)
%     r = d(i);
%     if isempty(sensorEdges{r})
%         sensorEdges{r} = c(i);
%     elseif numel(sensorEdges{r}) < MAX_ARITY
%         sensorEdges{r} = [sensorEdges{r} c(i)];
%     end
% end

edges = {receiverEdges, sensorEdges};


%% Do not use previous mask again, we might have reduced the arity

% receiver..
t2r_mask = false(nagents, nreceivers);
for a = 1:numel(receiverEdges)
    for r = receiverEdges{a}
        t2r_mask(r,a) = true;
    end
end

% sensor..
t2s_mask = false(nagents, nsensors);
for a = 1:numel(sensorEdges)
    for s = sensorEdges{a}
        t2s_mask(s,a) = true;
    end
end

transmitter_to_receiver(~t2r_mask) = inf;
transmitter_to_sensor(~t2s_mask) = inf;

fprintf('generated problem\n');
return;
% agentPos = {[10 23], [40 80], [65 45], [130 30]};
% receiverPos = {[15 60], [50 20], [80 80], [103 52], [128 80]};
% sensorPos = {[35 50], [92 20], [108 80]};

% Receiver edges - Sensor edges
% edges = {{1, [1 3], [2 3], [3 4], 4}, {[1 2 3], [3 4], 3}};

%% Plot the problem
clf;
scatter(tp(:,1), tp(:,2), 'ro')
hold on
scatter(rp(:,1), rp(:,2), 'bs')
scatter(sp(:,1), sp(:,2), 'kv')

for r = 1:nreceivers
    e = receiverEdges{r};
    for i = e
        plot([tp(i,1) rp(r,1)], [tp(i,2) rp(r,2)], 'k:');
    end
end

for s = 1:nsensors
    e = sensorEdges{s};
    for i = e
        plot([tp(i,1) sp(s,1)], [tp(i,2) sp(s,2)], 'k:');
    end
end

legend('Transmitter', 'Receiver', 'Sensor');

end

function edges = getEdgesFromDistances(transmitter_distances, MAX_ARITY)

% Every column will be 1 constraint
ordered = sort(transmitter_distances, 1);
max_dist = mean(ordered(MAX_ARITY,:));
for c = size(transmitter_distances, 2):-1:1
    % Find sensors that are closer than the average MAX_ARITY-closest
    nearbyTransmitters = find(transmitter_distances(:,c) < max_dist);
    if numel(nearbyTransmitters) > MAX_ARITY
        % If necessary cap to the MAX_ARITY-closest
        [~, order] = sort(transmitter_distances(nearbyTransmitters, c));
        nearbyTransmitters = nearbyTransmitters(order(1:MAX_ARITY));
    end
    edges{c} = nearbyTransmitters';
end

end

