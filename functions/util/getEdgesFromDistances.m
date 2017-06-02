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
