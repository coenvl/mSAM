function cellResults = fixSleepyLaptop(cellResults)
%FIXSLEEPYLAPTOP Sleeping mode of laptop influences experiment times...

algos = fieldnames(cellResults)';
for i = 1:numel(algos)
    % Fix the times matrix for cases when the computer went to sleep mode
    dt = cellfun(@diff, cellResults.(algos{i}).times, 'UniformOutput', false);
    k = cellfun(@(x) any(x > 10), dt);
    if any(k)
        a = find(k);
        
        % Assume actual time was average time of previous set of iterations
        dt_good = cellfun(@diff, cellResults.(algos{i}).times(~k), 'UniformOutput', false);
        t_new = mean(cellfun(@mean, dt_good));

        for j = 1:numel(a)
            b = find(dt{a(j)} > 10);
            t_err = diff(cellResults.(algos{i}).times{a(j)}(b + [0 1])) - t_new;
            cellResults.(algos{i}).times{a(j)}((b+1):end) = cellResults.(algos{i}).times{a(j)}((b+1):end) - t_err;
        end
    end
end


end

