k = findstr('.', options.constraint.type);
costType = options.constraint.type(k(end)+1:end);

filename = sprintf('results_%s_%s_%s_i%d_d%d_n%d_t%s.mat', class(exp), costType, func2str(settings.graphType), settings.numExps, settings.ncolors, settings.nagents, datestr(now, 30))

save(fullfile('data', settings.series, filename), 'settings', 'options', 'solvers', 'results');
