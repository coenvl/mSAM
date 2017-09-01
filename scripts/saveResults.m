k = strfind(options.constraint.type, '.');
costType = options.constraint.type(k(end)+1:end);

if isfield(settings, 'name')
    problemname = settings.name;
else
    problemname = class(exp);
end

filename = sprintf('results_%s_%s_%s_i%d_d%d_n%d_t%s.mat', ...
    problemname, costType, func2str(settings.graphType), ...
    settings.numExps, settings.ncolors, settings.nagents, datestr(now, 30))

path = fullfile('data', settings.series);
if ~exist(path, 'dir')
    mkdir(path);
end

savevars = intersect(who, ...
    {'settings', 'options', 'solvers', 'results', 'iterSolver', 'initSolver'});
save(fullfile(path, filename), savevars{:});
