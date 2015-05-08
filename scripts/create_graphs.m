%#ok<*FNDSB>

%% Make figure of solution and messages
% load data\exp_delaunayGraph_3_20150109T153015_i10_results.mat
% load data\exp_gametheory_delaunayGraph_3_20150109T164039_i10_results.mat

close all;
%%
% Settings:
font = 'times';
% font = 'courier';
% font = 'newcenturyschlbk';
% font = 'Helvetica';

% Do big image and have latex resize it. Looks nicer
titlesize = 18;
labelsize = 16;
legendsize = 14;
axissize = 14;
linewidth = 2;
decolinewidth = .25; %lines of box, legend, axis etc.
figwidth = 16; %centimeters
figheight = 7.5; %centimeters
titleweight = 'bold';

boxonoff = 'on';
gridonoff = 'on';
minorgrid = 'off';
minortick = 'on';
legendbox = 'boxoff';
outputfolder = 'C:\Data\Dropbox\papers\CCOP\images';
% expname = [folder filesep 'exp_e10_trandomSample'];
expname = [outputfolder filesep 'gametheory_exp_e10_trandomSample'];
format = 'eps';
printoptions = {'-transparent', '-q105'}; %, 'r600', '-m2'};
myalgo = 'CoCoA';
plotfun = @semilogy; %except for cost (always plot)
doExport = false;

algos = fieldnames(results);
algos = [myalgo; algos(find(~strcmp(algos, myalgo)))];

% colors = hsv(numel(algos) + 1);
% colors = gray(numel(algos) + 1);
colors = cubehelix(numel(algos) + 1, .5, -1.5, 3, 1);
styles = repmat({'o-', '--', '-.', '-', ':'},1,ceil(numel(algos)/5));
styles{strcmp(algos, myalgo)} = 'o-';

%% Plot 1 - Solution cost

figure(187)
clf
set(gcf, 'Units', 'centimeters', 'Position', [2 2 figwidth figheight]);
hold on

for i = 1:numel(algos)
    y = sum(results.(algos{i}).costs, 2) ./ numel(settings.numExps);
    plot(settings.nagents, y, styles{i}, 'linewidth', linewidth, 'color', colors(i,:));
%     hold on;
end

hl = legend(algos{:}, 'Location', 'NorthWest');
legend(legendbox);

set(hl, 'fontsize', legendsize, 'fontname', font, 'linewidth', decolinewidth);
set(gca, 'fontsize', axissize, 'fontname', font, 'linewidth', decolinewidth, ...
    'YMinorGrid', minorgrid, 'YMinorTick', minortick, 'XMinorGrid', minorgrid, 'XMinorTick', minortick);

box(boxonoff);
grid(gridonoff);

% ht = title('Solution cost', 'fontsize', titlesize, 'fontname', font, 'fontweight', titleweight);
xlabel('Problem size', 'fontsize', labelsize, 'fontname', font);
ylabel('Cost', 'fontsize', labelsize, 'fontname', font);

if doExport 
    export_fig(gcf, fullfile(outputfolder, sprintf('%s_costs.%s', expname, format)), printoptions{:}); 
end
%% Plot 2 - Cost function evaluations

figure(188)
clf;

for i = 1:numel(algos)
    y = sum(results.(algos{i}).evals, 2) ./ numel(settings.numExps);
    plotfun(settings.nagents, y, styles{i}, 'linewidth', linewidth, 'color', colors(i,:));
    hold on;
end
hl2 = legend(algos{:}, 'Location', 'NorthWest');
legend(legendbox);

set(hl2, 'fontsize', legendsize, 'fontname', font, 'linewidth', decolinewidth);
set(gca, 'fontsize', axissize, 'fontname', font, 'linewidth', decolinewidth, ...
    'YMinorGrid', minorgrid, 'YMinorTick', minortick, 'XMinorGrid', minorgrid, 'XMinorTick', minortick);

box(boxonoff);
grid(gridonoff);

% title('Solution speed', 'fontsize', titlesize, 'fontname', font, 'fontweight', titleweight);
xlabel('Problem size', 'fontsize', labelsize, 'fontname', font);
ylabel('# Cost function evaluations', 'fontsize', labelsize, 'fontname', font);

set(gcf, 'Units', 'centimeters', 'Position', [2 2 figwidth figheight]);

if doExport 
    export_fig(gcf, fullfile(outputfolder, sprintf('%s_evals.%s', expname, format)), printoptions{:}); 
end

%% Plot 3 - Number of Messages

figure(189)
clf;
 
algos = algos(find(~strcmp(algos, 'CFL'))); %Because CFL does not send messages 

for i = 1:numel(algos)
    y = sum(results.(algos{i}).msgs, 2) ./ numel(settings.numExps);
    plotfun(settings.nagents, y, styles{i}, 'linewidth', linewidth, 'color', colors(i,:));
    hold on;
end
hl3 = legend(algos{:}, 'Location', 'NorthWest');
legend(legendbox);

set(hl3, 'fontsize', legendsize, 'fontname', font, 'linewidth', decolinewidth);
set(gca, 'fontsize', axissize, 'fontname', font, 'linewidth', decolinewidth, ...
    'YMinorGrid', minorgrid, 'YMinorTick', minortick, 'XMinorGrid', minorgrid, 'XMinorTick', minortick);

box(boxonoff);
grid(gridonoff);

% title('Communication costs', 'fontsize', titlesize, 'fontname', font, 'fontweight', titleweight);
xlabel('Problem size', 'fontsize', labelsize, 'fontname', font);
ylabel('# Messages transmitted', 'fontsize', labelsize, 'fontname', font);

set(gcf, 'Units', 'centimeters', 'Position', [2 2 figwidth figheight]);

if doExport 
    export_fig(gcf, fullfile(outputfolder, sprintf('%s_msgs.%s', expname, format)), printoptions{:}); 
end
fprintf('Done!\n');
