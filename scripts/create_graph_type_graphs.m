%#ok<*FNDSB>

%% Make figure of solution and messages
% load data\exp_realgraphs_3_20150125T140326_i5_results.mat

close all;
%%
% Settings:
font = 'times';
% font = 'courier';
% font = 'newcenturyschlbk';
% font = 'Helvetica';

% Do big image and have latex resize it. Looks nicer
titlesize = 18;
labelsize = 18;
legendsize = 14;
axissize = 12;
linewidth = 2;
decolinewidth = .25; %lines of box, legend, axis etc.
figwidth = 20; %centimeters
figheight = 4; %centimeters
titleweight = 'bold';

boxonoff = 'on';
gridonoff = 'off';
minorgrid = 'off';
minortick = 'off';
ticklength = [0 0];
legendLoc = 'NorthOutside';
legendOrientation = 'Horizontal';
repeatLegend = false;

seplinewidth = .1;
seplinecolor = [0 0 0];
seplinestyle = ':';

outputfolder = 'C:\Data\Dropbox\papers\CCOP\images';
expname = 'realGraphs';
format = 'eps';
printoptions = {'-transparent', '-q105'}; %, 'r600', '-m2'};
myalgo = 'CoCoA';
doExport = true;

graphTypes = fieldnames(results);
algos = fieldnames(solvers);

% colors = hsv(numel(algos) + 1);
% colors = gray(numel(algos) + 1);
% colors = [215 25 28; 253 174 97; 255 255 191; 44 123 182; 171 217 233] ./ 255;
colors = cubehelix(numel(algos) + 1, .5, -1.5, 3, 1);

%% Plot 1 - Solution cost

figure(187);
clf;
set(gcf, 'Units', 'centimeters', 'Position', [2 2 figwidth figheight]);
hold on

for i = 1:numel(graphTypes)
    for j = 1:numel(algos)
        y(j) = sum(results.(graphTypes{i}).(algos{j}).costs, 2) ./ numel(settings.numExps);
    end
    
    y = y / max(y);
    
    for j = 1:numel(algos)
        bar(j + (i-1) * (numel(algos) + 1),y(j), 'FaceColor', colors(j,:));
        hold on;
    end
    
    if (i ~= numel(graphTypes))
        endx = i * (numel(algos) + 1);
        line([endx endx], [0 1], 'LineWidth', seplinewidth, ...
            'Color', seplinecolor, 'LineStyle', seplinestyle);
    end
end

hl = legend(algos{:}, 'Location', legendLoc, 'Orientation', legendOrientation);
set(hl, 'fontsize', legendsize, 'fontname', font, 'linewidth', decolinewidth);

set(gca, 'fontsize', axissize, 'fontname', font, 'linewidth', decolinewidth, ...
    'YMinorGrid', minorgrid, 'YMinorTick', minortick, ...
    'YTick', [0 1], 'XTick', (1:numel(graphTypes)) * (numel(algos) + 1) - (numel(algos)/2 + .5), ...
    'TickLength', ticklength, 'XTickLabel', graphTypes);

box(boxonoff);
grid(gridonoff);

% ht = title(graphTypes{i}, 'fontsize', titlesize, 'fontname', font, 'fontweight', titleweight);
% xlabel('Problem size', 'fontsize', labelsize, 'fontname', font);
% ylabel('Costs''', 'fontsize', labelsize, 'fontname', font);

if doExport 
    export_fig(gcf, fullfile(outputfolder, sprintf('%s_costs.%s', expname, format)), printoptions{:}); 
end

%% Plot 2 - Solution evals

figure(188);
clf;
set(gcf, 'Units', 'centimeters', 'Position', [2 2 figwidth figheight]);
hold on

for i = 1:numel(graphTypes)
    for j = 1:numel(algos)
        y(j) = sum(results.(graphTypes{i}).(algos{j}).evals, 2) ./ numel(settings.numExps);
    end
    
    y = y / max(y);
    
    for j = 1:numel(algos)
        bar(j + (i-1) * (numel(algos) + 1),y(j), 'FaceColor', colors(j,:));
        hold on;
    end
    
    if (i ~= numel(graphTypes))
        endx = i * (numel(algos) + 1);
        line([endx endx], [0 1], 'LineWidth', seplinewidth, ...
            'Color', seplinecolor, 'LineStyle', seplinestyle);
    end
end

if repeatLegend
    hl = legend(algos{:}, 'Location', legendLoc, 'Orientation', legendOrientation);
    set(hl, 'fontsize', legendsize, 'fontname', font, 'linewidth', decolinewidth);
else
   set(gcf, 'Position', [2 2 figwidth figheight*0.7]); 
end

set(gca, 'fontsize', axissize, 'fontname', font, 'linewidth', decolinewidth, ...
    'YMinorGrid', minorgrid, 'YMinorTick', minortick, ...
    'YTick', [0 1], 'XTick', (1:numel(graphTypes)) * (numel(algos) + 1) - (numel(algos)/2 + .5), ...
    'TickLength', ticklength, 'XTickLabel', graphTypes);

box(boxonoff);
grid(gridonoff);

% ht = title(graphTypes{i}, 'fontsize', titlesize, 'fontname', font, 'fontweight', titleweight);
% xlabel('Problem size', 'fontsize', labelsize, 'fontname', font);
% ylabel('#Evaluations normalized', 'fontsize', labelsize, 'fontname', font);

if doExport 
    export_fig(gcf, fullfile(outputfolder, sprintf('%s_evals.%s', expname, format)), printoptions{:}); 
end