function [ output_args ] = prettyExportFig(name, fig, ax, label_x, label_y, hl, graphOptions)
%PRETTYEXPORTFIG Summary of this function goes here
%   Detailed explanation goes here
    set(fig, 'Units', 'centimeters', 'Position', [3 3 graphOptions.figure.width graphOptions.figure.height]);
    
    if ~isempty(hl)
        set(hl, 'fontsize', graphOptions.legend.fontsize, 'fontname', graphOptions.legend.font, 'linewidth', ...
            graphOptions.legend.linewidth, 'Box', graphOptions.legend.box, 'Location', graphOptions.legend.location, 'Interpreter', 'none');
    end
    
    set(ax, 'fontsize', graphOptions.axes.fontsize, 'fontname', graphOptions.axes.font, 'linewidth', graphOptions.axes.linewidth, ...
        'YMinorGrid', graphOptions.axes.minorgrid, 'YMinorTick', graphOptions.axes.minortick, ...
        'XMinorGrid', graphOptions.axes.minorgrid, 'XMinorTick', graphOptions.axes.minortick, ...
        'Box', graphOptions.axes.box, 'YGrid', graphOptions.axes.grid, 'XGrid', graphOptions.axes.grid, ... 
        'YScale', graphOptions.axes.yscale,  'XScale', graphOptions.axes.xscale);
    
    set(label_x, 'fontsize', graphOptions.label.fontsize, 'fontname', graphOptions.label.font);
    set(label_y, 'fontsize', graphOptions.label.fontsize, 'fontname', graphOptions.label.font);
    
    filename = sprintf('%s.%s', name, graphOptions.export.format); 
    export_fig(fig, fullfile(graphOptions.export.folder, filename), graphOptions.export.arguments{:});
end

