function filename = prettyExportFig(name, fig, ax, label_x, label_y, hl, titles, graphOptions)
%PRETTYEXPORTFIG Summary of this function goes here
%   Detailed explanation goes here
    set(fig, 'Units', 'centimeters', 'Position', [3 3 graphOptions.figure.width graphOptions.figure.height]);
    
    if ~isempty(hl)
        set(hl, 'fontsize', graphOptions.legend.fontsize, 'fontname', graphOptions.legend.font, 'linewidth', ...
            graphOptions.legend.linewidth, 'Box', graphOptions.legend.box, 'Location', graphOptions.legend.location, 'Interpreter', 'none');
    end
    
    for i = 1:numel(ax)
    set(ax(i), 'fontsize', graphOptions.axes.fontsize, 'fontname', graphOptions.axes.font, 'linewidth', graphOptions.axes.linewidth, ...
        'YMinorGrid', graphOptions.axes.minorgrid, 'YMinorTick', graphOptions.axes.minortick, ...
        'XMinorGrid', graphOptions.axes.minorgrid, 'XMinorTick', graphOptions.axes.minortick, ...
        'Box', graphOptions.axes.box, 'YGrid', graphOptions.axes.grid, 'XGrid', graphOptions.axes.grid, ... 
        'YScale', graphOptions.axes.yscale,  'XScale', graphOptions.axes.xscale);
    end
    
    if ~isempty(label_x)
        set(label_x, 'fontsize', graphOptions.label.fontsize, 'fontname', graphOptions.label.font);
    end
    
    if ~isempty(label_y)
        set(label_y, 'fontsize', graphOptions.label.fontsize, 'fontname', graphOptions.label.font);
    end
    
    if ~isempty(titles)
        for i = 1:numel(titles)
        set(titles(i), 'fontsize', graphOptions.title.fontsize, 'fontname', graphOptions.title.font);
        end    
    end
        
    filename = sprintf('%s.%s', name, graphOptions.export.format); 
    export_fig(fig, fullfile(graphOptions.export.folder, filename), graphOptions.export.arguments{:});
end

