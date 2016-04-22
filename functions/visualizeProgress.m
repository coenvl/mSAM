%% visualizeProgress
% *Summary of this function goes here*
%
% Detailed explanation goes here
%
%% Copyright
% * *2016 - TNO*
% * *Author*: Coen van Leeuwen
% * *Since*: April 22, 2016
% 
%% See also:
%

%% Function Definition
function visualizeProgress(exp, solvername)

persistent handles;
persistent legendentries;

% Visualize data
ydata = exp.allcost;
xdata = exp.alltimes;

if numel(ydata) == 1
    style = {'LineStyle', 'none', 'Marker', 'o'};
else
    style = {'LineStyle', '-', 'Marker', 'none'};
end

if isempty(handles) || ~isfield(handles, 'fig') || ~ishandle(handles.fig)
    % Create figure
    handles.fig = figure(007);
    handles.ax = gca(handles.fig);
    hold(handles.ax, 'on');
    legendentries = {};
end

if ~isfield(handles, solvername) || ~ishandle(handles.(solvername))
    handles.(solvername) = plot(xdata, ydata, 'parent', handles.ax, style{:});
    legendentries = [legendentries solvername];
    %handles.legend = legend(handles.ax, solvertypes);
else
    set(handles.(solvername), 'XData', xdata, 'YData', ydata, style{:});
end

legend(handles.ax, legendentries);
drawnow;

end