%% GETGRAPHOPTIONS
% *Summary of this function goes here*
%
% Detailed explanation goes here
%
%% Copyright
% * *2015 - TNO*
% * *Author*: Coen van Leeuwen
% * *Since*: July 30, 2015
% 
%% See also:
%

%% Function Definition
function options = getGraphOptions()

scale_factor = 1;
font_scale_factor = 1;

% Figure options, do big image and have latex resize it. Looks nicer
options.figure.units = 'centimeters';
options.figure.width = scale_factor * 16;
options.figure.height = scale_factor * 12; %6.5;

% Label options
options.label.font = 'times';
options.label.fontsize = font_scale_factor * 16;

% Legend options
options.legend.font = 'times';
options.legend.fontsize = font_scale_factor * 12;
options.legend.box = 'off';
options.legend.linewidth = scale_factor * .25;
options.legend.location = 'NorthWest';

% Axes options
options.axes.font = 'times';
options.axes.fontsize = font_scale_factor * 14;
options.axes.box = 'on';
options.axes.grid = 'on';
options.axes.minorgrid = 'off';
options.axes.minortick = 'on';
options.axes.linewidth = scale_factor * .25;

% Export options
options.export.do = true;
options.export.format = 'png'; %'eps';
options.export.arguments = {'-transparent', '-q105', '-r600', '-m2'};
%options.export.folder = 'C:\Data\Dropbox\papers\CCOP\images';
options.export.folder = fullfile(pwd, 'figures');

% Plot options
options.plot.emphasize = 'CoCoA';
options.plot.y_fun = @(x)mean(x,2);
% options.plot.colors = cubehelix(6, .5, -1.5, 3, 1);
options.plot.linewidth = scale_factor * 2;

% On the error bar:
options.plot.errorlinewidth = scale_factor * 1;
sec = @(x)x(:,2); 
secmax = @(x)(sec(sort(x,2,'descend'))); 
secmin = @(x)(sec(sort(x,2,'ascend')));
options.plot.hi_error_fun = @(x) secmax(x);
options.plot.low_error_fun = @(x) secmin(x);
