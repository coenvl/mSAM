%% PRINTGRAPH
% *Summary of this function goes here*
%
% Detailed explanation goes here
%
%% Copyright
% * *2015 - TNO*
% * *Author*: Coen van Leeuwen
% * *Since*: July 03, 2015
% 
%% See also:
%

%% Function Definition
function filename = printGraph(edges, filename)

if nargin < 2
    filename = uiputfile('*.png', 'Select output image');
end

[path,name,~] = fileparts(filename);

gvfile = fullfile(path,[name '.gv']);

fid = fopen(gvfile', 'w');
fprintf(fid, 'graph {\n');
fprintf(fid, 'graph [K=0.03, overlap=false];\n');
fprintf(fid, 'node [style=filled];\n');

for i = 1:max(max(edges))
    fprintf(fid, ' %d [fillcolor=white];\n', i);
end

for i = 1:size(edges,1)
    fprintf(fid, ' %d -- %d;\n', edges(i,1), edges(i,2));
end
fprintf(fid, '}\n');
fclose(fid);

graphvizpath = getGraphvizLoc();

layout = fullfile(graphvizpath, 'sfdp.exe');
% layout = fullfile(graphvizpath, 'twopi.exe');

system(sprintf('%s %s -T png -o %s', layout, gvfile, filename));

function loc = getGraphvizLoc()

loc = getenv('path_graphviz');
if ~isempty(loc), return; end

if ispc()
    for bp = {'C:/Progra~1', 'C:/Progra~2'}
        a = dir(fullfile(bp{:}, 'Graphviz*'));
        if ~isempty(a) && exist(fullfile(bp{:}, a.name, 'bin', 'gvedit.exe'), 'file');
            loc = fullfile(bp{:}, a.name, 'bin');
            setenv('path_graphviz', loc);
            return;
        end
    end
end

error('Could not automatically find Graphviz location, please set the environment variables ''path_graphviz''');