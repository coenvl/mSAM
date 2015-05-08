%% POISSONSAMPLE
% *Select 2D points according to a poisson point process*
%
% Select a set of 2D points according to a poisson point process. In doing
% so, the average density in any subregion of the space is constant. More
% info see: http://en.wikipedia.org/wiki/Poisson_point_process. 
% 
% Iteratively selects previously generated points and generates K random
% points in a circle around it. Whenever any of the K points is more than 1
% unit distance from its nearest neighbor, it is added it to the set. When 
% none of the K meets this criterium, the selected point will not be
% selected in the future again.
%
%% Usage
% pos = poissonSample(n)
%   Generates N 2D points as a Nx2 matrix
%   
%% Copyright
% * *2015 - TNO*
% * *Author*: Coen van Leeuwen
% * *Since*: May 08, 2015

%% Function Definition
function pos = poissonSample(n)

% Maximum number of samples before rejection
k = 30; 

% Initialize
pos = [0 0];
queue = 1;

% Iterate until we have as many points as we need
while (size(pos, 1) < n && ~isempty(queue))
    
    % Pick the item from the queue which is nearest to origin
    %[~, i] = min(pdist2([0 0], pos(queue, :)));
    [~, i] = min(hypot(pos(queue, 1), pos(queue,2)));
    
    % Create k samples around data[i];
    found = false;
    for j = 1:k
        a = 2 * pi * rand;
        ri = rand + 1;
        x = pos(queue(i),1) + ri * cos(a);
        y = pos(queue(i),2) + ri * sin(a);
        
        %if (min(pdist2([x,y], pos)) > 1)
        if (min(hypot(pos(:,1) - x, pos(:,2) - y)) > 1)
            pos = [pos; x y];                           %#ok<AGROW>
            queue(end + 1) = length(pos);               %#ok<AGROW>
            found = true;
            break;
        end
    end
    
    % If no points can be found around the current item, remove from the
    % queue
    if ~found
        queue(i) = [];                                  %#ok<AGROW>
    end
end