%% POISSONSAMPLE
% *Select 3D points according to a poisson point process*
%
% Select a set of 3D points according to a poisson point process. In doing
% so, the average density in any subregion of the space is constant. More
% info see: http://en.wikipedia.org/wiki/Poisson_point_process. 
% 
% Iteratively selects previously generated points and generates K random
% points in a sphere around it. Whenever any of the K points is more than 1
% unit distance from its nearest neighbor, it is added it to the set. When 
% none of the K meets this criterium, the selected point will not be
% selected in the future again.
%
%% Usage
% pos = poissonSample3(n)
%   Generates N 3D points as a Nx3 matrix
%   
%% Copyright
% * *2015 - TNO*
% * *Author*: Coen van Leeuwen
% * *Since*: July 9, 2015

%% Function Definition
function pos = poissonSample3(n)

% Maximum number of samples before rejection
k = 30; 

% Initialize
pos = [0 0 0];
queue = 1;

% Iterate until we have as many points as we need
while (size(pos, 1) < n && ~isempty(queue))
    
    % Pick the item from the queue which is nearest to origin
    %[~, i] = min(pdist2([0 0], pos(queue, :)));
    [~, i] = min(pos(queue, 1).^2 + pos(queue,2).^2 + pos(queue,3).^2); % Actual distance not required...
    
    % Create k samples around data[i];
    found = false;
    for j = 1:k
        % Random point picking according to http://mathworld.wolfram.com/SpherePointPicking.html
        a = 2 * pi * rand;
        b = acos(2 * rand - 1);
        ri = rand + 1;
        
        x = pos(queue(i),1) + ri * cos(a) * sin(b);
        y = pos(queue(i),2) + ri * sin(a) * sin(b);
        z = pos(queue(i),3) + ri * cos(b);
        
        %if (min(pdist2([x,y], pos)) > 1)
        D = sqrt((pos(:,1) - x).^2 + (pos(:,2) - y).^2 + (pos(:,3) - z).^2);
        if (min(D) > 1)
            pos = [pos; x y z];                         %#ok<AGROW>
            queue(end + 1) = size(pos,1);               %#ok<AGROW>
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