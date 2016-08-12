function [ e ] = meetingSchedulingGraph(options)
%MEETINGSCHEDULINGGRAPH randomized meeting scheduling graph
%   Detailed explanation goes here
% 
% e = meetingSchedulingGraph(options)

nAgents = getSubOption(uint16(30), 'uint16', options, 'nAgents');
nMeetings = getSubOption(uint16(5), 'uint16', options, 'nMeetings');
meetingSizeFun = getSubOption(@(x) randi(x), 'function_handle', options, 'meetingSizeFun');

e = cell(1,nMeetings);

for m = 1:nMeetings
    meetingSize = max(2, meetingSizeFun(nAgents));
    rp = randperm(nAgents);
    attendingAgents = rp(1:meetingSize);
    
    % These agents will share a hard constraints
    [A, B] = meshgrid(attendingAgents);
    e_meeting = [A(:) B(:)];
    e_meeting = e_meeting(e_meeting(:,1) < e_meeting(:,2),:); %because undirected

    % Add nicely sorted rows
    e{m} = sortrows(e_meeting);
end

end

