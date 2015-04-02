function [ e ] = realGraph(filename)
%SCALEFREEGRAPH Summary of this function goes here
%   Detailed explanation goes here

if ~exist(filename, 'file')
    filename2 = fullfile(fileparts(mfilename('fullpath')), 'realdata', [filename '.adj']);
    if ~exist(filename2, 'file')
        error('Could not find file %s', filename)
    else
        filename = filename2;
    end
end

fid = fopen(filename);
X = textscan(fid, '%d %d %d');
fclose(fid);

e = [X{1} X{2}];
e = double(e);

e = e - min(min(e)) + 1; %To make sure we go from 1 to n;

end

