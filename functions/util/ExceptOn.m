function ExceptOn(condition,exceptionID,message,varargin)
%% ExceptOn 
%  Checks if the given condition is true.  If it is, throws an MException.
%
%% Description
%  ExceptOn checks if condition is true.  
% If true, this function throws an exception with indicated exceptionID
% and message.  The exception is generated with the address information of
% the calling routine.
% Both exceptionID and message are required. varargin are arguments to be
% used in the message.  ExceptOn is equivalent to the behavior of the 
% following code:
%
% <html>
%  <b><br />
%  if (condition)<br />
%     SNVE = MException(exceptionID, message, varargin);<br />
%     throw(SNVE);<br/>
%  end<br />
% </b>
% </html>
%
%%% Example: 
%  Checks if st is a structure.  If not, throw an exception with id
% Module:Desc, and message "The input 'st' must be a pure struct.".
% ExceptOn( ~isstruct(st), ...
%              'Module:Desc', ...
%              'The input ''st'' must be a %s struct. ', ...
%              'pure');
%
%  Copyright 2011 - TNO
%  Author: Julio Oliveira

if nargin < 3
    error('ExceptOn:BadInput','Too few arguments.');
end

if (~islogical(condition)) 
    error('ExceptOn:BadInput','Input argument ''condition'' must be logical.');
end

if (~ischar(exceptionID)) 
    error('ExceptOn:BadInput','Input argument ''exceptionID'' must be string.');
end

if (~ischar(message)) 
    error('ExceptOn:BadInput','Input argument ''message'' must be string.');
end

if (condition)
     SNVE = MException(exceptionID, ...
         message, ...
         varargin{:});
     throwAsCaller(SNVE);
end

end