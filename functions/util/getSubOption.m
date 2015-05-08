%% GETSUBOPTION
% *Get an option from a nested options struct*
%
% Get an option from a nested struct of multiple options. This function
% will check if the option exists, and if not will return the default value
% (if provided).
%
%% Usage:
% value = getSubOption(default, options, fieldname)
%   Check the options struct for the option with the field name. If it
%   doesn't exist, return the default value.
%
% value = getSubOption(default, type, options, fieldname)
%   Check for the option, and make sure it is of the correct type.
%
% value = getSubOption(default[, type], options, field, name, value ...)
%   Check for the option in the field options.field.name.value and return
%   the default value if it doesn't exist. If a type is provided, check if
%   it is of correct type.
%
%% Copyright
% * *2014 - TNO*
% * *Author*: Coen van Leeuwen
% * *Since*: August 27, 2014

%% Function Definition
function value = getSubOption(default, options, varargin)

% Option handling of THIS function
% The following exception will be thrown in case of incorrect usage of ME
usage = MException('GETSUBOPTION:USAGE', ...
    'Error using %s:\nUsage: value = %s(default, type, options, fieldnames...)', ...
    mfilename, mfilename);

% Minimum number of arguments
if nargin < 3
    usage.throwAsCaller();
end

% If the 'type' argument is provided, it should be the second argument
if ischar(options)
    type = options;
    options = varargin{1};
    fields = varargin(2:end);
else
    type = '';
    fields = varargin;
end

% Check is the options is a struct, and all fieldnames are chars
if ~isstruct(options) || ~all(cellfun(@ischar, fields))
    usage.throwAsCaller();
end

% The option name that the user wants to check
name = ['options' sprintf('.%s', fields{:})];

% Recursively get the field that the user wants to check
value = getSubStructValue(options, fields);

% Set the default value if it is not found
if isempty(value)
    fprintf('Using default value for %s\n', name);
    value = default;
end

% Optionally check the value type
if ~isempty(type)
    assert(isa(value, type), 'Option %s of invalid type %s (%s expected)', ...
        name, class(value), type);
end

end

%% GETSUBSTRUCTVALUE
% *Iteratively get the struct sub value*
%
% Returns the value of the struct.field1.field2.field3 or empty if it
% doesn't exist
function value = getSubStructValue(options, fields)

if ~isfield(options, fields{1})
    value = [];
else
    if numel(fields) == 1
        value = options.(fields{1});
    else
        value = getSubStructValue(options.(fields{1}), fields(2:end));
    end  
end

end
