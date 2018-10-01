function [ defaults ] = override_valid_fields( defaults, inputs, allowed_values)
%OVERRIDE_VALID_FIELDS replaces those fields of 'defaults' for which a valid
%value exists in 'inputs'. The value of each field is validated against a 
%list specified with the same field name in 'allowed_values' (if that field
%exists and is not empty)
%
%[ out ] = override_valid_fields( defaults, inputs, allowed_values)
%
%Based on original code by Mike Brookes
%
%Alastair Moore, November 2013

narginchk(2,3);

if nargin < 3
    allowed_values = [];
end
if ~isstruct(inputs)
    warning('override_valid_fields: Inputs is not a struct so all defaults retained')
    return
end

fn = fieldnames(inputs); %these are the fields we will try to copy
for ii=1:length(fn)
    if isfield(defaults,fn{ii}) %field must exist in defaults
        if ~isfield(allowed_values,fn{ii})...  %doesn't need to exist in allowed values
                || isempty(allowed_values.(fn{ii}))... %allowed values doesn't need to be specified
                || ismember(inputs.(fn{ii}),allowed_values.(fn{ii})) %but if it is, the input value must be allowed
            defaults.(fn{ii}) = inputs.(fn{ii}); %override the value
        else
            display(['Allowed values of '''  fn{ii} ''' are:'])
            display(allowed_values.(fn{ii}))
            error('Specified value is invalid. See above for allowed values.');
        end
    else
        warning([fn{ii} ' is unused']);
    end
end

