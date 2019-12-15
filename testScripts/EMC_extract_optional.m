function [STRUCTURE] = EMC_extract_optional(OPTIONAL, ONLY)
%
% Try to make parameters with default value easier.
%
% OPTIONAL (cell|struct):   If cell: {field1,value1 ; field2,value2 ; ...}
%                           If struct: return
%
% ONLY (cell):              {'field1', 'field2', etc.}
%                           Raise an error if OPTIONAL contains a field
%                           not in this cell.

if iscell(OPTIONAL)
    if ~isempty(OPTIONAL)
        fields = OPTIONAL(:, 1);
        STRUCTURE = cell2struct(OPTIONAL(:, 2), fields, 1);
        for idx = 1:numel(fields)
            if ~contains(ONLY, fields{idx})
                error('cell contains an unexpected field: %s', fields{idx})
            end
        end
    else
        STRUCTURE = struct();
    end
elseif isstruct(OPTIONAL)
    fields = fieldnames(OPTIONAL);
    for idx = 1:numel(fields)
        if ~contains(ONLY, fields{idx})
                error('structure contains an unexpected field: %s', fields{idx})
        end
    end
    STRUCTURE = OPTIONAL;
else
    error('OPTIONAL should be a cell or a structure, got %s', class(OPTIONAL))
end
end