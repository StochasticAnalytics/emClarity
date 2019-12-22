function [OPTION] = EMC_extract_option(OPTION, ONLY, FILTER)
%
% Try to make parameters with default value easier.
% Check that the OPTION cell|struct has the correct fields.
%
% OPTION (cell|struct):     If cell: {field1,value1 ; field2,value2 ; ...}
%
% ONLY (cell):              {'field1', 'field2', etc.}
%
% FILTER (bool):            If true: remove from OPTION the fields not in ONLY.
%                           If false: raise an error if OPTION contains a field
%                           that is not in ONLY.
%
%-------
% RETURN:                   Checked/filtered structure.

if iscell(OPTION)
    if ~isempty(OPTION)
        fields = OPTION(:, 1);
        OPTION = cell2struct(OPTION(:, 2), fields, 1);
        for idx = 1:numel(fields)
            if ~contains(ONLY, fields{idx})
                if (FILTER)
                    OPTION = rmfield(OPTION, fields{idx});
                else
                    error('cell contains an unexpected field: %s', fields{idx})
                end
            end
        end
    else
        OPTION = struct();
    end
elseif isstruct(OPTION)
    fields = fieldnames(OPTION);
    for idx = 1:numel(fields)
        if ~contains(ONLY, fields{idx})
            if (FILTER)
                OPTION = rmfield(OPTION, fields{idx});
            else
                error('structure contains an unexpected field: %s', fields{idx})
            end
        end
    end
else
    error('OPTION should be a cell or a structure, got %s', class(OPTION))
end

end
