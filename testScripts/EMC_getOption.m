function OPTION = EMC_getOption(OPTION, ONLY, FILTER)
%
% OPTION = EMC_getOption(OPTION, ONLY, FILTER)
% Check that the OPTION cell|struct has the correct fields.
%
% Input:
%   OPTION (cell|struct):   If cell: {param1, value1; param2, value2; ...}
%                                    param* should be non-empty character vectors.
%                                    size: [n, 2] with n being the number of parameters.
%                                    Can be empty.
%
%   ONLY (cell):            {param1, param2, ...} and param* should be character vectors.
%                           Can be empty.
%
%   FILTER (bool):          If true: remove from OPTION the fields not in ONLY.
%                           If false: raise an error if OPTION contains a field
%                           that is not in ONLY.
%
% Output:
%   OPTION (struct):        Checked/filtered structure.
%

% Created:  18Jan2020
% Version:  v.1.0   unittest (TF, 20Jan2020).
%           v.1.0.1 new error identifier convention (TF, 30Jan2020).
%           v.1.0.2 clearer error message when the cell has not the correct size.
%

if iscell(OPTION)
    if ~isempty(OPTION)
        if size(OPTION, 2) ~= 2
            error('EMC:OPTION', 'OPTION should be a nx2 cell, got %s cell', mat2str(size(OPTION)))
        end
        fields = OPTION(:, 1);
        OPTION = cell2struct(OPTION(:, 2), fields, 1);
        for idx = 1:numel(fields)
            if ~any(strcmp(ONLY, fields{idx}))
                if FILTER
                    OPTION = rmfield(OPTION, fields{idx});
                else
                    error('EMC:OPTION', 'OPTION cell contains an unexpected field: %s', fields{idx})
                end
            end
        end
    else
        OPTION = struct;
    end
elseif isstruct(OPTION)
    fields = fieldnames(OPTION);
    for idx = 1:numel(fields)
        if ~any(strcmp(ONLY, fields{idx}))
            if FILTER
                OPTION = rmfield(OPTION, fields{idx});
            else
                error('EMC:OPTION', 'OPTION structure contains an unexpected field: %s', fields{idx})
            end
        end
    end
else
    error('EMC:OPTION', 'OPTION should be a cell or a structure, got %s', class(OPTION))
end

end  % EMC_getOption
