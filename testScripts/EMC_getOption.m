function [OPTION] = EMC_getOption(OPTION, ONLY, FILTER)
%
% Try to make parameters with default value easier.
% Check that the OPTION cell|struct has the correct fields.
%
% Input:
%   OPTION (cell|struct): 	If cell: {param1, value1 ; param2, value2 ; ...} and param* should be 
%                                    non-empty character vectors or a string scalars.
%
%   ONLY (cell):           	{'param1', 'param2', ...}
%
%   FILTER (bool):        	If true: remove from OPTION the fields not in ONLY.
%                           If false: raise an error if OPTION contains a field
%                           that is not in ONLY.
%
% Output:
%   OPTION (struct):        Checked/filtered structure.
%
% Created:  18Jan2020
% Version:  v.1.0 unittest (TF, 20Jan2020).
%

if iscell(OPTION)
    if ~isempty(OPTION)
        fields = OPTION(:, 1);
        OPTION = cell2struct(OPTION(:, 2), fields, 1);
        for idx = 1:numel(fields)
            if ~any(strcmp(ONLY, fields{idx}))
                if FILTER
                    OPTION = rmfield(OPTION, fields{idx});
                else
                    error('EMC_getOption:OPTION', ...
                          'OPTION cell contains an unexpected field: %s', fields{idx})
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
                error('EMC_getOption:OPTION', ...
                      'OPTION structure contains an unexpected field: %s', fields{idx})
            end
        end
    end
else
    error('EMC_getOption:OPTION', ...
          'OPTION should be a cell or a structure, got %s', class(OPTION))
end

end
