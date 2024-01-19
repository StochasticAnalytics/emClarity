function [parameter_struct] = EMC_assert_deprecated_substitution(parameter_struct, default_value, current_field, deprecated_field)

    % Copying this stuct around is probaby not the most efficient way to do this
    % but accuracy is more important than speed here.

    % Handle type checks and default settings outside.
    if isfield(parameter_struct, deprecated_field)
        parameter_struct.(current_field) = parameter_struct.(deprecated_field);
        parameter_struct = rmfield(parameter_struct, deprecated_field);
    else
        if ~isfield(parameter_struct, current_field)
            parameter_struct.(current_field) = default_value;
        end
    end

end