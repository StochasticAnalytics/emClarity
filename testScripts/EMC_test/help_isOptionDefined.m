function ISTHERE = help_isOptionDefined(OPTIONCELL, PARAM)

ISTHERE = iscell(OPTIONCELL) && ~isempty(OPTIONCELL) && any(strcmp(OPTIONCELL(:, 1), PARAM));

end