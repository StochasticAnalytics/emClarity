function VALUE = help_getOptionParam(OPTIONCELL, PARAM)

VALUE = OPTIONCELL{strcmp(OPTIONCELL(:, 1), PARAM), 2};

end