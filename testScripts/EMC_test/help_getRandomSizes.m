function [SIZE] = help_getRandomSizes(NUMBER, RANGE, DIMENSIONS)

if strcmpi(DIMENSIONS, '3d')
    ndim = 3;
elseif strcmpi(DIMENSIONS, '2d')
    ndim = 2;
elseif strcmpi(DIMENSIONS, '1d')
    ndim = 1;
else
    error('number of dimensions not supported')
end

SIZE = num2cell(floor((RANGE(2)-RANGE(1)) .* rand(NUMBER, ndim) + RANGE(1)), 2);

end
