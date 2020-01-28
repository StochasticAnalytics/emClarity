function [IMG] = help_getInputOnes(METHOD, PRECISION, SIZE)

if strcmpi(METHOD, 'gpu')
    IMG = ones(SIZE, PRECISION, 'gpuArray');
else
    IMG = ones(SIZE, PRECISION);
end

end  % help_getInputOnes