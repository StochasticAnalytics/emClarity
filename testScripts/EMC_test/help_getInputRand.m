function [IMG] = help_getInputRand(METHOD, PRECISION, SIZE)

if strcmpi(METHOD, 'gpu')
    IMG = rand(SIZE, PRECISION, 'gpuArray');
else
    IMG = rand(SIZE, PRECISION);
end

end  % help_getInputRand