function [PRECISION, flgGPU, METHOD] = EMC_getClass(IMAGE)
%
% Get the METHOD and PRECISION of a given IMAGE
%

if isa(IMAGE, 'gpuArray')
    PRECISION = classUnderlying(IMAGE);
    flgGPU = true;
    METHOD = 'gpu';
else
    PRECISION = class(IMAGE);
    flgGPU = false;
    METHOD = 'cpu';
end

end
