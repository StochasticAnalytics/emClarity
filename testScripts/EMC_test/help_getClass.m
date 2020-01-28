function [METHOD, PRECISION] = help_getClass(IMG)

if isa(IMG, 'gpuArray')
    METHOD = 'gpu';
    PRECISION = classUnderlying(IMG);
else
    METHOD = 'cpu';
    PRECISION = class(IMG);
end

end  % help_getClass