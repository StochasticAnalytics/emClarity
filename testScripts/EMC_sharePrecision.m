function SHARE = EMC_sharePrecision(NUM1, NUM2)
% Are NUM1 and NUM2 have the same precision?
%

if isa(NUM1, 'gpuArray')
    NUM1_class = classUnderlying(NUM1);
else
    NUM1_class = class(NUM1);
end

if isa(NUM2, 'gpuArray')
    NUM2_class = classUnderlying(NUM2);
else
    NUM2_class = class(NUM2);
end

if strcmp(NUM2_class, NUM1_class)
    SHARE = true;
else
    SHARE = false;
end

end
