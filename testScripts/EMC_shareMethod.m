function SHARE = EMC_shareMethod(NUM1, NUM2)
% Are NUM1 and NUM2 sharing the same method?
%

if isa(NUM1, 'gpuArray')
    if isa(NUM2, 'gpuArray')
        SHARE = true;
    else
        SHARE = false;
    end
else
    if isa(NUM2, 'gpuArray')
        SHARE = false;
    else
        SHARE = true;
    end
end

end