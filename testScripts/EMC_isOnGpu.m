function ISONGPU = EMC_isOnGpu(NUM)

if isa(NUM, 'gpuArray')
    ISONGPU = true;
else
    ISONGPU = false;
end

end
