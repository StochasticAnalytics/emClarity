function NUM = EMC_setMethod(NUM, METHOD)
%
% Push NUM to device or gather NUM to host.
%

if strcmpi(METHOD, 'gpu')
    if ~isa(NUM, 'gpuArray')
        NUM = gpuArray(NUM);
    end
elseif strcmpi(METHOD, 'cpu')
    if isa(NUM, 'gpuArray')
        NUM = gather(NUM);
    end
else
    error('EMC:METHOD', "METHOD should be 'cpu' or 'gpu'")
end

end
