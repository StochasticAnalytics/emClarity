function IMAGE = EMC_medianFilter(IMAGE)
%
% IMAGE = EMC_medianFilter(IMAGE)
% Wrapper of medfilt2 and medfilt3 with default 3x3 or 3x3x3 kernel
%
% Input:
%   IMAGE (numeric):    2d or 3d numeric array
%
% Output:
%   IMAGE (numeric):    2d or 3d filtered image with the same precision and method.
%
% Note:
%   - As medfilt3 requires a lot of memory if IMAGE is on the GPU, it is gathered
%     to the host and later pushed back.

% Created:  23Feb2020, R2019a
% Version:  v.1.0.
%

if ndims(IMAGE) == 3
    if isa(IMAGE, 'gpuArray')
        IMAGE = gpuArray(medfilt3(gather(IMAGE)));
    else
        IMAGE = medfilt3(IMAGE);
    end
else
    IMAGE = medfilt2(IMAGE);
end

end  % EMC_medianFilter
