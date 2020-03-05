function IMAGE = EMC_convn(IMAGE, KERNEL)
%
% IMAGE = EMC_convn(IMAGE, KERNEL)
% Wrapper of the convn|conv2 functions to convolve 2d/3d IMAGEs with separable kernels.
%
% WARNING: If the KERNEL is a vector, this function thinks it is a separable kernel. This is
%          not equivalent to convn(IMAGE, KERNEL, 'same') in that case.
%
% Input:
%   IMAGE (2d|3d numeric):              2d or 3d numeric array to convolve.
%
%   KERNEL (2d|3d numeric | vector):    Kernel to convolve with the IMAGE.
%                                       If vector: use the kernel as a separable kernel.
%                                       If matrix: use the kernel as standard kernel.
%                                                  This is equivalent to convn(IMAGE, KERNEL, 'same').
%                                       NOTE: To preserve the scaling of the IMAGE, the kernel
%                                             should be correctly normalized. See EMC_gaussianKernel
%                                             for more details.
%
% Output:
%   IMAGE (2d|3d numeric):              Convolved image.
%
% Note:
%   - For the 2d case with separable kernels, convn(convn(img, k), k) is slower than
%     conv2(img, k, k) on CPU. Their are equivalent on GPU.
%
% Example:
%   - a = rand(1000,1000);
%     % b and c are equal, up to the single/double precision accuracy.
%     b = EMC_convn(a, EMC_gaussianKernel([5,5], 1, {}));  % classic kernel
%     c = EMC_convn(a, EMC_gaussianKernel([1,5], 1, {}));  % separable kernel
%
% See also EMC_gaussianKernel
%

% Created:  31Jan2020, R2019a
% Version:  v.1.0.  Allow row and column vectors (TF, 2Feb2020).
%           v.1.1.  Unittest (TF, 3Feb2020).
%

% let convn|conv2 do the checks for IMAGE.
if ~isnumeric(IMAGE) || isscalar(IMAGE)
    error('EMC:IMAGE', 'IMAGE should be a 2d or 3d numeric, got %s of size: %s', ...
          class(IMAGE), mat2str(size(IMAGE)))
elseif ~EMC_sharePrecision(IMAGE, KERNEL) || ~EMC_shareMethod(IMAGE, KERNEL)
    error('EMC:IMAGE', 'IMAGE should have same precision and method than KERNEL')
elseif isrow(KERNEL)
    if ndims(IMAGE) == 3
        IMAGE = convn(convn(convn(IMAGE, KERNEL', 'same'), KERNEL, 'same'), reshape(KERNEL,1,1,[]), 'same');
    else
        IMAGE = conv2(KERNEL, KERNEL, IMAGE, 'same');
    end
elseif iscolumn(KERNEL)
    if ndims(IMAGE) == 3
        IMAGE = convn(convn(convn(IMAGE, KERNEL, 'same'), KERNEL', 'same'), reshape(KERNEL,1,1,[]), 'same');
    else
        IMAGE = conv2(KERNEL, KERNEL, IMAGE, 'same');
    end
else
    IMAGE = convn(IMAGE, KERNEL, 'same');
end

end  % EMC_convn
