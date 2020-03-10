function IMAGE = EMC_irfftn(DFT, SIZE)
%
% IMAGE = EMC_irfftn(DFT, SIZE)
% N-D inverse Fast Fourier Transform wrapper, used for 2d/3d not-centered DFT of real images.
% Using the hermitian symmetry, recompute the full real-space image.
%
% Input:
%   DFT (numeric):      2d or 3d not-centered (zero-frequency first)
%                       half discrete Fourier transform (complex).
%
%   SIZE (vector):      Size of the output IMAGE.
%
% Output:
%   IMAGE (numeric):    Full inverse discrete Fourier transform of DFT (real).
%
% Note:
%   - Recomputing the full DFT has a non negligeable overload. This is only worth it
%     if a lot of operation is done on the half grid whether than the full grid.
%
%   - This function allocate in the heap the wrapper used to reconstruct the full DFT from the
%     half DFT (the input). As such if this function is called multiple time and if the SIZE does
%     not change, the same wrapper is used, speeding up the execution. If the SIZE changes, the
%     function automatically computes a new wrapper. To free the memory allocated by the wrapper,
%     use 'clear EMC_irfttn'.
%
% Example:
%   - img = rand(128,128);
%     out1 = EMC_irfftn(EMC_rfftn(img), size(img));
%     out2 = ifftn(fftn(img));
%     any(abs(out2-out1) > 1e-7, 'all')  % out1 == out2
%
% Other EMC-files required:
%   EMC_getClass, EMC_maskIndex
%
% See also EMC_maskIndex
%

% Created:  3Feb2020
% Version:  v.1.0.  unittest (TF, 8Feb2020).
%           v.1.1.  use EMC_maskIndex and a persistent wrapper to speed up
%                   Fourier coefficients wrapping (half to full grid) (TF, 14Feb2020).
%

[precision, isOnGpu, method] = EMC_getClass(DFT);
if isOnGpu
    IMAGE = zeros(SIZE, precision, 'gpuArray');
else
    IMAGE = zeros(SIZE, precision);
end

% h2f wrapper
persistent wrap
if isempty(wrap) || ~isequal(size(wrap), SIZE)
    wrap = EMC_maskIndex('nc2nc', SIZE, method, {});
end

cR = floor(SIZE/2) + 1;  % center receiver

switch ndims(DFT)
    case 3
        IMAGE(1:cR(1), :, :) = DFT;
        IMAGE = IMAGE(wrap);
        IMAGE(cR(1)+1:end, :, :) = conj(IMAGE(cR(1)+1:end, :, :));

    case 2
        IMAGE(1:cR(1),     :) = DFT;
        IMAGE = IMAGE(wrap);
        IMAGE(cR(1)+1:end, :) = conj(IMAGE(cR(1)+1:end, :));

    otherwise
        error('EMC:irfftn', 'DFT should be a 2d or 3d array, got %d', ndims(DFT))
end

IMAGE  = ifftn(IMAGE, 'symmetric');

end  % EMC_irfftn
