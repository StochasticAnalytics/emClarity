function IMAGE = EMC_irfftn(DFT, SIZE)
%
% N-D inverse Fast Fourier Transform wrapper, used for 2d/3d non-centered DFT of real images.
% Using the hermitian symmetry, recompute the full real-space image.
%
% Input:
%   DFT (numeric):      2d or 3d non-centered (zero-frequency first)
%                       half discrete Fourier transform (complex).
%
%   SIZE (vector):      Size of the output IMAGE.
%
% Output:
%   IMAGE (numeric):    Full inverse discrete Fourier transform of DFT (real).
%
% Note:
%   - Recomputing the full DFT has a non negligeable overload. For large arrays and
%     when a lot of computation if done on the half grids, it is usually worth it.
%
% Example:
%   - img = rand(10,10);       
%     out1 = EMC_irfftn(EMC_rfftn(img), size(img));
%     out2 = ifftn(fftn(img));
%     any(abs(out2-out1) > 1e-7, 'all')  % out1 == out2
%

% Created:  3Feb2020
% Version:  v.1.0.  unittest (TF, 8Feb2020).
%

[precision, isOnGpu] = EMC_getClass(DFT);
if isOnGpu
    IMAGE = zeros(SIZE, precision, 'gpuArray');
else
    IMAGE = zeros(SIZE, precision);
end

cR = floor(SIZE/2) + 1;  % center receiver
cD = ceil(SIZE/2);  % center donor

switch ndims(DFT)
    case 3
        IMAGE(1:cR(1),     :,           :)           = DFT;
        IMAGE(cR(1)+1:end, 1,           1)           = conj(DFT(cD(1):-1:2, 1,              1));

        IMAGE(cR(1)+1:end, 2:cR(2),     1)           = conj(DFT(cD(1):-1:2, end:-1:cD(2)+1, 1));
        IMAGE(cR(1)+1:end, cR(2)+1:end, 1)           = conj(DFT(cD(1):-1:2, cD(2):-1:2,     1));
        IMAGE(cR(1)+1:end, 1,           2:cR(3))     = conj(DFT(cD(1):-1:2, 1,              end:-1:cD(3)+1));
        IMAGE(cR(1)+1:end, 1,           cR(3)+1:end) = conj(DFT(cD(1):-1:2, 1,              cD(3):-1:2));

        IMAGE(cR(1)+1:end, 2:cR(2),     2:cR(3))     = conj(DFT(cD(1):-1:2, end:-1:cD(2)+1, end:-1:cD(3)+1));
        IMAGE(cR(1)+1:end, 2:cR(2),     cR(3)+1:end) = conj(DFT(cD(1):-1:2, end:-1:cD(2)+1, cD(3):-1:2));
        IMAGE(cR(1)+1:end, cR(2)+1:end, 2:cR(3))     = conj(DFT(cD(1):-1:2, cD(2):-1:2,     end:-1:cD(3)+1));
        IMAGE(cR(1)+1:end, cR(2)+1:end, cR(3)+1:end) = conj(DFT(cD(1):-1:2, cD(2):-1:2,     cD(3):-1:2));
    case 2
        IMAGE(1:cR(1),     :)           = DFT;
        IMAGE(cR(1)+1:end, 1)           = conj(DFT(cD(1):-1:2, 1));
        IMAGE(cR(1)+1:end, 2:cR(2))     = conj(DFT(cD(1):-1:2, end:-1:cD(2)+1));
        IMAGE(cR(1)+1:end, cR(2)+1:end) = conj(DFT(cD(1):-1:2, cD(2):-1:2));
    otherwise
        error('EMC:irfftn', 'DFT should be a 2d or 3d array, got %d', ndims(DFT))
end

IMAGE  = ifftn(IMAGE, 'symmetric');

end  % EMC_irfftn
