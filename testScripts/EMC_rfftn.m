function DFT = EMC_rfftn(IMAGE)
%
% N-D Fast Fourier Transform of a real 2d/3d real IMAGE.
% Use the hermitian symmetry of the IMAGE and return the non-redundant coefficients only.
%
% Input:
%   IMAGE (numeric):    2d or 3d numerical array (real).
%
% Output:
%   DFT (numeric):      Half discrete Fourier transform of IMAGE (complex).
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

% This is a real pain that MATLAB doesn't support this; compute the full spectrum
% and truncate to return non-redundant coefficients only.
DFT = fftn(IMAGE);
if ndims(DFT) == 3
    DFT = DFT(1:(floor(size(DFT, 1)/2) + 1), :, :);
else
    DFT = DFT(1:(floor(size(DFT, 1)/2) + 1), :);
end

end  % EMC_rfftn
