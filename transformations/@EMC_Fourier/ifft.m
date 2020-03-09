function IMAGE = ifft(obj, DFT)
%
% IMAGE = obj.ifft(DFT)
% N-D inverse Fast Fourier Transform wrapper, used for 2d/3d DFT of real images.
%
% Input:
%   DFT (numeric):      2d or 3d Discrete Fourier Transform (complex).
%
% Output:
%   IMAGE (numeric):    Inverse Discrete Fourier Transform of DFT (real).
%
% Property used:
%   obj.half
%   obj.half_wrap
%   obj.size_real
%   obj.centered
%   obj.index_ifftshift
%   obj.index_half2full
%
% Method used:
%   EMC_Fourier.half2full
%
% Note:
%   - Recomputing the full DFT has a non negligeable overload. This is only worth it
%     if a lot of operation is done on the half grid whether than the full grid.
%
% Example:
%   - >> ft = EMC_Fourier([128,128,128], 'gpu', {});
%     >> img = rand(128,128,128);
%     >> dft = ft.fft(img);
%     >> img_back = ft.ifft(dft);  % img_back == img
%

% Created:  3Feb2020
% Version:  v.1.0.  unittest (TF, 8Feb2020).
%           v.1.1.  use EMC_maskIndex and a persistent wrapper to speed up
%                   Fourier coefficients wrapping (half to full grid) (TF, 14Feb2020).
%           v.1.2.  EMC_irfftn is integrated into EMC_Fourier.ifft; can used
%                   half/full not-centered/centered DFTs (TF, 8Mar2020).
%

if obj.half
    IMAGE  = ifftn(EMC_Fourier.half2full(obj.half_wrap, DFT, obj.size_real, obj.index_half2full), 'symmetric');
elseif obj.centered
    IMAGE  = ifftn(DFT(obj.index_ifftshift), 'symmetric');
else
    IMAGE  = ifftn(DFT, 'symmetric');
end

end  % ifft
