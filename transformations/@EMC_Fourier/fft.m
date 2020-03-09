function DFT = fft(obj, IMAGE)
%
% DFT = obj.fft(IMAGE)
% N-D Fast Fourier Transform of a real 2d/3d IMAGE.
%
% Input:
%   IMAGE (numeric):    2d or 3d numerical array (real).
%
% Output:
%   DFT (numeric):      Half or full discrete Fourier transform of IMAGE (complex).
%
% Property used:
%   obj.half
%   obj.centered
%   obj.size_freq
%   obj.index_fftshift
%
% Example:
%   - >> ft = EMC_Fourier([128,128,128], 'gpu', {});
%     >> img = rand(128,128,128);
%     >> dft = ft.fft(img);
%

% Created:  3Feb2020
% Version:  v.1.0.  unittest (TF, 8Feb2020).
%           v.1.1.  EMC_rfftn is integrated into EMC_Fourier.fft; can return
%                   half/full not-centered/centered DFTs (TF, 8Mar2020).
%

% This is a real pain that MATLAB doesn't support this by default.
DFT = fftn(IMAGE);

% Return only non-redundant spectrum.
if obj.half
    DFT = DFT(1:obj.size_freq(1), :, :);
end

% Shift the zero frequency to the center of the spectrum.
if obj.centered
    DFT = DFT(obj.index_fftshift);
end

end  % fft
