function [ IMG ] = BH_bandLimitCenterNormalize_halfGrid( IMG, BANDPASS, ...
                                                         bhFourierTransformer, ...
                                                         PAD )
%Bandlimit, center, normalize and img.

% cufft does root(n) scaling both ways, the prefactor (p) then should be
% p^2 = n*root(n) or p = sqrt(sqrt(n^3)) = n^(3/4)
IMG = bhFourierTransformer.fwdFFT(BH_padZeros3d(IMG - mean(IMG(:)), 'fwd', PAD,'GPU', 'single'),0.75,0);

IMG = IMG .* BANDPASS;

IMG = IMG ./ (sqrt(sum(abs(IMG(:)).^2)));

end

