function [ IMG ] = BH_bandLimitCenterNormalize_2( IMG, BANDPASS, MASK, ...
                                                                 PAD, PRECISION )
%Bandlimit, mask, center, normalize and img.
%   Masks ROI in image, then bandlimits while also setting the mean to 0 and the
%   variance. Returns FT, if image is to be recovered must take ifft.
%
%   Input:
%
%   IMG = 3d image to operate on.
%
%   BANDPASS = bandpass filter to apply, obviously must be the same dimensions
%   as IMG.
%
%   MASK = binary mask for mean subtraction
%
%   PAD  = zero padding
%   
%   PRECISION = single, double, singleTaper, doubleTaper
%
%   Note: F(0,0,0) = prod(size(img)).*mean(img(:)) st we might set the img
%   to zero mean by manipulating the fft, however, the padding around the
%   image, or the masked out area should equal the mean, which must be
%   addressed prior to taking the fft.
%
%   use double precision to calculate the mean, and after masking.

applyBandpass = 0;
applyShift = 0;
if isnumeric(BANDPASS)
  applyBandpass = 1;
elseif strcmpi(BANDPASS, 'shift')
  % use the mask to shift the image prior to fft as needed for fourier
  % interp.
  applyShift = 1;
end
flgMeanWOMask = 1;
% if isnumeric(MASK) || islogical(MASK)
%   IMG(MASK) = IMG(MASK) - mean(double(IMG(MASK)));
%   flgMeanWOMask = 0;
% end

% When looping through a real stack and converting in place to fft, the 2nd
% and on will be complex, which is not convertable to a logical in padZeros
if ( applyShift )
  IMG = BH_padZeros3d(IMG, PAD(1,:), PAD(2,:),'GPU', PRECISION,real(mean(IMG(:))));
  IMG = fftn(IMG(MASK));
else
  IMG = fftn(BH_padZeros3d(IMG, PAD(1,:), PAD(2,:),'GPU', PRECISION,real(mean(IMG(:)))));
end
if (flgMeanWOMask)
  IMG(1) = 0;
end

if (applyBandpass)
  IMG = IMG .* BANDPASS;
end
  % set rms to 1 extra parenthesis intentional to save numel(tmpIMG) divisions
%%%if ~(isnumeric(MASK) || islogical(MASK))
  IMG = IMG ./ ((sqrt(sum(sum(sum(abs(IMG).^2)))) ./ numel(IMG)));
%%%end
clear  BANDPASS MASK PAD PRECISION
end

