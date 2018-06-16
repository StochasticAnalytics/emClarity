function [ IMG ] = BH_bandLimitCenterNormalize_cpu( IMG, BANDPASS, MASK, ...
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


if isnumeric(MASK) || islogical(MASK)
  IMG(MASK) = IMG(MASK) - mean(double(IMG(MASK)));
 % IMG(MASK) = IMG(MASK) ./rms(double(IMG(MASK)));
end


IMG = fftn(BH_padZeros3d(IMG, PAD(1,:), PAD(2,:),'cpu', PRECISION));

if isnumeric(BANDPASS)
  IMG = IMG .* BANDPASS;
end

  % set rms to 1 extra parenthesis intentional to save numel(tmpIMG) divisions
%%%if ~(isnumeric(MASK) || islogical(MASK))
  IMG = IMG ./ ((sqrt(sum(sum(sum(abs(IMG).^2)))) ./ numel(IMG)));
%%%end
clear  BANDPASS MASK PAD PRECISION
end

