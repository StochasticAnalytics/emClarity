function [ img, kurtosis ] = BH_movingRMS_2( img, featureSize, varargin )
%Calculate a moving average using circular convolution.
%   Detailed explanation goes here

if numel(featureSize) > 1
  error('movingAverage 2 works only on an isotropic sized image');
end

do_kurtosis = false;
if (nargin > 2)
  do_kurtosis = true;
end

padBy = max(ceil(featureSize./2));

startingSize = size(img);
img = padarray(img,padBy,'both','symmetric');

if (do_kurtosis)
  kurtosis = img;
else
  kurtosis = '';
end

rmsFilter =  EMC_gaussianKernel([1,featureSize],3,'gpu', {});
rmsFilter = rmsFilter ./ (sum(rmsFilter(:)));

img = EMC_convn(single(img.^2), rmsFilter);

[ LIMITS ] = EMC_limits(size(img), startingSize, {});

if (do_kurtosis)
  kurtosis = EMC_convn(single(kurtosis.^4), rmsFilter) ./ img.^2;
  [ kurtosis ] = EMC_resize(kurtosis, LIMITS, {'taper',false});
end

[ img ] = EMC_resize(img, LIMITS, {'taper',false});

clear LIMITS meanFilter padBy startingSize

img(img < 1e-9) = max(img(:));
img = sqrt(img);

end

