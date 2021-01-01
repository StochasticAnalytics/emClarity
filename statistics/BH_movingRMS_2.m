function [ img ] = BH_movingRMS_2( img, featureSize )
%Calculate a moving average using circular convolution.
%   Detailed explanation goes here

if numel(featureSize) > 1
  error('movingAverage 2 works only on an isotropic sized image');
end


padBy = max(ceil(featureSize./2));

startingSize = size(img);
img = padarray(img,padBy,'both','symmetric');

rmsFilter =  EMC_gaussianKernel([1,featureSize],3,'gpu', {});
rmsFilter = rmsFilter ./ (sum(rmsFilter(:)));

img = EMC_convn(single(img.^2), rmsFilter);

[ LIMITS ] = EMC_limits(size(img), startingSize, {});
[ img ] = EMC_resize(img, LIMITS, {'taper',false});

clear LIMITS meanFilter padBy startingSize

img(img < 1e-9) = max(img(:));
img = sqrt(img);

end

