function [ img, flgOOM ] = BH_movingAverage_2( img, featureSize )
%Calculate a moving average using circular convolution.
%   Detailed explanation goes here
flgOOM = 0;

if numel(featureSize) > 1
  error('movingAverage 2 works only on an isotropic sized image');
end

padBy = max(ceil(featureSize./2));

startingSize = size(img);
img = padarray(img,padBy,'both','symmetric');

meanFilter =  EMC_gaussianKernel([1,featureSize],3,'gpu', {});
meanFilter = meanFilter ./ sum(meanFilter(:));

img = EMC_convn(single(img), meanFilter);

[ LIMITS ] = EMC_limits(size(img), startingSize, {});
[ img ] = EMC_resize(img, LIMITS, {'taper',false});

clear LIMITS meanFilter padBy startingSize

end

