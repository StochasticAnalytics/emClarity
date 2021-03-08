function [ img ] = BH_movingRMS_3( img, featureSize, movingAverage )
%Calculate a moving average using circular convolution.
%   Detailed explanation goes here

if numel(featureSize) > 1
  error('movingAverage 2 works only on an isotropic sized image');
end



padBy = max(ceil(featureSize./2));

rmsFilter =  EMC_gaussianKernel([1,featureSize],featureSize./3,'gpu', {});
rmsFilter = rmsFilter ./ (sum(rmsFilter(:)));
mf = max(rmsFilter(:))
difference_term = mf.*(2.*img.*movingAverage + movingAverage.^2);

startingSize = size(img);
img = padarray(img,padBy,'both','symmetric');
img = EMC_convn(single(img.^2), rmsFilter);

[ LIMITS ] = EMC_limits(size(img), startingSize, {});


[ img ] = EMC_resize(img, LIMITS, {'taper',false});

img = img + difference_term;

clear LIMITS meanFilter padBy startingSize

img(img < 1e-9) = max(img(:));
img = sqrt(img);
% m = mean(img(:)) - std(img(:));
% img(img < m) = m;

end

