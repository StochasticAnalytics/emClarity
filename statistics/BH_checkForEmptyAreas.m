function [ stackOUT ] = BH_checkForEmptyAreas(stackIN)
%Check for regions with low local RMS due to prior resampling
%   Detailed explanation goes here

[d1,d2,d3] = size(stackIN);


  stackOUT = zeros([d1,d2,d3],'uint8');


for iPrj = 1:d3

  iProjection = gpuArray(stackIN(:,:,iPrj));


  localMean = BH_movingAverage(iProjection,[16,16]);
  localRMS  = BH_movingRMS(iProjection-localMean,[16,16]);
  localMean = BH_movingAverage(localRMS > 0.1.*std(localRMS(:)),[16,16]);

  if (abs(1-mean(localMean(:))) > 0.01)
    stackOUT(:,:,iPrj) = uint8(gather(localMean > mean(localMean(:)) - 3*std(localMean(:))));
  else
    stackOUT(:,:,iPrj) = uint8(gather(localMean));
  end
  
end

clear iProjection localMean localRMS 

