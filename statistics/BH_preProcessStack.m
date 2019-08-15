function [ stackIN, samplingMask ] = BH_preProcessStack( stackIN, gradientAliasMask,varargin)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

[d1,d2,d3] = size(stackIN);

isGPU = isa(stackIN(1),'gpuArray');

doMedianFilter = 0;
if nargin > 2
  doMedianFilter = varargin{1};
  % When median filtering, do gradient, median, bandpass
  bandPassFilter = gradientAliasMask{2};
  gradientAliasMask = gradientAliasMask{1};
end


% Optionally take a preconstructed filter
if numel(gradientAliasMask) == 1
  if (gradientAliasMask)
    gradientAliasMask = BH_bandpass3d(1.*[d1,d2,1],0,0,-0.235,'GPU','nyquistHigh');
  else
    gradientAliasMask = BH_bandpass3d(1.*[d1,d2,1],0,0,0,'GPU','nyquistHigh');
  end
end

[ samplingMask ] = BH_checkForEmptyAreas(stackIN);

for iPrj = 1:d3

  iProjection = gpuArray(stackIN(:,:,iPrj));
  iProjection = real(ifftn(fftn(iProjection).*gradientAliasMask));
  iProjection(samplingMask(:,:,iPrj)<0.9) = mean(mean(iProjection(samplingMask(:,:,iPrj)>=0.9)));
  
  try
    iProjection = BH_padZeros3d(iProjection,[0,0],[0,0],'GPU','singleTaper',mean(iProjection(:)));
  catch
    figure, imshow3D(stackIN(:,:,iPrj));
    iPrj
    figure, imshow3D(gather(real(ifftn(fftn(iProjection).*gradientAliasMask))));
    figure, imshow3D(gather(samplingMask(:,:,iPrj)));
    error('sdf')
  end
  
 
  largeOutliersMean= mean(iProjection(:));
  largeOutliersSTD = std(iProjection(:));
  largeOutliersIDX = (abs(iProjection(:)) > largeOutliersMean + 5*largeOutliersSTD);
  iProjection(largeOutliersIDX) = (2*largeOutliersSTD).*randn([gather(sum(largeOutliersIDX(:))),1],'single','gpuArray');
 
  if doMedianFilter
    iProjection = medfilt2(iProjection,[1,1].*doMedianFilter);
    iProjection = real(ifftn(fftn(iProjection).*bandPassFilter));
  end
  
  if isGPU
    stackIN(:,:,iPrj) = iProjection - mean(iProjection(:));
  else
    stackIN(:,:,iPrj) = gather(iProjection - mean(iProjection(:)));
  end
  
  
  

end

