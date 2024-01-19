function [ img ] = BH_movingRMS( img, featureSize )
%Calculate a moving average using circular convolution.
%   Detailed explanation goes here
flgOOM = 0;
if ndims(img) ~= numel(featureSize)
  fprintf('ndims img %d\n',ndims(img));
  fprintf('img size %d\n',size(img));
  fprintf('numel featureSize %d\n',ndims(featureSize));
  fprintf('featureSize %d\n', featureSize);
  error('image must have same dimension as number of dims given in feature size (movingRMS).')
end

padBy = max(ceil(featureSize./2));
totalSize = padBy.*2 + size(img);

% Pad the image by 1/2 the feature size, and the mask to the corresponding
% total size.
padImg = BH_multi_padVal(size(img),totalSize);
padKernel = BH_multi_padVal(featureSize,totalSize);
trimVal = BH_multi_padVal(totalSize,size(img));

if numel(featureSize) == 2
  flg2d = true;
  padBy = padBy.*[1,1];
else
  padBy = padBy.*[1,1,1];
  flg2d = false;
end

if (flg2d)
  rmsFilter = BH_mask3d('sphere',featureSize,featureSize./2,[0,0],'2d');
else
  rmsFilter = BH_mask3d('sphere',featureSize,featureSize./2,[0,0,0]);
end


rmsFilter = BH_padZeros3d(rmsFilter,padKernel(1,:),padKernel(2,:),'GPU','single');
rmsFilter = rmsFilter ./ (sum(rmsFilter(:))-1);

img = padarray(img,padBy,'both','symmetric');


img = BH_padZeros3d(fftshift(real(ifftn(fftn(rmsFilter).*fftn(img.^2)))),...
  trimVal(1,:),trimVal(2,:),'GPU','single');

clearvars -except img



img(img < 1e-9) = max(img(:));
img = sqrt(img);
clearvars -except img
end

