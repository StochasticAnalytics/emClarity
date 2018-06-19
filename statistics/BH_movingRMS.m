function [ RMS ] = BH_movingRMS( img, featureSize )
%Calculate a moving average using circular convolution.
%   Detailed explanation goes here

if ndims(img) ~= numel(featureSize)
  error('image must have same dimension as number of dims given in feature size.')
end

if isa(img,'gpuArray')
  useGPU = 1;
else
  useGPU = 0;
end

[d1,d2,d3] = size(img);

% Pad the image by 1/2 the feature size, and the mask to the corresponding
% total size.

padImg = [floor(featureSize./2);ceil(featureSize./2)];
padAvg = [floor(size(img)./2);ceil(size(img)./2)];
totalSize = size(img)+sum(padImg,1);

% use the same type as the image
if isa(gather(img(1)), 'single')
  precision = 'single';
elseif isa(gather(img(1)), 'double')
  precision = 'double';
else
  error('Img must be float or double!\n')
end

  if numel(featureSize) ==  2
    rmsFilterCore = BH_mask3d('sphere',featureSize,featureSize./2,[0,0],'2d');
  else
    rmsFilterCore = BH_mask3d('sphere',featureSize,featureSize./2,[0,0,0]);
  end
  
if (useGPU)
  rmsFilter = zeros(totalSize, precision, 'gpuArray');
  paddedImg = zeros(totalSize, precision,'gpuArray');
else
  rmsFilter = zeros(totalSize, precision);
  rmsFilterCore = gather(rmsFilterCore);
  paddedImg = zeros(totalSize, precision);
end

img = img.^2;

if numel(featureSize) == 2
  
  rmsFilter(padAvg(1,1)+1:end-padAvg(2,1), ...
             padAvg(1,2)+1:end-padAvg(2,2)) = rmsFilterCore;

%   rmsFilter(padAvg(1,1)+1:end-padAvg(2,1), ...
%              padAvg(1,2)+1:end-padAvg(2,2)) = ...
%                                 zeros(featureSize) + (prod(featureSize)-1)^-1;
  rmsFilter = rmsFilter ./ (sum(rmsFilter(:))-1);

  paddedImg(padImg(1,1)+1:end-padImg(2,1), ...
            padImg(1,2)+1:end-padImg(2,2)) = img;
          
  % pad the image like 'replicate' corners first        
  paddedImg(1:padImg(1,1),1:padImg(2,1)) = ...
        img(1:padImg(1,1),1:padImg(2,1));
  paddedImg(end-padImg(1,1)+1:end,1:padImg(2,1)) = ...
        img(end-padImg(1,1)+1:end,1:padImg(2,1));
  paddedImg(1:padImg(1,1),end-padImg(2,1)+1:end) = ...
        img(1:padImg(1,1),end-padImg(2,1)+1:end);
  paddedImg(end-padImg(1,1)+1:end,end-padImg(2,1)+1:end) = ...
        img(end-padImg(1,1)+1:end,end-padImg(2,1)+1:end);
      
  paddedImg(1:padImg(1,1),padImg(2,1)+1:end-padImg(2,1)) = ...
        img(1:padImg(1,1),:); 
  paddedImg(end-padImg(1,1)+1:end,padImg(2,1)+1:end-padImg(2,1)) = ...
        img(end-padImg(1,1)+1:end,:);
  paddedImg(1+padImg(1,1):end-padImg(1,1),1:padImg(2,1)) = ...
        img(:,1:padImg(2,1));
  paddedImg(1+padImg(1,1):end-padImg(1,1),end-padImg(2,1)+1:end) = ...
        img(:,end-padImg(2,1)+1:end);
      
  
  RMS = fftshift(real(ifftn(fftn(rmsFilter).*fftn(paddedImg))));
  RMS = RMS(padImg(1,1)+1:end-padImg(2,1), ...
            padImg(1,2)+1:end-padImg(2,2));
          
else
  
    rmsFilter(padAvg(1,1)+1:end-padAvg(2,1), ...
             padAvg(1,2)+1:end-padAvg(2,2), ...
             padAvg(1,3)+1:end-padAvg(2,3)) = rmsFilterCore;
%   rmsFilter(padAvg(1,1)+1:end-padAvg(2,1), ...
%              padAvg(1,2)+1:end-padAvg(2,2), ...
%              padAvg(1,3)+1:end-padAvg(2,3)) = ...
%                                 zeros(featureSize) + prod(featureSize)^-1;
  rmsFilter = rmsFilter ./ (sum(rmsFilter(:))-1);

  paddedImg(padImg(1,1)+1:end-padImg(2,1), ...
            padImg(1,2)+1:end-padImg(2,2), ...
            padImg(1,3)+1:end-padImg(2,3)) = img;
          
  RMS = fftshift(real(ifftn(fftn(rmsFilter).*fftn(paddedImg))));
  RMS = RMS(padImg(1,1)+1:end-padImg(2,1), ...
            padImg(1,2)+1:end-padImg(2,2), ...
            padImg(1,3)+1:end-padImg(2,3));
end                                

% deal with any negative numbers at the edges, or zeros which should only
% occur if there is a homogenous field due to padding by a constant.
RMS(RMS < 1e-9) = max(RMS(:));
RMS = sqrt(RMS);
clearvars -except RMS
end

