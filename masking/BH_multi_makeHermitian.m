function [ fullXform ] = BH_multi_makeHermitian(halfXform, originalSize, scalingFactor)
% Make the hermitian symmetric pair and optionally rescale


% For now just recreate the Hermitian pair
padConv = 0;
switch numel(originalSize)
  case 1
    padConv = 1;
    flg2d = false;
  case 2
    flg2d = true;
    originalSize = [originalSize, 1];
  case 3
    flg2d = false;
end

  
% pull the halfXform to the gpu incase it is not already there
halfXform = gpuArray(halfXform);

if (padConv)
  [nX,nY,nZ] = size(halfXform);
  nX = nX + originalSize(1);

  oX = originalSize(1)+1;
  isOdd = 1;
  fullXform = zeros([nX,nY,nZ],'single','gpuArray');
  
else
  nX = originalSize(1);
  nZ = originalSize(3);
  oX = floor(nX/2)+1;
  isOdd = mod(nX,2);

  fullXform = zeros(originalSize,'single','gpuArray');   
end

if (flg2d)
  fullXform(oX-1+isOdd:end,:,:) = halfXform;
  fullXform(1:oX-2+isOdd,:) = halfXform(oX-1+isOdd:-1:2,:);
else
  fullXform(oX-1+isOdd:end,:,:) = halfXform;
  fullXform(1:oX-2+isOdd,:,:) = halfXform(oX-1+isOdd:-1:2,:,nZ:-1:1);
end

clear halfXform

if (abs(scalingFactor - 1.0) > 1e-3)
  fullXform =  BH_reScale3d(fullXform,'',scalingFactor,'GPU');
end

end
