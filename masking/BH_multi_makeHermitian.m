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
  nY = originalSize(2);
  nZ = originalSize(3);
  oX = floor(nX/2)+1;
  isOdd = mod(nX,2);

  fullXform = zeros(originalSize,'single','gpuArray');   
end

hermitianMask = EMC_maskIndex('c2c', originalSize, 'GPU',{});
fullXform(1:size(halfXform,1), :, :) = halfXform;
fullXform = fullXform(hermitianMask);
      

% % % % NOTE this is now used for things that are abs()^2 but if used more
% % % % generally an option to multiply the hermitian conjugate side by -1 will
% % % % be needed.
% % % if (flg2d)
% % %   % Shift the mirrored vectors if even to keep the origin in the correct
% % %   % place.
% % %   sY = (1-mod(nY,2));
% % %   fullXform(oX-1+isOdd:end,:,:) = halfXform;
% % %   fullXform(1:oX-2+isOdd,1+sY:end) = halfXform(oX-1+isOdd:-1:2,nY:-1:1+sY);
% % % else
% % %   sY = (1-mod(nY,2));
% % %   sZ = (1-mod(nZ,2));
% % %   fullXform(oX-1+isOdd:end,:,:) = halfXform;
% % %   fullXform(1:oX-2+isOdd,1+sY:end,1+sZ:end) = halfXform(oX-1+isOdd:-1:2,nY:-1:1+sY,nZ:-1:1+sZ);
% % % end

clear halfXform

if (abs(scalingFactor - 1.0) > 1e-3)
  fullXform =  BH_reScale3d(fullXform,'',scalingFactor,'GPU');
end

end
