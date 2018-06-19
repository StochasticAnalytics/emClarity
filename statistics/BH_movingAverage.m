function [ AVG, flgOOM ] = BH_movingAverage( img, featureSize )
%Calculate a moving average using circular convolution.
%   Detailed explanation goes here
flgOOM = 0;
if ndims(img) ~= numel(featureSize)
  error('image must have same dimension as number of dims given in feature size.')
end

if isa(img,'gpuArray')
  useGPU = 1;
else
  useGPU = 0;
  featureSize = gather(featureSize);
end

[d1,d2,d3] = size(img);

% Pad the image by 1/2 the feature size, and the mask to the corresponding
% total size.

padImg = [floor(featureSize./2);ceil(featureSize./2)];
padAvg = [floor(size(img)./2);ceil(size(img)./2)];
totalSize = gather(size(img)+sum(padImg,1));

% use the same type as the image
if isa(gather(img(1)), 'single')
  precision = 'single';
elseif isa(gather(img(1)), 'double')
  precision = 'double';
else
  error('Img must be float or double!\n')
end

try
  if numel(featureSize) == 2
    meanFilterCore = BH_mask3d('sphere',featureSize,featureSize./2,[0,0],'2d');
  else
    meanFilterCore = BH_mask3d('sphere',featureSize,featureSize./2,[0,0,0]);
  end
  
  if (useGPU)
    paddedImg = zeros(totalSize, precision, 'gpuArray');     
    meanFilter = zeros(totalSize, precision, 'gpuArray');
      
  else
    meanFilter = zeros(totalSize, precision);
    paddedImg = zeros(totalSize, precision);
    meanFilterCore = gather(meanFilterCore);
  end
  
  
catch
  error('error in movingAvg, %f %f totalSize, %s precision\n',totalSize,precision);
end

if numel(featureSize) == 2
  
 meanFilter(padAvg(1,1)+1:end-padAvg(2,1), ...
            padAvg(1,2)+1:end-padAvg(2,2)) = meanFilterCore;
  
%  meanFilter(padAvg(1,1)+1:end-padAvg(2,1), ...
%             padAvg(1,2)+1:end-padAvg(2,2)) = ...
%                                zeros(featureSize) + prod(featureSize)^-1;

%   meanFilter(padAvg(1,1)+1:end-padAvg(2,1), ...
%              padAvg(1,2)+1:end-padAvg(2,2)) = ...
%   meanFilter(padAvg(1,1)+1:end-padAvg(2,1), ...
%              padAvg(1,2)+1:end-padAvg(2,2)) + prod(featureSize)^-1;
                                 
  meanFilter = meanFilter ./ sum(meanFilter(:));
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
  clear img    
  
  AVG = fftshift(real(ifftn(fftn(meanFilter).*fftn(paddedImg))));
  AVG = AVG(padImg(1,1)+1:end-padImg(2,1), ...
            padImg(1,2)+1:end-padImg(2,2));
          
else
  
  meanFilter(padAvg(1,1)+1:end-padAvg(2,1), ...
             padAvg(1,2)+1:end-padAvg(2,2), ...
             padAvg(1,3)+1:end-padAvg(2,3)) = meanFilterCore;
%   meanFilter(padAvg(1,1)+1:end-padAvg(2,1), ...
%              padAvg(1,2)+1:end-padAvg(2,2), ...
%              padAvg(1,3)+1:end-padAvg(2,3)) = ...
%                                 zeros(featureSize) + prod(featureSize)^-1;
  meanFilter = meanFilter ./ sum(meanFilter(:));

  paddedImg(padImg(1,1)+1:end-padImg(2,1), ...
            padImg(1,2)+1:end-padImg(2,2), ...
            padImg(1,3)+1:end-padImg(2,3)) = img;
  clear img
  
  % The following results in oom in 16 but not 17b?
  try
    AVG = fftshift(real(ifftn(fftn(meanFilter).*fftn(paddedImg))));
  catch

    AVG = fftshift(real(ifftn(fftn(gather(meanFilter)).*fftn(gather(paddedImg)))));
    clear meanFilter paddedImg
    fprintf(['\n\nin calculating local statistics, ran out of memory.\n',...
             ' problem seems to only happen in versions older than 2017b\n',...
             ' so performing this on the cpu at a small penality in speed\n\n.'])
    flgOOM = 1;

  end
  
  AVG = AVG(padImg(1,1)+1:end-padImg(2,1), ...
            padImg(1,2)+1:end-padImg(2,2), ...
            padImg(1,3)+1:end-padImg(2,3));
end                                
clear ans
clearvars -except AVG flgOOM
end

