function [filtIMG,lVar,lMean,particleMask] = BH_localWiener2d( IMAGE, bandPass,...
  flgZeroMean,...
  flgParticle,varargin)
%Calculate a locally adaptive wiener filter and smooth image
%
% Use a gaussian disc with 2*stdDev = radius of area for calculating stats.
% Estimate local noise variance as the average of all local variances

% For now it seems to be pretty consistent that a radius of ~ 6-8 A works
% well.

% Need to add a check for Nans and/or just pre-normalize things or better
% yet also include a whitening step.

pixelSize = bandPass(3);
% RADIUS is used for the calculation of local statistics, dividing by 3
% sets the gaussian for this process to fall to e^-3 by the particle
% radius.
if numel(bandPass) < 4
  RADIUS = (110./3)/pixelSize;
else
  RADIUS = (bandPass(4)/3)/pixelSize;
end
flgIterate = 1;
reduceSmoothing = 1;

if strcmp(class(bandPass),'gpuArray')
  useGPU = 1;
else
  useGPU = 0;
end

if ismatrix(IMAGE)
  [d1,d2] = size(IMAGE);
  flgMovie = 0;
else
  [d1,d2,d3] = size(IMAGE);
  flgMovie = 1;
end

if nargin == 5
  totalDose = varargin{1};
else
  totalDose = 38;
end


if ( useGPU )
  bandPassFilt = BH_bandpass3d([d1,d2,1],10^-6,bandPass(1),bandPass(2), ...
    'GPU',pixelSize);
else
  bandPassFilt = BH_bandpass3d([d1,d2,1],10^-6,bandPass(1),bandPass(2), ...
    'cpu',pixelSize);
end

if ( flgMovie )
  if ( useGPU )
    imgSUM = real(ifftn(fftn(gpuArray(sum(IMAGE,3))).*bandPassFilt));
  else
    imgSUM = real(ifftn(fftn(sum(IMAGE,3)).*bandPassFilt));
  end
else
  imgSUM = real(ifftn(fftn(IMAGE).*bandPassFilt));
end


% Although outliers should have already been cleaned up, there is a chance
% for outliers. For now, try checking the kurtosis, and given a deviation

K = kurtosis(imgSUM(:));
if abs(K - 3) > 0.75
  
  
  medKernel = min(9,max(3,ceil(sqrt(abs(K-3)))));
  % medfilt wants an odd kernel size
  medKernel = (medKernel + (1-mod(medKernel,2)));
  
  % For some stupid effing reason, medfilt2 returns an error if I pass the
  % medKernel directly, but will take a number of the same class
  %   fprintf(['Found strong outliers, kurtosis = %3.3f,',...
  %            'running median filter size %d\n'],K,medKernel(1));
  
  switch medKernel
    case 3
      imgSUM = medfilt2(imgSUM,[3,3]);
    case 5
      imgSUM = medfilt2(imgSUM,[5,5]);
    case 7
      imgSUM = medfilt2(imgSUM,[7,7]);
    case 9
      imgSUM = medfilt2(imgSUM,[9,9]);
  end
  
  
end


% from normal (k=3) run a median filter
% the 4th root localizes the kernal more.
% The kernel is already normalized, st the area underneath is already sum=1
gFT = BH_multi_gaussian2d(-1.*[d1,d2], RADIUS,0);
[g1x,g1y] = BH_multi_gaussian2d(-1.*[d1,d2], RADIUS,1);

if ( useGPU )
  gFT = gpuArray(gFT);
end

% Calculate the local statistics
lMean = real(ifftn(fftn(imgSUM).*gFT));
imgSUM = imgSUM - lMean;
lVar  = real(ifftn(fftn((imgSUM).^2).*gFT)).^1;

% lVardY = gather(real(ifftn(fftn((imgSUM).^2).*g1y).^1));
% lVardX = gather(real(ifftn(fftn((imgSUM).^2).*g1x).^1));



% Estimate the noise variance simply as the average from the local
% variances as matlab does. This could likely be improved by considering a
% histogram and some cutoff -- this may be necessary with carbon in the
% image.
% lNoiseVar = mean2(lVar);
%lNoiseVar = median(lVar(:))
inc = (max(lVar(:))-min(lVar(:)))/1000;
x = min(lVar(:)):inc:max(lVar(:))-inc;
h = hist(lVar(:),x);
% figure, plot(x,h)
% It would be nice to have the option to exclude high variance outliers in
% a reasonable way
% if nargin == 5
%   figure, plot(x,h)
%   figure, plot(x,convn(h,fspecial('gaussian',[1,128],32),'same'))
%   figure, plot(x(1:end-2),diff(diff(convn(h,fspecial('gaussian',[1,128],32),'same'))))
%   figure, imshow3D(gather(lVar>.02))
%   return
% end
[maxR,maxC]= max(h(:));

% Where the histogram drops to 96% its max AFTER reaching it.
try
  lNoiseVar = x(-1.*find(h(:)>(.96)*maxR,1,'first')+2*(maxC));
catch
  save('wienerFail.mat','inc','x','h','maxR','maxC','lVar','lMean','imgSUM','IMAGE','gFT');
  error('wienerFiltFail on hist max')
end
% lNotExtremeVar = reshape(lVar(:) < 3.5*std(lVar(:))+mean(lVar(:)),size(lVar));


if ( flgZeroMean )
  wienerFilt = (max(lVar-lNoiseVar,0)./(lVar)).^2;
else
  wienerFilt = lMean + (max(lVar-lNoiseVar,0))./lVar;
end


if ( flgParticle || flgIterate )
  [ particleMask ] = calc_particleMask(wienerFilt.*imgSUM, useGPU,totalDose);
else
  particleMask = '';
end

if ( flgIterate )
  maxIter = 10;
  lastRatio = lNoiseVar;
  iIter = 1;
  while lastRatio > 1.001 && iIter <= maxIter
    if iIter == 1
      oldBestEst = lNoiseVar;
    else
      oldBestEst = newBestEst;
    end
    newBestEst = mean(lVar(particleMask<0.05));
    fprintf('Current noise variance estimates are %3.3e and now %3.3e\n',...
      lNoiseVar, newBestEst);
    
    if ( flgZeroMean )
      wienerFilt = (max(lVar-newBestEst,0)./(lVar)).^2;
    else
      wienerFilt = lMean + (max(lVar-newBestEst,0))./lVar;
    end
    
    [ particleMask ] = calc_particleMask(wienerFilt.*imgSUM, useGPU,totalDose);
    iIter = iIter + 1;
    lastRatio = newBestEst/oldBestEst;
  end
  
  
end

if ( flgMovie )
  
  filtIMG = zeros([d1,d2,d3],'single');
  
  for iPrj = 1:d3
    if iPrj == 1
      iIMG = real(ifftn(fftn(sum(IMAGE(:,:,2:end),3)).*bandPassFilt));
    elseif iPrj < d3
      iIMG = real(ifftn(fftn(sum(IMAGE(:,:,[1:iPrj-1,iPrj+1:end]),3)).*bandPassFilt));
    else
      iIMG = real(ifftn(fftn(sum(IMAGE(:,:,1:end-1),3)).*bandPassFilt));
    end
    % Scale the mean
    filtIMG(:,:,iPrj) = gather(wienerFilt.*(iIMG-(lMean.*((d3-1)/d3))));
    filtIMG(:,:,iPrj) = filtIMG(:,:,iPrj) ./ rms(filtIMG(:,:,iPrj));
  end
else
  filtIMG = wienerFilt.*IMAGE;
  filtIMG = filtIMG - mean(filtIMG(:));
  filtIMG = filtIMG ./ rms(filtIMG(:));
end

clear IMAGE imgSUM wienerFilt iIMG
end

function [ particleMask ] = calc_particleMask(filtIMG, useGPU,totalDose)

% Use the filtered image to make a simplified mask
%%% ADD adaptive cutoff based on number of initial pixels wanted (1% for
%%% example)
filtIMG = filtIMG ./ rms(filtIMG(:));
%     intCutoff = 3.0;

intCutoff = 0.5*sqrt(totalDose);
[gx] = BH_multi_gaussian2d(2.*[5,5],2.0,0);


if ( useGPU )
  gx = gpuArray(gx);
end

% avoid disconnected densitys by setting high initial cutoff
dilatedIMG = convn(filtIMG,gx,'same') > 0.01; %real(ifftn(fftn(filtIMG).*gx));
binaryCut = (abs(filtIMG) .* dilatedIMG) > intCutoff;
sum(binaryCut(:))./numel(binaryCut);


% relax the cutoff and iteratively dilate the region
intCutoff = intCutoff./3 ;
growthFactor = 1.2;
while growthFactor > 1.05
  nStart = sum(binaryCut(:));
  
  dilatedIMG = convn(binaryCut,gx,'same') > 0.01; %real(ifftn(fftn(binaryCut).*gx));
  binaryCut  = ( abs(filtIMG) .* dilatedIMG ) > intCutoff;
  growthFactor = sum(binaryCut(:))/nStart;
end


% relax the cutoff and iteratively dilate the region
intCutoff = intCutoff./5;

growthFactor = 1.2;
while growthFactor > 1.05
  nStart = sum(binaryCut(:));
  
  dilatedIMG = convn(binaryCut,gx,'same') > 0.01; %real(ifftn(fftn(binaryCut).*gx));
  binaryCut  = ( abs(filtIMG) .* dilatedIMG) > intCutoff;
  growthFactor = sum(binaryCut(:))/nStart;
end
%   particleMask = real(ifftn(fftn(binaryCut).*gx.^2));
[gx] = BH_multi_gaussian2d(7.*[1,1],3.0,0);
particleMask = convn(convn(binaryCut,gx,'same'),gx,'same');
particleMask = particleMask ./ max(particleMask(:));
end
