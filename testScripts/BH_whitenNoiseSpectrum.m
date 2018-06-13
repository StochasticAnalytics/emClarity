function [ reWeight, radialAvg, radialVect, rFit ] = BH_whitenNoiseSpectrum( imgIN, maskIN, pixelSize, flgNorm )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Needs to know where the noise is - either provide a binary mask or create one
% with localWiener2d.

if ismatrix(imgIN)
  [d1,d2] = size(imgIN);
else
  [d1,d2] = size(sum(imgIN,3));
end


imgIN = gpuArray(imgIN);

if isnumeric(maskIN)
  if size(maskIN) ~= [d1,d2]
    error('mask and img size do not match')
  end
  maskIN = gpuArray(maskIN);
elseif numel(pixelSize) > 1
 bandpass = pixelSize;
 pixelSize = pixelSize(3);
 [~,~,~,maskIN] = BH_localWiener2d( imgIN,bandpass,1,1);
else
  maskIN = 0;
end

% Invert the mask and get power spectra
% % % PS = abs(fftn(imgIN.*(1-maskIN)));
if isreal(imgIN(1))
  PS = abs(fftn(imgIN.*(1-maskIN))).^2;
else
  PS = abs(imgIN).^2;
end


[radialGrid] = BH_multi_gridCoordinates([d1,d2],'Cartesian',...
                                                         'GPU',{'none'},1,0,1);

radialGrid = radialGrid ./ pixelSize;                                                       
[gKernel] = BH_multi_gaussian2d(-1.*[d1,d2],3,0);                                                       
nMax = 4;                                                
% Divide bins into roughly equal areas, with min ~ 2 at nyquist
oX = ceil((d1+1)/2);
radialMax = oX;
radialVect = radialGrid(1:radialMax,1);
% radialAvg = zeros([d1,d2],'single','gpuArray');


minPixPerBin = ceil(nMax*2*pi*radialMax);

i = 2;
n = floor(0.005*oX);
idxVect = n;
while n < oX
  idxVect = [idxVect (n+floor((idxVect(i-1).^.25)))];
  n = idxVect(i);
  i = i + 1;
end

whiteningCutoff = 1; % Ang

idxVect = idxVect(1:end-1);
idxVect(end) = oX;
radialAvg = zeros([length(idxVect)-1,1],'single','gpuArray');
radialSampledAt = zeros([length(idxVect)-1,1],'single','gpuArray');

nRing = 1;
nShrink = 0;
for iRing = 1:length(idxVect)-1
  if (radialVect(idxVect(iRing+1)) < 1/whiteningCutoff)
    avgMask = (radialGrid < radialVect(idxVect(iRing+1)) & radialGrid >= radialVect(idxVect(iRing)));
    currentSum = sum(PS(avgMask))/sum(avgMask(:));
    radialAvg(nRing) = currentSum;
  else
    radialAvg(nRing) = currentSum*(1+.0001*nShrink);
    nShrink = nShrink + 1;
  end
  radialSampledAt(nRing) = radialVect(idxVect(iRing+1));
    nRing = nRing +1;
end



clear avgMask PS


 rFit = fit(gather(double(radialSampledAt)),gather(double(radialAvg).^0.5),'spline');
%   figure, plot(radialVect,rFit(radialVect))
%   figure, plot(radialVect,1./rFit(radialVect));

% radialAvg = radialAvg.^-.5 .* radialMask;
% figure, imshow3D(gather(radialAvg))
radialAvg = reshape(rFit(radialGrid),d1,d2);
reWeight = real(ifftn(fftn(imgIN)./radialAvg));
% reWeight = reWeight ./ std(reWeight(maskIN < 0.01))^1;
if (flgNorm(1))
  reWeight = reWeight - mean(reWeight(:));
  reWeight = reWeight ./ rms(reWeight(:));
end


end

