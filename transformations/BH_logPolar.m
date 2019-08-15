

% stackIn = getVolume(MRCImage('tilt01.fixed'),-1,-1,[7:8]);

pixelSize = 1.0;
imgRot = gpuArray(stackIn(:,:,2));
imgIn = gpuArray(stackIn(:,:,1));
tiltAxis = -85;
                                
% abs, sqrt, log - spectrum to use for fitting
spectrumToUse = 'log'

minRes = 800;
maxRes = 12
computeSize = 1024;

angMax = pi;
angMin = 0;


tiltAxisSearch = 81.5:1:105;
tiltAxisCCC = zeros(length(tiltAxisSearch),1,'gpuArray');
tiltAngles = [1:7,23:30];
nTilts = length(tiltAngles);
mScores = 0.*tiltAngles;

%%%%%%%%%%% search for the tilt axis
[d1,d2,d3] = size(stackIn); 
oX = floor(d1/2) + 1;
oY = floor(d2/2) + 1;
% First crop the image to kill high frequencies, then expand to
bp = BH_bandpass3d([d1,d2],0,minRes,maxRes,'GPU',pixelSize);

% % % mask = (BH_mask3d('sphere',[d1,d2],min([oX,oY] - 7).*[1,1],[0,0],'2d'));
% % % for iAng = 1:nTilts
% % %   imgRot = gpuArray(stackIn(:,:,tiltAngles(iAng)));
% % % for i = 1:length(tiltAxisSearch)
% % % 	tmp = (BH_resample2d(imgRot.*mask,[tiltAxisSearch(i),0,0],[0,0,0],'Bah','GPU','inv',1,[d1,d2]));
% % %   tmp = fftn(tmp(oX-192:oX+191,:));
% % %   tmp = tmp(bp);
% % %   tiltAxisCCC(i) = sum(sqrt(abs(tmp(:))));
% % % end
% % % [m,c] = max(tiltAxisCCC)
% % % tiltAxis = tiltAxisSearch(c)
% % %   mScores(iAng) = tiltAxis
% % % end
% % % % figure, plot(tiltAxisSearch,tiltAxisCCC)
% % % 
% % % 
% % % return

targetScaleFactor = 2*pixelSize/maxRes;
if (targetScaleFactor > 1)
  error('range error in maxRes')
end

% TODO make this a nice FT size since the res limit doesn't need to be
% exact
trimSize = floor(targetScaleFactor.*size(imgIn))
scaleFactor = trimSize ./size(imgIn)

tX = floor(trimSize(1)/2) + 1;
tY = floor(trimSize(2)/2) + 1;

trimVal = BH_multi_padVal(size(imgIn),trimSize);
padVal  = BH_multi_padVal(trimSize,computeSize.*[1,1]);

imgIn = gpuArray(imgIn);
imgIn = fftshift(fftn(imgIn).*bp);
imgIn = BH_padZeros3d(imgIn,'fwd',trimVal,'GPU','single');
imgIn = real(ifftn(ifftshift(imgIn)));
imgIn = BH_padZeros3d(imgIn,'fwd',padVal,'GPU','single');

imgRot = gpuArray(imgRot);
imgRot = fftshift(fftn(imgRot).*bp);
imgRot = BH_padZeros3d(imgRot,'fwd',trimVal,'GPU','single');
imgRot = real(ifftn(ifftshift(imgRot)));
imgRot = BH_padZeros3d(imgRot,'fwd',padVal,'GPU','single');

SAVE_IMG(MRCImage(gather(imgIn)),'~/tmp/imgIn_pre.mrc');
SAVE_IMG(MRCImage(gather(imgRot)),'~/tmp/imgRot_pre.mrc');

[d1,d2] = size(imgIn); 
oX = floor(d1/2) + 1;
oY = floor(d2/2) + 1;

% [gridX, gridY ] = BH_multi_gridCoordinates([d1,d2],'Cartesian','GPU', ...
%                   {'single',BH_defineMatrix([tiltAxis,0,0],'Bah','forward'),[0,0,0]','forward',1,1},0,1,0);
% figure, imshow3D(gather(gridX))
% return

mask_FT = (BH_mask3d('sphere',[d1,d2],min([oX,oY] - 7).*[1,1],[0,0],'2d'));
mask = (BH_mask3d('sphere',[d1,d2],min([tX,tY] - 7).*[1,1],[0,0],'2d'));


imgIn = imgIn .* mask;
imgRot = imgRot.*mask;
% imgRot = mask.*BH_resample2d(imgIn,[-23,0,0],[0,0,0],'Bah','GPU','forward',0.9,size(imgIn));

imgIn_FT = fftshift(fftn(imgIn));
imgRot_FT = fftshift(fftn(imgRot));

switch spectrumToUse
  case 'abs'
    imgIn_FT = abs(imgIn_FT);
    imgRot_FT = abs(imgRot_FT);
  case 'sqrt'
    imgIn_FT = sqrt(abs(imgIn_FT));
    imgRot_FT = sqrt(abs(imgRot_FT));
  case 'log'
    imgIn_FT = log(abs(imgIn_FT));
    imgRot_FT = log(abs(imgRot_FT));
  otherwise
    error('abs sqrt log on spectrum');
end

radMax = log(oX);
radMin = 0;
radInc = log(oX)./oX;
angInc = (angMax - angMin)./d2;



[rad,ang] = ndgrid(gpuArray(radMin:radInc:radMax), gpuArray(angMin:angInc:angMax));

[x,y] = BH_multi_gridCoordinates([d1,d2],'Cartesian','GPU',{'none'},0,-1,0);

xInterp = exp(rad).* (tan(ang).^2  + 1).^-(1/2);
yInterp = exp(rad).* (tan(ang).^-2 + 1).^-(1/2);

% Because we have tan^2 the quadrants need to be addressed
angMask = ang > pi/2 & ang <= pi;
xInterp(angMask) = -1.*xInterp(angMask);

angMask = ang > pi & ang <= 1.5*pi;
xInterp(angMask) = -1.*xInterp(angMask);
yInterp(angMask) = -1.*yInterp(angMask);

angMask = ang > 1.5*pi & ang <= 2*pi;
yInterp(angMask) = -1.*yInterp(angMask);

xInterp = xInterp + oX;
yInterp = yInterp + oY;

logPolarImg = interpn(x,y,imgIn_FT.*mask_FT,xInterp,yInterp,'linear',0);

figure, imshow3D(gather(imgIn_FT));
figure, imshow3D(gather(logPolarImg));

% Just for visualization, rescale the mask as well
% mask = BH_resample2d(mask,[20,0,0],[0,0,0],'Bah','GPU','forward',0.9,size(imgIn_FT));

figure, imshow3D(gather(imgRot_FT));

logPolarRot =  interpn(x,y,imgRot_FT.*mask_FT,xInterp,yInterp,'linear',0);
figure, imshow3D(gather(logPolarRot))

logPolarImg = logPolarImg - mean(logPolarImg(:));
logPolarImg = logPolarImg ./ rms(logPolarImg(:));

logPolarRot = logPolarRot - mean(logPolarRot(:));
logPolarRot = logPolarRot ./ rms(logPolarRot(:));

iCCC = fftshift(real(ifftn(fftn(logPolarRot).*conj(fftn(logPolarImg)))));
iCCC = iCCC ./ std(iCCC(:));

[m,c] = max(iCCC(:));
[i,j] = ind2sub(size(iCCC),c);

[bx,by] = ndgrid(-1:1,-1:1);
comBox = iCCC(i-1:i+1,j-1:j+1);
bx = bx(:);
by = by(:);
comBox = comBox(:);

com = [sum(bx.*comBox), sum(by.*comBox)] ./ sum(comBox(:));
figure, imshow3D(iCCC); hold on
plot(j,i,'bo');


sx = i - (floor(size(iCCC,1)/2) + 1.0 )
sy = j - (floor(size(iCCC,2)/2) + 1.0 )

% Needs to be inverted because we are looking at a power spectru
magShift = 1./exp(((com(1)+sx).*radInc))
angShift = (com(2) + sy).*angInc*180/pi

% Calc the translational shift without considering rotation
iCCC = fftshift(real(ifftn(fftn(imgRot).*conj(fftn(imgIn)))));
iCCC = iCCC ./ std(iCCC(:));

[m,c] = max(iCCC(:));
[i,j] = ind2sub(size(iCCC),c);

[bx,by] = ndgrid(-1:1,-1:1);
comBox = iCCC(i-1:i+1,j-1:j+1);
bx = bx(:);
by = by(:);
comBox = comBox(:);

com = [sum(bx.*comBox), sum(by.*comBox)] ./ sum(comBox(:));
figure, imshow3D(iCCC); hold on
plot(j,i,'bo');

sx = (i - (floor(size(iCCC,1)/2) +1)) ./ scaleFactor(1)
sy = (j - (floor(size(iCCC,2)/2) +1)) ./ scaleFactor(2)

% Repeat but with the rotated image
% Calc the translational shift without considering rotation
imgRot = BH_resample2d(imgRot,[angShift,0,0],[0,0,0],'Bah','GPU','forward',magShift,size(imgIn_FT));


% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % 
% % % % Only 6 FFTs are needed.
% % % refFT = conj(fftn(imgIn));
% % % imgFT = fftn(imgRot);
% % % mskFT = fftn(mask);
% % % ref_FT_SQ = conj(fftn(imgIn.^2));
% % % img_FT_SQ = fftn(imgRot.^2);
% % % 
% % % iF1F2 = real(ifftn(imgFT.*refFT));
% % % iF1M = real(ifftn(imgFT.*conj(mskFT)));
% % % iF2M = real(ifftn(refFT.*mskFT));
% % % iMM = real(ifftn(mskFT.*conj(mskFT)));
% % % iF1SQ = real(ifftn(ref_FT_SQ.*conj(mskFT)));
% % % iF2SQ = real(ifftn(img_FT_SQ.*mskFT));
% % % 
% % % 
% % % denom_1 = (iF1SQ - iF1M.^2./iMM).*(iF2SQ - iF2M.^2./iMM);
% % % denom_1(~isfinite(denom_1)) = 1;
% % % 
% % % denom_1(denom_1 < 1e-6) = 1;
% % % 
% % % denom_1 = sqrt(denom_1);
% % % 
% % % iCCC = (iF1F2 - iF1M.*iF2M./iMM ) ./ denom_1;
% % % min(iCCC(:))
% % % max(iCCC(:))
% % % std(iCCC(:))
% % % iCCC = iCCC ./ std(iCCC(:));
% % % figure, imshow3D(gather(iCCC))
iCCC = real(ifftn(fftn(imgRot).*conj(fftn(imgIn))));
iCCC = fftshift(iCCC);
figure, imshow3D(gather(iCCC))

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% iCCC = fftshift(real(ifftn(fftn(imgRot).*conj(fftn(imgIn)))));
% iCCC = iCCC ./ std(iCCC(:));
% 
% maskBinary = mask > 0.01;
% mask = fftshift(real(ifftn(fftn(mask).*conj(fftn(mask)))));
% mask = mask ./ max(mask(:));
% mask = 1./mask.* maskBinary;
% mask(maskBinary) = sqrt(mask(maskBinary));
% figure, imshow3D(gather(mask));
% iCCC = iCCC .* mask;

[m,c] = max(iCCC(:));
[i,j] = ind2sub(size(iCCC),c);

[bx,by] = ndgrid(-1:1,-1:1);
comBox = iCCC(i-1:i+1,j-1:j+1);
bx = bx(:);
by = by(:);
comBox = comBox(:);

com = [sum(bx.*comBox), sum(by.*comBox)] ./ sum(comBox(:));
figure, imshow3D(iCCC); hold on
plot(j,i,'bo');

sx = (i - (floor(size(iCCC,1)/2) +1)) ./ scaleFactor(1)
sy = (j - (floor(size(iCCC,2)/2) +1)) ./ scaleFactor(2)

imgRot = BH_resample2d(imgRot,[0,0,0],scaleFactor.*[sx,sy,],'Bah','GPU','inv',1,size(imgIn_FT));

SAVE_IMG(MRCImage(gather(imgIn)),'~/tmp/imgIn.mrc');
SAVE_IMG(MRCImage(gather(imgRot)),'~/tmp/imgRot.mrc');
