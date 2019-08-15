clear

imName = [matlabroot '/toolbox/images/imdata/sherlock.jpg'];

imgIn = imread(imName);
imgIn = single(imgIn(:,:,1));

angMax = 2.*pi;
angMin = 0;


[d1,d2] = size(imgIn); 
oX = floor(d1/2) + 1;
oY = floor(d2/2) + 1;

mask = gather(BH_mask3d('sphere',[d1,d2],min([oX,oY] - 7).*[1,1],[0,0],'2d'));
imgIn = imgIn .* mask;

radMax = log(oX);
radMin = 0;
radInc = log(oX)./oX;
angInc = (angMax - angMin)./d2;



[rad,ang] = ndgrid(radMin:radInc:radMax, angMin:angInc:angMax);
bp = BH_bandpass3d(size(rad),0,800,4,'cpu',1);

[x,y] = BH_multi_gridCoordinates([d1,d2],'Cartesian','cpu',{'none'},0,-1,0);

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


logPolarImg = interpn(x,y,imgIn,xInterp,yInterp,'linear',0);

figure, imshow3D(imgIn.*mask)
figure, imshow3D(logPolarImg);

imgRot = BH_resample2d(imgIn,[-23,0,0],[0,0,0],'Bah','cpu','forward',1.2,size(imgIn));
% Just for visualization, rescale the mask as well
mask = BH_resample2d(mask,[20,0,0],[0,0,0],'Bah','cpu','forward',0.9,size(imgIn));

figure, imshow3D(imgRot.*mask);

logPolarRot =  interpn(x,y,imgRot,xInterp,yInterp,'linear',0);
figure, imshow3D(logPolarRot)

logPolarImg = logPolarImg - mean(logPolarImg(:));
logPolarImg = logPolarImg ./ rms(logPolarImg(:));

logPolarRot = logPolarRot - mean(logPolarRot(:));
logPolarRot = logPolarRot ./ rms(logPolarRot(:));

iCCC = fftshift(real(ifftn(fftn(logPolarRot).*conj(fftn(logPolarImg).*bp))));
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

magShift = exp(((com(1)+sx).*radInc))
angShift = (com(2) + sy).*angInc*180/pi


