function [ SF3D ] = BH_weightMaskMex(SIZE, SAMPLING, TLT, ...
                                     xyzSubTomo,reconGeometry)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


doHalfMask = false;
doSqCTF = true;

SIZE = gather(uint32(SIZE));

TLT = gather(single(sortrows(TLT,1)));
tiltAngles = TLT(:,4);
% reconGeometry = subTomoMeta.('reconGeometry').(tomoList{iTomo});

pixelSize = gather(single(TLT(:,16).*SAMPLING.*10^10));
reconShift = reconGeometry(2,:); % already at the appropriate sampling rate.
originVol = ceil((reconGeometry(1,1:3)+1)./2);

prjVector = xyzSubTomo - originVol + reconShift;
        
iCs = single(TLT(:,17).*10^3);
iWavelength = single(TLT(:,18).*10^10);
iPhaseShift = TLT(:,19)  ;
iDefocus = single(-TLT(:,15).*10^10)  ;
iddF = single(TLT(:,12).*10^10);
idPHI = single(TLT(:,13).*180.0/pi);

nTilts =   uint32(size(TLT,1));
exposure = gather(single(TLT(:,11)));
        
% Need a defocus offset based on XYZ position in the tomogram
for iPrj = 1:nTilts
  rTilt =  BH_defineMatrix(TLT(iPrj,4),'TILT','forwardVector') ;
  prjCoords = rTilt * prjVector';
  iDefocus(iPrj) = iDefocus(iPrj)-(prjCoords(3).*pixelSize(iPrj));
end
  
iDefocus = gather(single(iDefocus));
% TODO figure out how to get this from the data
iThickness = 75; % nm
fractionOfDose = gather(single(TLT(:,14)./mean(TLT(:,14))));
fractionOfElastics = exp(-1.*iThickness./( cosd(TLT(:,4)).*400 ));
fractionOfElastics = fractionOfElastics ./ max(fractionOfElastics(:));


 [SF3D, WGT] = mexSF3D(doHalfMask,doSqCTF,SIZE,pixelSize,iWavelength,iCs, ...
                gather(single(iDefocus + iddF)), ...
                gather(single(iDefocus - iddF)), ...
                idPHI,iPhaseShift,nTilts,tiltAngles, ...
                exposure,fractionOfElastics.*fractionOfDose,int16(1));


SF3D = SF3D ./ (WGT+0.01);
SF3D = SF3D - (min(SF3D(:)) + 1e-4);
SF3D = SF3D ./ max(SF3D(:));

end

