function [ SF3D ] = BH_weightMaskMex(SIZE, SAMPLING, TLT, subtomo_origin_in_tomo_frame, reconGeometry, wiener_constant)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


doHalfMask = false;
doSqCTF = true;

SIZE = gather(uint32(SIZE));

TLT = gather(single(sortrows(TLT,1)));
tiltAngles = TLT(:,4);

% We are only using these coordinates to figure out a change in defocus, so there is no need to worry about sampling rate.
pixelSize_angstrom = gather(single(TLT(:,16).*10^10));
tomo_origin_wrt_tilt_origin = [reconGeometry.dX_specimen_to_tomo, ...
                                reconGeometry.dY_specimen_to_tomo, ...
                                reconGeometry.dZ_specimen_to_tomo];              
tomo_origin_in_tomo_frame = emc_get_origin_index([reconGeometry.NX, ...
                                                  reconGeometry.NY, ...
                                                  reconGeometry.NZ]);

subtomo_origin_in_specimen_frame = subtomo_origin_in_tomo_frame - tomo_origin_in_tomo_frame + tomo_origin_wrt_tilt_origin;

iCs = single(TLT(:,17).*10^3);
iWavelength = single(TLT(:,18).*10^10);
iPhaseShift = TLT(:,19)  ;
iDefocus = single(abs(TLT(:,15)).*10^10)  ;
iddF = single(TLT(:,12).*10^10);
idPHI = single(TLT(:,13).*180.0/pi);

nTilts =   uint32(size(TLT,1));
exposure = gather(single(TLT(:,11)));

% Need a defocus offset based on XYZ position in the tomogram
for iPrj = 1:nTilts
  % If there are local alignments, then this isn't quite right, but they additional rotations about Z are < 1 degree, which
  % shouldn't have a major impact here.
  rTilt =  BH_defineMatrix(TLT(iPrj,4),'TILT','fwdVector') ;
  prjCoords = rTilt * subtomo_origin_in_specimen_frame';
  iDefocus(iPrj) = iDefocus(iPrj)-(prjCoords(3).*pixelSize_angstrom(iPrj));
end

iDefocus = gather(single(iDefocus));
% TODO figure out how to get this from the data
iThickness = 75; % nm
fractionOfDose = gather(single(TLT(:,14)./mean(TLT(:,14))));
fractionOfElastics = exp(-1.*iThickness./( cosd(TLT(:,4)).*400 ));
fractionOfElastics = fractionOfElastics ./ max(fractionOfElastics(:));


[SF3D] = mexSF3D(doHalfMask,doSqCTF,SIZE,pixelSize_angstrom * SAMPLING,iWavelength,iCs, ...
  gather(single(iDefocus + iddF)), ...
  gather(single(iDefocus - iddF)), ...
  idPHI,iPhaseShift,nTilts,tiltAngles, ...
  exposure,fractionOfElastics.*fractionOfDose,int16(1), ...
  gather(single(wiener_constant)));


% SF3D = SF3D ./ (WGT+0.01);
SF3D = SF3D - (min(SF3D(:)) - 1e-4);
% figure, hist(SF3D(:))
% SF3D = SF3D ./ max(SF3D(:));

end

