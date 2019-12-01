function [ reconGeom ] = BH_offsets( rCoords, tiltName, SAMPLING)
%UNTITLED Summary of this function goes here
%      rCoords = subTomoMeta.mapBackGeometry.(tiltName).coords(tomoNumber,:);

rCoords = rCoords ./ SAMPLING;
rCoords(1:4) = fix(rCoords(1:4));

if exist(tiltName,'file')
  header = getHeader(MRCImage(tiltName,0));  
else
  error('tiltName %s does not exist\n', tiltName);
end

    oY = ceil((header.nY +1)./2);
    nY = rCoords(3) - rCoords(2) + 1;
    dY = floor(rCoords(2)+nY/2) -oY;  
    reconGeom = zeros(2,3);
    reconGeom(1,1:3) = [rCoords(1), nY, rCoords(4)];
    % value specify location of origin, but SHIFT in IMOD's tilt takes the
    % location to shift the origin too, so multiply oX by -1. The notion for Z is
    % flipped since imod does reconstruction on a rotated frame. I.e. a positive
    % number shifts the recon "up" which when rotated to the microscope frame is
    % actually "down" (in Z) so no need to multipliy oZ by -1.
    reconGeom(2,1:3)   = round([-1*rCoords(5),dY,rCoords(6)]);

end



