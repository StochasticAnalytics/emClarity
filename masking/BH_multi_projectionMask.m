function [ evalMask, deltaZ ] = BH_multi_projectionMask( SIZE, TLT, METHOD, zTOLERANCE )
%Calculate the maximum valid area in a slap projection, and z-coords there.
%   Reduce computational burden in CTF correction and get defocus gradient
%   values for a number of circumstaces.
%   As with other programs, now assuming the tilt-axis = Y

% SIZE(1,1:3) are the dimensions of the stack of projections
% SIZE(2,1:3) are the X/Y dimensions of max reconstuction and estimated
% thicknes in Z of the specimen

if nargin == 4
  zShift  = zTOLERANCE;
else
  zShift = 0;
end
d1 = SIZE(1,1); d2 = SIZE(1,2); d3 = SIZE(1,3);
r1 = SIZE(2,1); r2 = SIZE(2,2); r3 = SIZE(2,3)/2;


% Do calc on GPU either way, and if cpu pull results
[X1,Y1,~,~,~,~] = BH_multi_gridCoordinates([r1,r2],'Cartesian', ...
  'GPU',{'none'},0,1,0);



clear deltaZ
if strcmp(METHOD,'GPU')
  deltaZ   = zeros([d1,d2,d3], 'int16', 'gpuArray');
  evalMask = false([d1,d2,d3], 'gpuArray');
  Z1 = zeros([r1,r2],'single','gpuArray');
elseif strcmp(METHOD,'cpu')
  deltaZ   = zeros([d1,d2,d3], 'int16');
  evalMask = false([d1,d2,d3]);
  Z1 = zeros([r1,r2],'single');
else
  error('METHOD must be GPU or cpu\n');
end

for iPrj = 1:d3
  
  R = BH_defineMatrix([TLT(iPrj,6),TLT(iPrj,4),TLT(iPrj,6)], 'Bah','fwdVector');
  
  rInv = R';
  
  % Z-starts at zero so it doesn't contribute
  X2 = X1.*R(1) + Y1.*R(4);
  Y2 = X1.*R(2) + Y1.*R(5);
  
  % To account for thickness, the X-dimension must be expanded
  Xmin = X2 + (Z1 - R(7)*r3);
  Xmax = X2 + (Z1 + R(7)*r3);
  
  
  % Now the gridvectors have values on a centered grid, shift back to
  % index coordinates, where the origin is at lower left = 1,1
  x3min = (round(Xmin + d1/2 ));
  x3max = (round(Xmax + d1/2));
  y3 = (round(Y2 + d2/2 ));
  
  % Remove any coordinates that are out of the original frame, the values here
  % are the x,y pairs that are sampled in the projection frame.
  x3max( x3max < 1  | x3max > d1) = 1;
  x3min( x3min < 1  | x3min > d1) = 1;
  y3( y3 < 1  | y3 > d2) = 1;
  
  
  % uniqI/J are coordinate indices in the unshifted frame that are sampled
  linInd = unique([sub2ind([d1,d2], x3min(:),y3(:));sub2ind([d1,d2], x3max(:),y3(:))]);
  [uniqI,uniqJ] = ind2sub([d1,d2], linInd);
  
  xP = uniqI - d1/2;
  yP = uniqJ - d2/2;
  
  zP = -1/rInv(9).*(xP.*rInv(3)+yP.*rInv(6));
  
  clear interpZ maskZ
  if strcmp(METHOD,'GPU')
    maskZ = false([d1,d2],'gpuArray');
    interpZ = zeros([d1,d2],'int16','gpuArray');
  else
    maskZ = false([d1,d2]);
    interpZ =  zeros([d1,d2],'int16');
    zP = gather(zP);
  end
  
  interpZ(linInd) = int16(round(zP));
  
  if ( any(zShift) )
    if strcmp(METHOD,'GPU')
      maskZ(linInd) = (zShift(1) - zShift(2) < interpZ(linInd) & ...
        zShift(1) + zShift(2) > interpZ(linInd));
    else
      maskZ(linInd) = gather((zShift(1) - zShift(2) < interpZ(linInd) & ...
        zShift(1) + zShift(2) > interpZ(linInd)));
    end
  else
    maskZ(linInd) = true;
  end
  
  evalMask(:,:,TLT(iPrj,1)) = maskZ;
  deltaZ(:,:,TLT(iPrj,1)) = interpZ;
  
end
clear X1 Y1 X2 Y2 Z2 rSample rElevation rTilt rInPlane R l x1 y1 xShift yShift
clear interpZ x3 y3 xNew yNew zNew linInd



end

