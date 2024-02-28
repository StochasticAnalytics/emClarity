% setVolumeAndHeaderFromVolume
%
% This file is part of PEET (Particle Estimation for Electron Tomography).
% Copyright 2000-2014 The Regents of the University of Colorado & BL3DEMC:
%           The Boulder Laboratory For 3D Electron Microscopy of Cells.
% See PEETCopyright.txt for more details.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  $Author: John Heumann $
%
%  $Date: 2014/01/13 20:00:38 $
%
%  $Revision: 6b413b88334c $
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mRCImage = setVolumeAndHeaderFromVolume(mRCImage, volume)

if isa(volume, 'uint8')
  mRCImage.header.mode = 0;
elseif isa(volume, 'int16')
  if isreal(volume)
    mRCImage.header.mode = 1;
  else
    mRCImage.header.mode = 3;
  end
elseif isa(volume, 'single')
  if isreal(volume)
    mRCImage.header.mode = 2;
  else
    mRCImage.header.mode = 4;
  end
elseif isa(volume, 'double')
  if ~isreal(volume)
    PEETError('Double precision complex volumes are not supported!');
  end
  mRCImage.header.mode = 2;
elseif isa(volume, 'uint16')
  if isreal(volume)
    mRCImage.header.mode = 6;
  else
    PEETError('Complex uint16 volumes are not supported!');
  end
elseif isa(volume, 'half')
  if isreal(volume)
    mRCImage.header.mode = 12;
  else
    PEETError('Complex half-precision volumes are not supported!');
  end
else
  PEETError('Volume must be uint8, int16, unit16, half, single, or double!');
end


if isa(volume, 'double')
  mRCImage.volume = single(volume); % Store float*64 as float*32 
else
  mRCImage.volume = volume;
end
%  TODO: is this the correct/best way to set the values
mRCImage.header.nX = size(mRCImage.volume, 1);
mRCImage.header.nY = size(mRCImage.volume, 2);
mRCImage.header.nZ = size(mRCImage.volume, 3);
mRCImage.header.nXStart = 0;
mRCImage.header.nYStart = 0;
mRCImage.header.nZStart = 0;
mRCImage.header.mX = size(mRCImage.volume, 1);
mRCImage.header.mY = size(mRCImage.volume, 2);
mRCImage.header.mZ = size(mRCImage.volume, 3);
mRCImage.header.cellDimensionX = size(mRCImage.volume, 1);
mRCImage.header.cellDimensionY = size(mRCImage.volume, 2);
mRCImage.header.cellDimensionZ = size(mRCImage.volume, 3);
mRCImage.header.cellAngleX = 90;
mRCImage.header.cellAngleY = 90;
mRCImage.header.cellAngleZ = 90;

mRCImage = setStatisticsFromVolume(mRCImage);

mRCImage.type = 'BL3DFS';
mRCImage.header.spaceGroup = 1;
mRCImage.header.nSymmetryBytes = 0;

mRCImage.header.imodStamp = defaultIMODStamp();
writeBytesAsSigned = getWriteBytesAsSigned(mRCImage);
mRCImage.header.imodFlags = writeBytesAsSigned;

mRCImage.header.creatorID = int16(0);
mRCImage.header.extraInfo1 = char(zeros(1, 30, 'uint8'));
mRCImage.header.nBytesPerSection = int16(0);
mRCImage.header.serialEMType = int16(0);
mRCImage.header.extraInfo2 = char(zeros(1, 20, 'uint8'));
mRCImage.header.idtype = 0;
mRCImage.header.lens = 0;
mRCImage.header.ndl = 0;
mRCImage.header.nd2 = 0;
mRCImage.header.vdl = 0;
mRCImage.header.vd2 = 0;
mRCImage.header.tiltAngles = single([0 0 0 0 0 0]);

mRCImage.header.extra = zeros(100, 1);
mRCImage.header.xOrigin = 0;
mRCImage.header.yOrigin = 0;
mRCImage.header.zOrigin = 0;
mRCImage.header.map = 'MAP ';
mRCImage.header.machineStamp = uint8([68 65 0 0]); % assume little endian
%mRCImage.header.machineStamp = uchar([17 17 0 0]); % assume big endian
  
mRCImage.flgVolume = 1;
mRCImage.header.nBytesExtended = 0;
mRCImage.header.nLabels = 0;

mRCImage.extended = [];
mRCImage.dataIndex = 1024;
