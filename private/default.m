%default        Default MRCImage file structure
%
%   mRCImage = default
%
%   mRCImage    The MRCImage object structure.
%
%   Bugs: none known
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

function mRCImage = default

mRCImage.fid = [];
mRCImage.filename = [];
mRCImage.endianFormat = 'ieee-le';
mRCImage.type = 'BL3DFS';
mRCImage.version = '1.0';

mRCImage.dataIndex = -Inf;
mRCImage.volume = [];
mRCImage.flgVolume = 0;

mRCImage.header.nX = -Inf;
mRCImage.header.nY = -Inf;
mRCImage.header.nZ = -Inf;
mRCImage.header.mode = -Inf;
mRCImage.header.nXStart = -Inf;
mRCImage.header.nYStart = -Inf;
mRCImage.header.nZStart = -Inf;
mRCImage.header.mX = -Inf;
mRCImage.header.mY = -Inf;
mRCImage.header.mZ = -Inf;
mRCImage.header.cellDimensionX = -Inf;
mRCImage.header.cellDimensionY = -Inf;
mRCImage.header.cellDimensionZ = -Inf;
mRCImage.header.cellAngleX = -Inf;
mRCImage.header.cellAngleY = -Inf;
mRCImage.header.cellAngleZ = -Inf;
mRCImage.header.mapColumns = 1;
mRCImage.header.mapRows = 2;
mRCImage.header.mapSections = 3;
mRCImage.header.minDensity = -Inf;
mRCImage.header.maxDensity = -Inf;
mRCImage.header.meanDensity = -Inf;
mRCImage.header.spaceGroup = -Inf;
mRCImage.header.nSymmetryBytes = -Inf;
mRCImage.header.nBytesExtended = -Inf;
mRCImage.header.creatorID = -Inf;
mRCImage.header.nBytesPerSection = -Inf;
mRCImage.header.serialEMType = -Inf;
mRCImage.header.imodStamp = -Inf;
mRCImage.header.imodFlags = -Inf;
mRCImage.header.idtype = -Inf;
mRCImage.header.lens = -Inf;
mRCImage.header.ndl = -Inf;
mRCImage.header.nd2 = -Inf;
mRCImage.header.vdl = -Inf;
mRCImage.header.vd2 = -Inf;
mRCImage.header.tiltAngles = [];
mRCImage.header.extra = [];
mRCImage.header.xOrigin = -Inf;
mRCImage.header.yOrigin = -Inf;
mRCImage.header.zOrigin = -Inf;
mRCImage.header.map = '';
mRCImage.header.machineStamp = '';
mRCImage.header.densityRMS = -Inf;
mRCImage.header.nLabels = -Inf;
mRCImage.header.labels = [blanks(80); blanks(80); 
                    blanks(80); blanks(80); 
                    blanks(80); blanks(80);
                    blanks(80); blanks(80);
                    blanks(80); blanks(80) ];

mRCImage.extended = [];

mRCImage.forceWriteByteMode = [];

return 
