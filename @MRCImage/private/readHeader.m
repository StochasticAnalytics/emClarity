%readHeader     Read in the MRC file header
%
%   mRCImage = readHeader(mRCImage, debug)
%
%   mRCImage    The constructed MRCImage object
%
%   debug       OPTIONAL: Set to non-zero to get debugging output (default:0)
%
%   Read the header into the returned MRCImage object
%
%   Bugs: none known
%
% This file is part of PEET (Particle Estimation for Electron Tomography).
% Copyright 2000-2014 The Regents of the University of Colorado & BL3DEMC:
%           The Boulder Laboratory For 3D Electron Microscopy of Cells.
% See PEETCopyright.txt for more details.

% TODO:
% - implement Extra section correctly for BL3DFS
% - implement MRC reading

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  $Author: John Heumann $
%
%  $Date: 2014/01/13 20:00:38 $
%
%  $Revision: 6b413b88334c $
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mRCImage = readHeader(mRCImage, debug)

if nargin < 2
  debug = 0;
end
debugFD = 2;

if debug
  fprintf(debugFD, 'Reading first 3 variables\n');
end
  
%  Read in the dimensions of the data
mRCImage.header.nX = fread(mRCImage.fid, 1, 'int32');
mRCImage.header.nY = fread(mRCImage.fid, 1, 'int32');
mRCImage.header.nZ = fread(mRCImage.fid, 1, 'int32');

if debug
  fprintf(debugFD, 'Checking big versus little endian format...\n');
end
%  Check to see if the bytes need to be swapped, reopen the file in the
%  opposite format
if  mRCImage.header.nX < 0 || mRCImage.header.nY < 0 ||                ...
    mRCImage.header.nZ < 0 || (mRCImage.header.nX > 65535 &&           ...
    mRCImage.header.nY > 65535 && mRCImage.header.nZ > 65535)
  if debug
    fprintf(debugFD, 'Non-native endian format, swapping!\n');
  end
  [fname, perm, fileEndianFormat] = fopen(mRCImage.fid); %#ok<ASGLU>
  if strcmp('ieee-be', fileEndianFormat) == 1
    mRCImage.endianFormat = 'ieee-le';
elseif strcmp('ieee-be.l64', fileEndianFormat) == 1
    mRCImage.endianFormat = 'ieee-le.l64';
  elseif strcmp('ieee-le', fileEndianFormat) == 1
    mRCImage.endianFormat = 'ieee-be';
  elseif strcmp('ieee-le.l64', fileEndianFormat) == 1
    mRCImage.endianFormat = 'ieee-be.l64';
  else
    mRCImage.endianFormat = 'ieee-be';
  end
  if debug
    fprintf(debugFD, 'File format %s, attempting to open as %s\n',     ...
            fileEndianFormat, mRCImage.endianFormat);
  end
  fclose(mRCImage.fid);
  mRCImage.fid = fopen(mRCImage.filename, 'r', mRCImage.endianFormat);

  % reread the data dimensions in the correct endian format
  mRCImage.header.nX = fread(mRCImage.fid, 1, 'int32');
  mRCImage.header.nY = fread(mRCImage.fid, 1, 'int32');
  mRCImage.header.nZ = fread(mRCImage.fid, 1, 'int32');
else 
  if debug
    fprintf(debugFD, 'Native endian format.\n');
  end
end

if debug
  fprintf(debugFD, 'nX: %d\n', mRCImage.header.nX);
  fprintf(debugFD, 'nY: %d\n', mRCImage.header.nY);
  fprintf(debugFD, 'nZ: %d\n', mRCImage.header.nZ);
end

mRCImage.header.mode = fread(mRCImage.fid, 1, 'int32');
mRCImage.header.nXStart = fread(mRCImage.fid, 1, 'int32');
mRCImage.header.nYStart = fread(mRCImage.fid, 1, 'int32');
mRCImage.header.nZStart = fread(mRCImage.fid, 1, 'int32');
mRCImage.header.mX = fread(mRCImage.fid, 1, 'int32');
mRCImage.header.mY = fread(mRCImage.fid, 1, 'int32');
mRCImage.header.mZ = fread(mRCImage.fid, 1, 'int32');
mRCImage.header.cellDimensionX = fread(mRCImage.fid, 1, 'float32');
mRCImage.header.cellDimensionY = fread(mRCImage.fid, 1, 'float32');
mRCImage.header.cellDimensionZ = fread(mRCImage.fid, 1, 'float32');
mRCImage.header.cellAngleX = fread(mRCImage.fid, 1, 'float32');
mRCImage.header.cellAngleY = fread(mRCImage.fid, 1, 'float32');
mRCImage.header.cellAngleZ = fread(mRCImage.fid, 1, 'float32');
mRCImage.header.mapColumns = fread(mRCImage.fid, 1, 'int32');
mRCImage.header.mapRows = fread(mRCImage.fid, 1, 'int32');
mRCImage.header.mapSections = fread(mRCImage.fid, 1, 'int32');
mRCImage.header.minDensity = fread(mRCImage.fid, 1, 'float32');
mRCImage.header.maxDensity = fread(mRCImage.fid, 1, 'float32');
mRCImage.header.meanDensity = fread(mRCImage.fid, 1, 'float32');
mRCImage.header.spaceGroup = fread(mRCImage.fid, 1, 'int32');
% nBytesExtended is called nsymbt in the 2014 MRC Standard
mRCImage.header.nBytesExtended = fread(mRCImage.fid, 1, 'int32');
if debug
  fprintf(debugFD, 'nBytesExtended %d\n', mRCImage.header.nBytesExtended);
end

% MRC EXTRA section
mRCImage.header.creatorID = fread(mRCImage.fid, 1, 'int16');
mRCImage.header.extraInfo1 = fread(mRCImage.fid, 30, 'uchar');
mRCImage.header.nBytesPerSection = fread(mRCImage.fid, 1, 'int16');
if debug
  fprintf(debugFD, 'nBytesPerSection %d\n', mRCImage.header.nBytesPerSection);
end
mRCImage.header.serialEMType = fread(mRCImage.fid, 1, 'int16');
if debug
  fprintf(debugFD, 'serialEMType %d\n', mRCImage.header.serialEMType);
end
mRCImage.header.extraInfo2 = fread(mRCImage.fid, 20, 'uchar');
  
% IMOD stamp / flags for deciding whether to read / write mode 0 files
% as signed or unsigned bytes (added in PEET 1.8.0)
mRCImage.header.imodStamp = fread(mRCImage.fid, 1, 'int32');
mRCImage.header.imodFlags = fread(mRCImage.fid, 1, 'int32');

mRCImage.header.idtype = fread(mRCImage.fid, 1, 'int16');
mRCImage.header.lens = fread(mRCImage.fid, 1, 'int16');
mRCImage.header.ndl = fread(mRCImage.fid, 1, 'int16');
mRCImage.header.nd2 = fread(mRCImage.fid, 1, 'int16');
mRCImage.header.vdl = fread(mRCImage.fid, 1, 'int16');
mRCImage.header.vd2 = fread(mRCImage.fid, 1, 'int16');
mRCImage.header.tiltAngles(1) = fread(mRCImage.fid, 1, 'float32');
mRCImage.header.tiltAngles(2) = fread(mRCImage.fid, 1, 'float32');
mRCImage.header.tiltAngles(3) = fread(mRCImage.fid, 1, 'float32');
mRCImage.header.tiltAngles(4) = fread(mRCImage.fid, 1, 'float32');
mRCImage.header.tiltAngles(5) = fread(mRCImage.fid, 1, 'float32');
mRCImage.header.tiltAngles(6) = fread(mRCImage.fid, 1, 'float32');
mRCImage.header.xOrigin = fread(mRCImage.fid, 1, 'float32');
mRCImage.header.yOrigin = fread(mRCImage.fid, 1, 'float32');
mRCImage.header.zOrigin = fread(mRCImage.fid, 1, 'float32');
mRCImage.header.map = char(fread(mRCImage.fid, 4, 'uchar')); %#ok<*FREAD>
mRCImage.header.machineStamp = char(fread(mRCImage.fid, 4, 'uchar'));
mRCImage.header.densityRMS = fread(mRCImage.fid, 1, 'float32');

% Read in the label data and skip blank labels 
mRCImage.header.nLabels = fread(mRCImage.fid, 1, 'int32');
if debug
  fprintf(debugFD, 'Reading %d labels\n', mRCImage.header.nLabels);
end

for iLabel = 1:mRCImage.header.nLabels
  mRCImage.header.labels(iLabel,:) = fread(mRCImage.fid, 80, 'char');
  if debug
    fprintf(debugFD, '%s\n', mRCImage.header.labels(iLabel,:));
  end
end

for iJunk = mRCImage.header.nLabels+1:10
  junk = fread(mRCImage.fid, 80, 'uchar=>uchar');  %#ok<NASGU>
end

mRCImage.extended = ...
  fread(mRCImage.fid, mRCImage.header.nBytesExtended, 'uchar=>uchar');

% Set the image data starting point
mRCImage.dataIndex = 1024 + mRCImage.header.nBytesExtended;

if debug
  % display the nubmer of bytes in the file
  nFileBytes = getFileNBytes(mRCImage);
  fprintf(debugFD, 'File contains %d bytes\n', nFileBytes);
  fprintf(debugFD, 'MRC Header length: 1024\n');
  fprintf(debugFD, 'Extended header length: %d\n',                     ...
    mRCImage.header.nBytesExtended);
  dataBytes = mRCImage.header.nX  * mRCImage.header.nY                 ...
      * mRCImage.header.nZ * getModeBytes(mRCImage);
  fprintf(debugFD, 'Data bytes: %d\n', dataBytes);
  fprintf(debugFD, 'Difference: %d\n' , nFileBytes - 1024 - dataBytes - ...
          mRCImage.header.nBytesExtended);
end
