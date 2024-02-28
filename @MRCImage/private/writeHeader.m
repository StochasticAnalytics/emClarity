%writeHeader    Write out the MRC file header
%
%   mRCImage = writeHeader(mRCImage, debug)
%
%   mRCImage    The MRCImage object
%
%   debug       OPTIONAL: Set to non-zero to get debugging output (default:0)
%
%   Write out the header contained in the returned MRCImage object to the file
%   specified by that object.
%
%   Calls: none
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
%  $Date: 2014/02/10 18:30:37 $
%
%  $Revision: c1b2c05695d8 $
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mRCImage = writeHeader(mRCImage, debug)

if nargin < 2
  debug = 0;
end
debugFD = 2;

if debug
  fprintf(debugFD, 'Writing first 3 variables\n');
end

% Check the permissions on the file ID, reopen it if needed

mRCImage.fid = openWritable(mRCImage);

% Return to the beginning of the file
status = fseek(mRCImage.fid, 0, 'bof');
if status
  disp('Could not move the file pointer to the begining ');
  PEETError('Could not seek to beginning of file id %d', mRCImage.fid);
end

% Write out the dimensions of the data
writeAndCheck(mRCImage.fid, mRCImage.header.nX, 'int32');
writeAndCheck(mRCImage.fid, mRCImage.header.nY, 'int32');
writeAndCheck(mRCImage.fid, mRCImage.header.nZ, 'int32');
writeAndCheck(mRCImage.fid, mRCImage.header.mode, 'int32');
writeAndCheck(mRCImage.fid, mRCImage.header.nXStart, 'int32');
writeAndCheck(mRCImage.fid, mRCImage.header.nYStart, 'int32');
writeAndCheck(mRCImage.fid, mRCImage.header.nZStart, 'int32');
writeAndCheck(mRCImage.fid, mRCImage.header.mX, 'int32');
writeAndCheck(mRCImage.fid, mRCImage.header.mY, 'int32');
writeAndCheck(mRCImage.fid, mRCImage.header.mZ, 'int32');
writeAndCheck(mRCImage.fid, mRCImage.header.cellDimensionX, 'float32');
writeAndCheck(mRCImage.fid, mRCImage.header.cellDimensionY, 'float32');
writeAndCheck(mRCImage.fid, mRCImage.header.cellDimensionZ, 'float32');
writeAndCheck(mRCImage.fid, mRCImage.header.cellAngleX, 'float32');
writeAndCheck(mRCImage.fid, mRCImage.header.cellAngleY, 'float32');
writeAndCheck(mRCImage.fid, mRCImage.header.cellAngleZ, 'float32');
writeAndCheck(mRCImage.fid, mRCImage.header.mapColumns, 'int32');
writeAndCheck(mRCImage.fid, mRCImage.header.mapRows, 'int32');
writeAndCheck(mRCImage.fid, mRCImage.header.mapSections, 'int32');
writeAndCheck(mRCImage.fid, mRCImage.header.minDensity, 'float32');
writeAndCheck(mRCImage.fid, mRCImage.header.maxDensity, 'float32');
writeAndCheck(mRCImage.fid, mRCImage.header.meanDensity, 'float32');

writeAndCheck(mRCImage.fid, mRCImage.header.spaceGroup, 'int32');
writeAndCheck(mRCImage.fid, mRCImage.header.nBytesExtended, 'int32');
% MRC EXTRA section
writeAndCheck(mRCImage.fid, mRCImage.header.creatorID, 'int16');
writeAndCheck(mRCImage.fid, mRCImage.header.extraInfo1, 'uchar');
  
writeAndCheck(mRCImage.fid, mRCImage.header.nBytesPerSection, 'int16');
writeAndCheck(mRCImage.fid, mRCImage.header.serialEMType, 'int16');
writeAndCheck(mRCImage.fid, mRCImage.header.extraInfo2, 'uchar');
  
mRCImage.header.imodStamp = defaultIMODStamp();
writeAndCheck(mRCImage.fid, mRCImage.header.imodStamp, 'int32');
if getWriteBytesAsSigned(mRCImage)
  mRCImage.header.imodFlags =                                          ...
    int32(bitor(uint32(mRCImage.header.imodFlags), 1));
end
writeAndCheck(mRCImage.fid, mRCImage.header.imodFlags, 'int32');

writeAndCheck(mRCImage.fid, mRCImage.header.idtype, 'int16');
writeAndCheck(mRCImage.fid, mRCImage.header.lens, 'int16');
writeAndCheck(mRCImage.fid, mRCImage.header.ndl, 'int16');
writeAndCheck(mRCImage.fid, mRCImage.header.nd2, 'int16');
writeAndCheck(mRCImage.fid, mRCImage.header.vdl, 'int16');
writeAndCheck(mRCImage.fid, mRCImage.header.vd2, 'int16');
writeAndCheck(mRCImage.fid, mRCImage.header.tiltAngles, 'float32');

% Now support only "newer" BL3DFS format
writeAndCheck(mRCImage.fid, mRCImage.header.xOrigin, 'float32');
writeAndCheck(mRCImage.fid, mRCImage.header.yOrigin, 'float32');
writeAndCheck(mRCImage.fid, mRCImage.header.zOrigin, 'float32');
writeAndCheck(mRCImage.fid, mRCImage.header.map, 'uchar');
writeAndCheck(mRCImage.fid, mRCImage.header.machineStamp, 'uchar');
writeAndCheck(mRCImage.fid, mRCImage.header.densityRMS, 'float32');

% Write out the label data and skip blank labels 
if mRCImage.header.nLabels < 10
  writeAndCheck(mRCImage.fid, mRCImage.header.nLabels + 1, 'int32');
else
  writeAndCheck(mRCImage.fid, mRCImage.header.nLabels, 'int32');
end

% Write any pre-existing labels
for iLabel = 1:mRCImage.header.nLabels
  writeAndCheck(mRCImage.fid, mRCImage.header.labels(iLabel,:), 'uchar');
end

% If there's room, add a label indicating writing by PEET
if mRCImage.header.nLabels < 10
msg = ['Written by PEET / MatTomo ' datestr(now)];
writeAndCheck(mRCImage.fid, msg, 'uchar');
writeAndCheck(mRCImage.fid, char(blanks(80 - length(msg))), 'uchar');
mRCImage.header.nLabels = mRCImage.header.nLabels + 1;
end

% Use blank messages for the remainder
for iJunk = mRCImage.header.nLabels+1:10
  writeAndCheck(mRCImage.fid, char(zeros(1, 80)), 'uchar');
end


writeAndCheck(mRCImage.fid, mRCImage.extended, 'uchar');

return

% Simple error checking write
function writeAndCheck(fid, matrix, precision)
  nElements = numel(matrix);
  if strcmp(precision, 'half')
    count = fwrite(fid, matrix, 'uint16');
  else
    count = fwrite(fid, matrix, precision);
  end
  if count ~= nElements
    error('Matrix contains %d elements, but only wrote %d',nElements, count);
  end
return
