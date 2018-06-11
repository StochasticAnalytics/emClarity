%writeVTK       Write out the MRC stack in VTK Structured Points 
%
%   writeVTK(mRCImage, filename)
%
%   mRCImage    The MRCImage object.
%
%   filename    The name of the file to write.
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

function writeVTK(mRCImage, filename, dataTypeOut, scaleOffset)

if nargin < 3
  dataTypeOut = 'same';
end
if nargin < 4
  scaleOffset = 0;
end
  
fid = fopen(filename, 'w+');
if fid == -1
  disp(msg)
  PEETError(['Unable to open ' filename ' w+']);
end

fprintf(fid, '# vtk DataFile Version 3.0\n');
fprintf(fid, 'MRC Image stack\n');
fprintf(fid, 'BINARY\n');
fprintf(fid, 'DATASET STRUCTURED_POINTS\n');

fprintf(fid, 'DIMENSIONS %d %d %d\n', ...
        getNX(mRCImage), getNY(mRCImage), getNZ(mRCImage));
fprintf(fid, 'SPACING %f %f %f\n', ...
        getCellX(mRCImage) / getNX(mRCImage), ...
        getCellY(mRCImage) / getNY(mRCImage), ...
        getCellZ(mRCImage) /getNZ(mRCImage));

fprintf(fid, 'ORIGIN %f %f %f\n', ...
        getNXStart(mRCImage), getNYStart(mRCImage), getNZStart(mRCImage));

fprintf(fid, 'POINT_DATA %d\n', prod(getDimensions(mRCImage)));

%  Convert the mode into the correct string
switch  mRCImage.header.mode
 case 0
  dataType = 'unsigned_char';
 case 1
  dataType = 'unsigned_short';
 case 2
  dataType = 'float';
 otherwise
  PEETError(['Mode ' int2str(mRCImage.header.mode) 'unsupported!']);
end
if isempty(mRCImage.filename)
  mRCImage.filename = 'MRCData';
end
[path name] = fileparts(mRCImage.filename);

fprintf(fid, 'SCALARS %s %s 1\n', name, dataType);
fprintf(fid, 'LOOKUP_TABLE default\n');

% TODO: finish scaling and datatype conversion
% Fit the data to the expect range
if length(scaleOffset) == 1 && scaleOffset ~= 0
  header = getHeader(mRCImage);
  shift = header.minDensity
  range = header.maxDensity - header.minDensity
end
for iZ = 1:getNZ(mRCImage)
  im = getImage(mRCImage, iZ);
  im = im(:);
  
  fwrite(fid, im, getModeString(mRCImage));
end

fclose(fid);
