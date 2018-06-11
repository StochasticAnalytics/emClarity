%writeEM         Write the file in EM format
%
%   writeEM(mRCImage, filename)
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

function writeEM(mRCImage, filename)

fid = fopen(filename, 'w+');
if fid == -1
  disp(msg)
  PEETError(['Unable to open ' filename ' w+']);
end
%
% Write the 512 byte EM header
%
count = fwrite(fid, [6 0 0], 'char');
if count ~= 3
  PEETError(['Unable to write to : ' filename]);
end
switch  mRCImage.header.mode
 case 0
  count = fwrite(fid, 1, 'char');
 case 1
  count = fwrite(fid, 2, 'char');
 case 2
  count = fwrite(fid, 5, 'char');
 case 4
  count = fwrite(fid, 8, 'char');
 otherwise
  PEETError(['Unimplemented mode in MRCImage: ' getModeString(mRCImage)]);
end
if count ~= 1
  PEETError(['Unable to write to : ' filename]);
end

count = fwrite(fid, [mRCImage.header.nX mRCImage.header.nY             ...
                     mRCImage.header.nZ], 'int32');
if count ~= 3
  PEETError(['Unable to write to : ' filename]);
end

%  Comment section
count = fwrite(fid, zeros(80,1), 'char');
if count ~= 80
  PEETError(['Unable to write to : ' filename]);
end

% User defined section
count = fwrite(fid, zeros(40,1), 'int32');
if count ~= 40
  PEETError(['Unable to write to : ' filename]);
end

% User data section
count = fwrite(fid, zeros(256,1), 'char');
if count ~= 256
  PEETError(['Unable to write to : ' filename]);
end


% Write out the data now with X dimension as the fastest, then Y, then Z
nElem = mRCImage.header.nX * mRCImage.header.nY;
modeString = getModeString(mRCImage);
for iZ = 1:mRCImage.header.nZ
  count = fwrite(fid, getImage(mRCImage, iZ), modeString);
  if count ~= nElem
    PEETError(['problem writing data to : ' filename]);
  end
end
