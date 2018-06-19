%zero           Zero the volume data in a MRCImage
%
%   mRCImage = zero(mRCImage)
%
%   mRCImage    The MRCImage object.
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

function mRCImage = zero(mRCImage)

modeStr = getModeString(mRCImage);
if strcmp(modeStr, 'int16*2') || strcmp(modeStr, 'float32*2')
  flgComplex = true;
  modeStr = modeStr(1 : end - 2);
else
  flgComplex = false;
end

%  Zero the volume data
if mRCImage.flgVolume
  z = zeros(mRCImage.header.nX, mRCImage.header.nY, mRCImage.header.nZ);
  if flgComplex
    z = complex(z, z);
  end
  if strcmp(modeStr, 'float32')
    mRCImage.volume = single(z);
  else
    mRCImage.volume = feval(modeStr, z);
  end
else
  mRCImage.volume = [];
  % Check the permissions on the file ID, reopen it if needed
  mRCImage.fid = openWritable(mRCImage);
  
% % %   %  Move the data pointer to the start of the image data section of the
% % %   %  file
% % %   status = fseek(mRCImage.fid, mRCImage.dataIndex, 'bof');
% % %   if status
% % %     PEETError(ferror(mRCImage.fid))
% % %   end
% % % 
% % %   % Write out the correct number of zeros in image size chunks
% % %   nImageElements = mRCImage.header.nX * mRCImage.header.nY;
% % %   if flgComplex
% % %     img = zeros(2 * nImageElements, 1);
% % %     for iSlice = 1:mRCImage.header.nZ
% % %       count = fwrite(mRCImage.fid, img, modeStr);
% % %       if count ~= 2 * nImageElements
% % %         PEETError(['Expected ' int2str(nImageElements) ' elements, wrote ' ...
% % %           int2str(count / 2) '\n  ' ferror(mRCImage.fid)]);
% % %       end
% % %     end
% % %   else
% % %     img = zeros(nImageElements, 1);
% % %     for iSlice = 1:mRCImage.header.nZ
% % %       count = fwrite(mRCImage.fid, img, modeStr);
% % %       if count ~= nImageElements
% % %         PEETError(['Expected ' int2str(nImageElements) ' elements, '   ...
% % %           'wrote ' int2str(count) '\n  ' ferror(mRCImage.fid)]);
% % %       end
% % %     end
% % %   end
end

