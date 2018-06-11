%getImage       Get the specified image out of the MRCImage
%
%   img = getImage(mRCImage, index)
%
%   img         The specified image as a 2D array
%
%   mRCImage    The MRCImage object containing the stack of interest
%
%   index       The index of the image to extract, indices start at 1.
%
%
%   Return the specifed image from the MRCImage object.
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

function img = getImage(mRCImage, index)

if index < 1 || index > mRCImage.header.nZ
  PEETError('Image index out of range!');
end

index = double(index); % avoid overflow in idxDataStart computation

if mRCImage.flgVolume
  img = mRCImage.volume(:,:, index);
else
  modeStr = getModeString(mRCImage);
  if strcmp(modeStr, 'int16*2') || strcmp(modeStr, 'float32*2')
    modeStr = modeStr(1 : end - 2);
    flgComplex = true;
  else
    flgComplex = false;
  end
  nImageElements = mRCImage.header.nX * mRCImage.header.nY;
  
  % Move the file pointer to the requested image
  idxDataStart = mRCImage.dataIndex + (index - 1) * nImageElements *   ...
    getModeBytes(mRCImage);
  precision = [modeStr '=>' modeStr];
  fseek(mRCImage.fid, idxDataStart, 'bof');
  
  if flgComplex      % handle reading a complex image
    % Read in the image from the MRC file
    [temp, count] = fread(mRCImage.fid, 2 * nImageElements, precision);
    if count ~= 2 * nImageElements
      PEETError(['Expected ' int2str(nImageElements) ' elements, read ' ...
        int2str(count / 2)]);
    end
    img = complex(temp(1:2:end-1), temp(2:2:end));
  else               % normal (not complex) image
    % Read in the image from the MRC file
    [img, count] = fread(mRCImage.fid, nImageElements, precision);
    if count ~= nImageElements
      PEETError(['Expected ' int2str(nImageElements) ' elements, read ' ...
        int2str(count)]);
    end
  end
  
  if mRCImage.header.mode == 0 && getReadBytesAsSigned(mRCImage)
    % We just read a byte image as unsigned (which was the pre 1.8.0
    % convention), but signed bytes appear to be what was intended. Remap
    % values accordingly to preserve the correct ordering. Note that in any
    % case, img will be represented internally as unsigned bytes 0..255.
    topHalf = img >= 0 & img < 127;
    img(topHalf) = img(topHalf) + 128;
    img(~topHalf) = img(~topHalf) - 128;
    mRCImage = setStatisticsFromVolume(mRCImage);
  end

  % Reshape the image to match nX and nY
  img = reshape(img, mRCImage.header.nX, mRCImage.header.nY);
end
