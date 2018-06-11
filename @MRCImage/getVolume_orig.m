%getVolume      Extract a volume from the MRC image
%
%   vol = getVolume(mRCImage, iRange, jRange, kRange)
%
%   vol         The extracted volume
%
%   mRCImage    The opened MRCImage object
%
%   iRange      The indices of the i (first, columns) dimension to be
%               extracted as [iMin iMax], an empty array specifies all.
%
%   jRange      The indices of the j (second, rows) dimension to be
%               extracted as [jMin jMax], an empty array specifies all.
%
%   kRange      The indices of the k (third, planes) dimension to be
%               extracted as [kMin kMax], an empty array specifies all.
%
%   MRCImage.getVolume extracts a three dimensional region from an MRCImage
%   object.
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

function vol = getVolume(mRCImage, iRange, jRange, kRange)

if nargin < 2 || isempty(iRange)
  iIndex = 1:getNX(mRCImage);
  iRange=1;
else
  if length(iRange) == 1 && iRange ~= -1
    iIndex = iRange;
  elseif iRange == -1
    iIndex = 1:getNX(mRCImage);
  else
    iIndex = iRange(1):iRange(2);
  end
end

if nargin < 3 || isempty(jRange)
  jIndex = 1:getNY(mRCImage);
    jRange=1;
else
  if length(jRange) == 1 && jRange ~= -1
    jIndex = jRange;
  elseif jRange == -1
    jIndex = 1:getNY(mRCImage);
  else
    jIndex = jRange(1):jRange(2);
  end
end

if nargin < 4 || isempty(kRange)
  kIndex = 1:getNZ(mRCImage);
else
  if length(kRange) == 1
    kIndex = kRange;
  else
    kIndex = kRange(1):kRange(2);
  end
end

% If the volume is already loaded return the selected indices
if mRCImage.flgVolume

  if length(iRange)==1 && length(jRange) ==1 && iRange ==-1 && jRange == -1

    vol = mRCImage.volume(iIndex, jIndex, kIndex);
  else
    vol = mRCImage.volume;
  end
  % If the below is used, it copies the data, which is a huge bummer for full tomograms.
  %vol = mRCImage.volume(iIndex, jIndex, kIndex);
  return
end


modeStr = getModeString(mRCImage);

%  Allocate the output matrix - NOTE: always single precision
vol = zeros(length(iIndex), length(jIndex), length(kIndex), 'single');
if strcmp(modeStr, 'int16*2') || strcmp(modeStr, 'float32*2')
  vol = complex(vol, vol);
  flgComplex = true;
  modeStr = modeStr(1 : end - 2);
else
  flgComplex = false;
end  

%  Walk through the images
l = 1;
precision = [modeStr '=>' modeStr];
nImageElements = length(iIndex);
nModeBytes = getModeBytes(mRCImage);
mode =  mRCImage.header.mode;
readBytesAsSigned = getReadBytesAsSigned(mRCImage);
nX = mRCImage.header.nX;
nY = mRCImage.header.nY;
nXYModeBytes = nX * nY * nModeBytes;
for k = kIndex
  if k < 1 || k > mRCImage.header.nZ
    PEETError('Image index out of range');
  end
  idxSectionStart = mRCImage.dataIndex + (k - 1) * nXYModeBytes;

  m=1;
  for j = jIndex
    idxDataStart = idxSectionStart + ((j - 1) * nX +                   ...
      iIndex(1) - 1 )* nModeBytes;
    if (flgComplex)   % handle reading complex image
      % Move the file pointer to the next line
      idxDataStart = idxSectionStart + ((j - 1) * nX +                 ...
        iIndex(1) - 1 ) * nModeBytes;
      fseek(mRCImage.fid, idxDataStart, 'bof');
      
      % Read in the line from the MRC file
      [temp, count] = fread(mRCImage.fid, 2 * nImageElements, precision);
      if count ~= 2 * nImageElements
        PEETError(['Expected ' int2str(nImageElements)                 ...
          ' elements, read ' int2str(count / 2)]);
      end
      img = complex(temp(1:2:end-1), temp(2:2:end));      
    else              % normal (not complex) image
      % Move the file pointer to the next line
      fseek(mRCImage.fid, idxDataStart, 'bof');
      
      % Read in the line from the MRC file
      [img, count] = fread(mRCImage.fid, nImageElements, precision);
      if count ~= nImageElements
        PEETError(['Expected ' int2str(nImageElements)                 ...
          ' elements, read ' int2str(count)]);
      end
    end
    
    if mode == 0 && readBytesAsSigned
      % We just read a byte image as unsigned (which was the pre 1.8.0
      % convention), but signed bytes appear to be what was intended. Remap
      % values accordingly to preserve the correct ordering. In any case,
      % img will be represented internally as unsigned bytes 0..255.
      topHalf = img >= 0 & img < 127;
      img(topHalf) = img(topHalf) + 128;
      img(~topHalf) = img(~topHalf) - 128;
      mRCImage = setStatisticsFromVolume(mRCImage);
    end
  
    vol(:, m, l) = img;
    m = m + 1; % next line
  end
  l = l + 1; %next slice
end
