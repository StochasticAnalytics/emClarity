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
%  This has been significantly modified to reduce memory burden and speed
%  up reads for subvolumes. A subregion read is 10x faster on a small
%  volume, much more from a large volume by removing a loop over columns
%  and using the "skip" behavior in fread.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function vol = getVolume(mRCImage, varargin)

% TODO fix this cluster fck
flgCloseFile = true;
if nargin > 1

  if length(varargin) >= 1
    if strcmpi(varargin{1}, 'keep')
      flgCloseFile = false;
      iIndex = 1:getNX(mRCImage);
      iRange=1;
    else
      iRange = varargin{1};
      if isempty(iRange)
        iIndex = 1:getNX(mRCImage);
        iRange=1;
      else
        iRange = varargin{1};
        if length(iRange) == 1 && iRange ~= -1
          iIndex = iRange;
        elseif iRange == -1
          iIndex = 1:getNX(mRCImage);
        else
          iIndex = iRange(1):iRange(2);
        end
      end
    end 
  end   
 
  if length(varargin) >= 2
    if strcmpi(varargin{2}, 'keep')
      flgCloseFile = false;
      jIndex = 1:getNY(mRCImage);
      jRange=1;
    else
      jRange = varargin{2};
      if isempty(jRange)
        jIndex = 1:getNY(mRCImage);
        jRange=1;
      else
        jRange = varargin{2};
        if length(jRange) == 1 && jRange ~= -1
          jIndex = jRange;
        elseif jRange == -1
          jIndex = 1:getNY(mRCImage);
        else
          jIndex = jRange(1):jRange(2);
        end
      end
    end
  end

  if length(varargin) >= 3
    if strcmpi(varargin{3}, 'keep')
      flgCloseFile = false;
      kIndex = 1:getNZ(mRCImage);
      kRange = 1;
    else
      kRange = varargin{3};
      if isempty(kRange)
        kIndex = 1:getNZ(mRCImage);
        kRange = 1;
      else
        if length(kRange) == 1
          kIndex = kRange;
        else
          kIndex = kRange(1):kRange(2);
        end
      end
    end
  end


  if length(varargin) == 4
    if strcmpi(varargin{4}, 'keep')
      flgCloseFile = false;
    else
      error('Fourth arg to getVolume must be "keep" to keep MRCobject or nothing.')
    end
  end

else
      iIndex = 1:getNX(mRCImage);
      iRange=1;
      jIndex = 1:getNY(mRCImage);
      jRange=1;
      kIndex = 1:getNZ(mRCImage);
      kRange=1;
end




% If the volume is already loaded return the selected indices
if mRCImage.flgVolume

  % FIXME
  % I don't remember why this conditional was set, and it was probably
  % stupid. Commenting it out and watch for a break : /
  % It looks like replacing any -1 with [] in calls to getVolume should fix
  % ?
  
  
%   if length(iRange)==1 && length(jRange) ==1 && iRange ==-1 && jRange == -1
    vol = mRCImage.volume(iIndex, jIndex, kIndex);
%   else
%     vol = mRCImage.volume;
%   end
  % If the below is used, it copies the data, which is a huge bummer for full tomograms.
  %vol = mRCImage.volume(iIndex, jIndex, kIndex);
  return
end


modeStr = getModeString(mRCImage);



if strcmp(modeStr, 'int16*2') || strcmp(modeStr, 'float32*2')
  % Changed these to 1/0 rather than true/false so they can be used directly to 
  % scale the word length vi calc not bool - BAH 2017-11-22
  flgComplex = 1;
  modeStr = modeStr(1 : end - 2);
else
  flgComplex = 0;
end  
% fread doesn't yet recognize "half"
% read as uint16 and then typecast to half later
if strcmp(modeStr, 'half')
  modeStr = 'uint16';
end

% Changed counters and index variables to be clear to me BAH 2017-11-22
% l --> nSlice, k --> iSlice
% m --> nRow, j --> iCol
%
%  Walk through the images


% Number of elements garaunteed to be linearly continuos. - BAH 2017-11-22
nImageElements = length(iIndex);
nModeBytes = getModeBytes(mRCImage);
mode =  mRCImage.header.mode;
readBytesAsSigned = getReadBytesAsSigned(mRCImage);
nX = mRCImage.header.nX;
nY = mRCImage.header.nY;
nXYModeBytes = nX * nY * nModeBytes;

nPixelsBetween =  nX - nImageElements;
headerOffset = mRCImage.dataIndex;
yOffset = (jIndex(1) - 1)  * nX * nModeBytes;
xOffset = (iIndex(1) - 1 ) * nModeBytes;

wordLength = (flgComplex*1 + 1) .* [nImageElements];
flgReCast = 1;

readX = length(iIndex);
readY = length(jIndex);
imgSize = readX*readY;

switch mode
  % Make a string to tell fread how long each "value is. BAH 2017-11-22
  case 0
    precisionString = 'int8';
    nToRead = sprintf('%d*int8=>single',wordLength);
    nToSkip = 1*nPixelsBetween * (flgComplex*1 + 1);
    %  Allocate the output matrix - NOTE: always single precision at end
    % 2x as fast to recast the whole array, than to either read in as
    % single or to recast each slice.

    
  case 1
    precisionString = 'int16';
    nToRead = sprintf('%d*int16=>int16',wordLength);
    nToSkip = 2*nPixelsBetween * (flgComplex*1 + 1);
    %  Allocate the output matrix - NOTE: always single precision at end
    % 2x as fast to recast the whole array, than to either read in as
    % single or to recast each slice.
    
  case 2
     precisionString = 'single';
    flgReCast = 0;
    nToRead = sprintf('%d*single=>single',wordLength(1));
    nToSkip = 4*nPixelsBetween * (flgComplex*1 + 1);
     %  Allocate the output matrix - NOTE: always single precision at end
    % 2x as fast to recast the whole array, than to either read in as
    % single or to recast each slice.
  case 6
    precisionString = 'uint16';
    nToRead = sprintf('%d*uint16=>uint16',wordLength);
    nToSkip = 2*nPixelsBetween * (flgComplex*1 + 1);

  case 12
    precisionString = 'uint16';
    nToRead = sprintf('%d*uint16=>uint16',wordLength);
    nToSkip = 2*nPixelsBetween * (flgComplex*1 + 1);
  otherwise
    error('did not recognize mode value %d\n', mode)
end

%  Allocate the output matrix - NOTE: always single precision at end
% 2x as fast to recast the whole array, than to either read in as
% single or to recast each slice.
if (flgComplex)
  for iVol = 1:2
    vol{iVol} = zeros(length(iIndex), length(jIndex), length(kIndex), precisionString);
  end
else
  vol = zeros(length(iIndex), length(jIndex), length(kIndex), precisionString);
end

nSlice = 1;

for iSlice = kIndex
  if iSlice < 1 || iSlice > mRCImage.header.nZ
    PEETError('Image index out of range');
  end

  idxSectionStart = ((iSlice - 1) * nXYModeBytes) + xOffset + yOffset + headerOffset; 

  % Move the file pointer to the next slice  BAH 2017-11-22
  fseek(mRCImage.fid, idxSectionStart, 'bof');
  [img, iCount] = fread(mRCImage.fid, [readX, readY], nToRead, nToSkip);

  if iCount ~= imgSize
    PEETError(['Expected ' int2str(imgSize) ' elements, read ' int2str(iCount)]);
  end

  if ( flgComplex )
    vol{1}(:,:,nSlice)  = reshape(img(1:2:end),readX,readY);
    vol{2}(:,:,nSlice)  = reshape(img(2:2:end),readX,readY);
  else
    vol(:,:,nSlice) = (img);  % BAH 2017-11-22
  end
  nSlice = nSlice + 1;
end


% No need to do this within the read loop BAH 2017-11-22
%if mode == 0 && readBytesAsSigned
  % We just read a byte image as unsigned (which was the pre 1.8.0
  % convention), but signed bytes appear to be what was intended. Remap
  % values accordingly to preserve the correct ordering. In each case,
  % img will be represented internally as unsigned bytes 0..255.
%  topHalf = vol >= 0 & vol < 127;
%  vol(topHalf) = vol(topHalf) + 128;
%  vol(~topHalf) = vol(~topHalf) - 128;
%end

if ( flgReCast )
  if ( flgComplex )
    vol = complex(single(vol{1}),single(vol{2}));
  else
    if mRCImage.header.mode == 12
      vol = half.typecast(vol);
    else
      vol = single(vol);
    end
  end
else
  if ( flgComplex )
    vol = complex(vol{1},vol{2});
  end
end

% When dealing with many subtomograms, too many files are open (OS dependent)
% always close MRCImage object unless specified.
if (flgCloseFile)
  fclose(mRCImage.fid);
  mRCImage.fid = [];
end



end
  
