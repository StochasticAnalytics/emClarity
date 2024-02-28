%save           Save the MRCImage object
%
%   mRCImage = save(mRCImage, newFilename);
%
%   mRCImage    The MRCImage object
%
%   newFilename OPTIONAL: Change the filename associated with this object and
%               save it to the new filename.
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

%  There is no need (within BAH usage) to return mRCImage, and if requested
%  this binds memory that I can't seem to otherwise free up.
function [ ] = SAVE_IMG(varargin)
%function mRCImage = SAVE_IMG(varargin)

%args
% mRCImage, newFilename, pixelSize

fixHeader = false;
mRCImage = varargin{1};
  
if nargin > 1
  newFilename = varargin{2};

  mRCImage = close(mRCImage);

  mRCImage.filename = newFilename;
  % Check to see if the file exists
  result = dir(mRCImage.filename);
  if size(result, 1) > 0
    if result.isdir
      PEETError('%s is a directory!', newFilename);
    else
      % if the file exists it needs to move out of the way
      delete(mRCImage.filename);
    end
  end
end

% Bah add to optionally include pixel size.
if nargin > 2
  pixelSize = varargin{3};
  if numel(pixelSize) == 1
    pixelSize = pixelSize .* [1,1,1];
  end

 % Also set origin to -nx/2 which results in any xforms from an alignment in 
 % Chimera be used directly.
  if nargin > 3
    if numel(varargin{4}) == 1
      % center the origin
      mRCImage.header.xOrigin = -1*mRCImage.header.cellDimensionX/2;
      mRCImage.header.yOrigin = -1*mRCImage.header.cellDimensionY/2;
      mRCImage.header.zOrigin = -1*mRCImage.header.cellDimensionZ/2;
    else
      mRCImage.header.xOrigin = varargin{4}(1);
      mRCImage.header.yOrigin = varargin{4}(2);
      mRCImage.header.zOrigin = varargin{4}(3);  
    end
  end
  % Note that the following assume cellDimension == array size from matlab
  mRCImage.header.cellDimensionX = mRCImage.header.cellDimensionX * pixelSize(1);
  mRCImage.header.cellDimensionY = mRCImage.header.cellDimensionY * pixelSize(2);
  mRCImage.header.cellDimensionZ = mRCImage.header.cellDimensionZ * pixelSize(3);

 
end
% Write out the header
mRCImage = writeHeader(mRCImage);

% Write out the volume if it is not already on the disk
if mRCImage.flgVolume
  modeStr = getModeString(mRCImage);
  if strcmp(modeStr, 'int16*2') || strcmp(modeStr, 'float32*2')
    modeStr = modeStr(1 : end - 2);
    flgComplex = true;
  else
    flgComplex = false;
  end
  % if isa(mRCImage.volume, 'half')
  %   mrcImage.volume = typecast(mRCImage.volume, 'uint16');
  % end
  if strcmp(modeStr, 'half')
    modeStr = 'uint16';
  end
  nElements = numel(mRCImage.volume);
  
  if mRCImage.header.mode == 0 && getWriteBytesAsSigned(mRCImage)
    % Volume is represented internally as unsigned (0..255).
    % Reorder to preserve correct ordering if interpreted as signed, and
    % re-write the header (with shifted density statistics).
    topHalf = mRCImage.volume >= 128;
    mRCImage.volume(topHalf) = mRCImage.volume(topHalf) - 128;
    mRCImage.volume(~topHalf) = mRCImage.volume(~topHalf) + 128;
    mRCImage.header.minDensity = mRCImage.header.minDensity - 128;
    mRCImage.header.maxDensity = mRCImage.header.maxDensity - 128;
    mRCImage.header.meanDensity = mRCImage.header.meanDensity - 128;
    mRCImage = writeHeader(mRCImage);
  end
  
  % Handle writing complex arrays when requested  
  if flgComplex
    nX = mRCImage.header.nX;
    nY = mRCImage.header.nY;
    nZ = mRCImage.header.nZ;
    if modeStr(1) == 'i'
      temp = int16(zeros(2*nX, nY, nZ));
    else
      temp = single(zeros(2*nX, nY, nZ));
    end
    temp(1:2:2*nX-1, :, :) = real(mRCImage.volume);
    temp(2:2:2*nX, :, :) = imag(mRCImage.volume);
    count = fwrite(mRCImage.fid, temp, modeStr);
    if count ~= 2 * nElements
      fprintf('Matrix contains %d but only wrote %d elements\n', ...
              nElements, count / 2);
      PEETError('Failed writing matrix!');
    end
  else % normal (not complex) data

    if (mRCImage.header.minDensity == 0.0 && mRCImage.header.maxDensity ==0 )
      if numel(mRCImage.volume) < 768^3
        mRCImage.header.minDensity = min(min(min(mRCImage.volume)));
        mRCImage.header.maxDensity = max(max(max(mRCImage.volume)));
    
        mRCImage.header.meanDensity = mean(mean(mean(single(mRCImage.volume))));
        mRCImage.header.densityRMS = std(std(std(single(mRCImage.volume))));
        mRCImage = writeHeader(mRCImage);
        fixHeader = 0;
      else
        fixHeader = 1;
      end
    end
    
    count = fwrite(mRCImage.fid, mRCImage.volume, modeStr);
    if count ~= nElements
      fprintf('Matrix contains %d but only wrote %d elements\n', ...
              nElements, count);
      PEETError('Failed writing matrix!');
    end
  end
end

% Set the cell angles to 90 and space group to 1 (vol not stack) by default
  mRCImage.header.cellAngleX = 90.0;
  mRCImage.header.cellAngleY = 90.0;
  mRCImage.header.cellAngleZ = 90.0;
  mRCImage.header.spaceGroup = 1;
  mRCImage = writeHeader(mRCImage);
  
%clear mRCImage.volume;
mRCImage.volume = [];
close(mRCImage);

%% TODO this is really slow, like 20s sloW! sET IN THE MRCIMAGE object
if (fixHeader)
  system(sprintf('alterheader -mmm  %s > /dev/null',newFilename))
end
