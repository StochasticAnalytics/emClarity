%putImage       Replace the specified image of a MRCImage object
%
%   mRCImage = putImage(mRCImage, img, index)
%
%   img         The specified image as a 2D array
%
%   mRCImage    The MRCImage object containing the stack of interest
%
%   index       The index of the image to replace, indices start at 1.
%
%
%   Replaces the specifed image in the MRCImage object.
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

function mRCImage = putImage(mRCImage, img, index)

% Error checking
if index < 1
  PEETError('Image index out of range!');
end

index = double(index); % avoid overflow in idxDataStart computation

[nX, nY] = size(img);
if (nX ~= mRCImage.header.nX || nY ~= mRCImage.header.nY)
  fprintf('Image size: %d x %d\n', nX, nY);
  fprintf('Stack image size: %d x %d\n', ...
          mRCImage.header.nX,  mRCImage.header.nY);
  PEETError('Image is not the same size as the stack!');
end

if mRCImage.flgVolume
  mRCImage.volume(:,:, index) = img;
else
  % Check the permissions on the file ID, reopen it if needed
  mRCImage.fid = openWritable(mRCImage);
    
  modeStr = getModeString(mRCImage);
  if strcmp(modeStr, 'int16*2') || strcmp(modeStr, 'float32*2')
    modeStr = modeStr(1 : end - 2);
    flgComplex = true;
  else
    flgComplex = false;
  end
  nImageElements = mRCImage.header.nX * mRCImage.header.nY;
  % Move the file pointer to the beginning requested image
  nImageBytes = nImageElements * getModeBytes(mRCImage);
  idxDataStart = mRCImage.dataIndex + (index - 1) * nImageBytes;
  moveOrExtendFP(mRCImage.fid, idxDataStart, idxDataStart + nImageBytes);
  
  % Write out the image from the MRC file
  if flgComplex      % handle complex image
    if modeStr(1) == 'i'
      temp = int16(zeros(2 * nX, nY));
    else
      temp = single(zeros(2 * nX, nY));
    end
    temp(1:2:2*nX-1, :) = real(img);
    temp(2:2:2*nX, :) = imag(img);
    count = fwrite(mRCImage.fid, double(temp), modeStr);
    if count ~= 2 * nImageElements
      PEETError(['Expected ' int2str(nImageElements) ' elements, '     ...
        'wrote ' int2str(count / 2) '\n  ' ferror(mRCImage.fid)]);
    end
  else               % normal (not complex) image
    count = fwrite(mRCImage.fid, double(img), modeStr);
    if count ~= nImageElements
      PEETError(['Expected ' int2str(nImageElements) ' elements, '     ...
        'wrote ' int2str(count) '\n  ' ferror(mRCImage.fid)]);
    end
  end
  
  % Update the header if slice is greater than the current number of
  % slice in the header
  if index > mRCImage.header.nZ
    mRCImage.header.nZ = index;
    mRCImage = writeHeader(mRCImage);
  end
end


function moveOrExtendFP(fid, index, totalSize)

status = fseek(fid, index, 'bof');

if status
  [message, errno] = ferror(fid); %#ok<ASGLU>
  if errno == -27
    %  Seek to the end of the file
    stat2 = fseek(fid, 0, 'eof');
    if stat2
      PEETError(['Could not seek to the end of the existing file\n  '  ...
        ferror(fid)]);
    end
    
    % Write out enough zeros to extend the file
    pos = ftell(fid);
    nZeros = totalSize - pos;
    count = fwrite(fid, zeros(nZeros, 1), 'uchar');
    if count ~= nZeros
      PEETError(['Unable to extend file!\n  ' ferror(fid)]);
    end
    
    % Attempt to move the file pointer again
    stat3 = fseek(fid, index, 'bof');
    if stat3
      PEETError(['Could not move the file pointer after extending file' ...
        '\n  ' ferror(fid)]);
    end
  else
    PEETError(ferror(fid));
  end
end
