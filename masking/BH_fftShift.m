function [shiftMask] = BH_fftShift(windowRadius,imgSize,useGPU,varargin)
%Calculate the linear indices to get the shifted center of a mask (for CCC)
%   Test in 2D, extend later - for 4k aout 40x faster than fftshift then
%   center mask. About 2x as fast for full size, but the current setup only
%   would work for odd sizes.
%
%   Pass a zero to windowRadius to return full mask - should go through
%   some test cases to confirm working for sub mask, as well as 3d case.


if any(imgSize < 0)
  % make an ifft mask instead;
  flgIfft = 1;
  imgSize = abs(imgSize);
else
  flgIfft = 0;
end

doHalfGrid = 0;
if nargin == 4
  if strcmp(varargin{1},'halfgrid')
    doHalfGrid = 1;
  end
end


if numel(imgSize) == 2
  imgSize = [imgSize];
  flg2D = 1;
  shiftVect = [2,2];
else
  flg2D = 0;
  shiftVect = [2,2,2];
end



if (useGPU)
  METHOD = 'GPU';
else
  METHOD = 'cpu';
end

% if (useGPU)
%   shiftMask = zeros(windowRadius.*shiftVect+1,'single','gpuArray');
% else
%   shiftMask = zeros(windowRadius.*shiftVect+1,'single');
% end


% At somepoint, it would make sense to calculate directly rather than
% relying on ifftshift, but this is the safest "quick" approach.
if (flg2D)
  % the negative value instructs origin at corner (1:nx) rather than
  % -nx/2:nx/2 and then -1 is unshifted, while -2 shifts those coordinates
  % to a "centered" frame.
  if (doHalfGrid)
    [X,Y] = BH_multi_gridCoordinates([imgSize,1],'Cartesian',METHOD,{'none'},0,-1,0,{'halfgrid'});
    
    % TODO fix the -2 option in gridCoordinates to do this
    sY = size(Y,2);
    oY = floor(sY/2) + 1;
    isOddY = mod(sY,2);
    flgIfft = flgIfft * isOddY;
    
    tmpY = Y(:,oY+isOddY-flgIfft:end);
    Y(:,oY+flgIfft:end) = Y(:,1:oY-1+isOddY-flgIfft);
    Y(:,1:oY-1+flgIfft) = tmpY; clear tmpY
    
  else
    [X,Y] = BH_multi_gridCoordinates([imgSize,1],'Cartesian',METHOD,{'none'},0,-2,0);
  end
  
  
else
  
  if (doHalfGrid)
    imgSize
    [X,Y,Z] = BH_multi_gridCoordinates([imgSize],'Cartesian',METHOD,{'none'},0,-1,0,{'halfgrid'});
    
    % TODO fix the -2 option in gridCoordinates to do this
    sZ = size(Z,3);
    oZ = floor(sZ/2) +1;
    isOddZ = mod(sZ,2);
    flgIfftZ = flgIfft * isOddZ;
    
    sY = size(Y,2);
    oY = floor(sY/2) + 1;
    isOddY = mod(sY,2);
    flgIfftY = flgIfft * isOddY;
    
    tmpY2 = Y(:,oY+isOddY-flgIfftY:end,oZ+isOddZ-flgIfftZ:end);
    tmpY1 = Y(:,oY+isOddY-flgIfftY:end,1:oZ-1+isOddZ-flgIfftZ);
    % Put in 2 from 0
    Y(:,oY:end,oZ:end) = Y(:,1:oY-1+isOddY,1:oZ-1+isOddZ);
    % Put in 1 from 3
    Y(:,oY:end,1:oZ-1) = Y(:,1:oY-1+isOddY,oZ+isOddZ:end);
    % Put in 0 from tmpY2
    Y(:,1:oY-1+flgIfftY,1:oZ-1+flgIfftZ) = tmpY2;
    % Put in 3 from 1
    Y(:,1:oY-1+flgIfftY,oZ+flgIfftZ:end) = tmpY1; clear tmpY2 tmpY1
    
    tmpZ2 = Z(:,oY+isOddY-flgIfftY:end,oZ+isOddZ-flgIfftZ:end);
    tmpZ1 = Z(:,oY+isOddY-flgIfftY:end,1:oZ-1+isOddZ-flgIfftZ);
    % Put in 2 from 0
    Z(:,oY:end,oZ:end) = Z(:,1:oY-1+isOddY,1:oZ-1+isOddZ);
    % Put in 1 from 3
    Z(:,oY:end,1:oZ-1) = Z(:,1:oY-1+isOddY,oZ+isOddZ:end);
    % Put in 0 from tmpY2
    Z(:,1:oY-1+flgIfftY,1:oZ-1+flgIfftZ) = tmpZ2;
    % Put in 3 from 1
    Z(:,1:oY-1+flgIfftY,oZ+flgIfftZ:end) = tmpZ1; clear tmpY2 tmpY1
    
  else
    [X,Y,Z] = BH_multi_gridCoordinates(imgSize,'Cartesian',METHOD,{'none'},0,-2,0);
  end
  
end


% Cut out the grid coordinates.
if ( gather(windowRadius(1)) )
  % Set the extraval to zero, will it break?
  if windowRadius(1) < 0
    % Tread as a diameter
    padVal = BH_multi_padVal(imgSize, abs(windowRadius));
  else
    padVal = BH_multi_padVal(imgSize, windowRadius.*shiftVect+0);
  end
  
  
  X = BH_padZeros3d(X,padVal(1,:),padVal(2,:),METHOD,'single');
  Y = BH_padZeros3d(Y,padVal(1,:),padVal(2,:),METHOD,'single');
  
end

if (flg2D)
  shiftMask = sub2ind(size(X),X,Y);
else
  if ( windowRadius )
    Z = BH_padZeros3d(Z,padVal(1,:),padVal(2,:),METHOD,'single');
  end
  shiftMask = sub2ind(size(X),X,Y,Z);
end

% sub2ind forces double - copy and modify code when bored. For now, waste a
% little mem temporarily then convert if possible.

if numel(shiftMask) < 2^16 - 1
  
  shiftMask = uint16(shiftMask);
  
elseif numel(shiftMask) < 2^32 - 1
  
  shiftMask = uint32(shiftMask);
  
else
  
  shiftMask = uint64(shiftMask);
  
end

% Although they are generally avoided where this mask would be used, since
% volume sizes are curretnly optimized for FFT (and are then even in dim)
% explicilty create an ifftshift mask that handles odd volumes if
% requested.
if (flgIfft && ~doHalfGrid)
  shiftMask = ifftshift(ifftshift(shiftMask));
end

% It would make sense to just calculate this straight away.
if (doHalfGrid)
  
end

end

