function [ MASK ] = BH_mask3d( SHAPE, SIZE, RADIUS, CENTER )
%Create a mask for real space 3d images.
%
%   Input variables:
%
%   SHAPE = 'sphere' or 'cylinder' : string
%   SIZE  = [x, y, z] dimension of mask : vector, int
%   RADIUS= [rx, ry, rz] radius in each dimension : vector, float
%
%   CENTER= [cx, cy, cz] center of the mask : vector, float
%
%   Output variables:
%
%   MASK  = 3d MRC image file, single precision float.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Goals & Restrictions:
%
%   Many image processing packages allow for the user to define the
%   apodization in real space, or similarly the fall off in reciprocal
%   space. Any sharp transitions in either can influence alignment
%   negatively, and even result in falsely inflating reslution estimations.
%
%   This is avoided by --
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   TODO
%   - complete the Goals & Restrictions outline.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


pixelFallOff = 6;

% Make sure that the begining and end of the taper happens at a predictable
% place after convolution.
convCutLow  = 0.025;

%gaussian = @(x,m,APOSIZE) (exp( -1.*(x-m).^2 ./ (2.*APOSIZE.^2) ));
taper = 0.5+0.5.*cos((((1:pixelFallOff)).*pi)./(length((1:pixelFallOff+1))));

% Check that input is approprate.
[mShape, mSize, mRadius, mCenter, binaryMask] = ...
  parseVariables( SHAPE, SIZE, RADIUS, CENTER);





METHOD = 'cpu';

clear mWindow
if strcmpi(mShape, 'rectangle')
  mWindow(2*mRadius(1)+14-1, 2*mRadius(2)+14-1, 2*mRadius(3)+14-1) = (single(0));
  mWindow = mWindow + 1;
else
  mWindow(mSize(1), mSize(2), mSize(3)) = (single(0));
  mWindow = mWindow + 1;
end

[ gaussKernel ] = (BH_multi_gaussian3d(5, 0.5 ));



if strcmpi(mShape, 'sphere')
  
  [ G1,G2,G3,~,~,~ ] = BH_multi_gridCoordinates( mSize, 'Cartesian', METHOD, ...
    {'single',...
    [1,0,0;0,1,0;0,0,1],...
    mCenter','forward',1,1}, ...
    0, 1, 0 );
  
  ellipsoid = (G1./mRadius(1)).^2 + (G2./mRadius(2)).^2 + (G3./mRadius(3)).^2 ;
  fullMask = (ellipsoid <= 1) ;
  mWindow = mWindow .* fullMask;
  for iShell = 1:pixelFallOff
    ellipsoid = (G1./(mRadius(1)+iShell)).^2 + ...
      (G2./(mRadius(2)+iShell)).^2 + ...
      (G3./(mRadius(3)+iShell)).^2 ;
    ellipsoid = (ellipsoid <= 1) - mWindow;
    mWindow = mWindow + (ellipsoid .* taper(iShell));
  end
  clear borderMask  G1 G2 G3
  mWindow = convn(mWindow, gaussKernel, 'same') ;
  mWindow = mWindow ./ max(mWindow(:));
  % Set the edges at mRadius back to 1 after convolution
  mWindow(fullMask) = 1;
  mWindow(mWindow < convCutLow)  = 0;
end

% z filtering for cylinder
if strcmpi(mShape, 'cylinder')
  
  [ G1,G2,G3,~,~,~ ] = BH_multi_gridCoordinates( mSize, 'Cartesian', METHOD, ...
    {'single',...
    [1,0,0;0,1,0;0,0,1],...
    mCenter','forward',1,1}, ...
    0, 1, 0 );
  G3 = abs(G3);
  ellipsoid = (G1./mRadius(1)).^2 + (G2./mRadius(2)).^2;
  fullMask = (ellipsoid <= 1) & (G3 <= mRadius(3));
  mWindow = mWindow .* fullMask;
  for iShell = 1:pixelFallOff
    ellipsoid = (G1./(mRadius(1)+iShell)).^2 + ...
      (G2./(mRadius(2)+iShell)).^2;
    
    ellipsoid = ((ellipsoid <= 1).*(G3 <= mRadius(3)+iShell)) - mWindow;
    
    mWindow = mWindow + (ellipsoid .* taper(iShell));
    
    mWindow = convn(mWindow, gaussKernel, 'same') ;
    mWindow = mWindow ./ max(mWindow(:));
    % Set the edges at mRadius back to 1 after convolution
    mWindow(fullMask) = 1;
    mWindow(mWindow < convCutLow)  = 0;
  end
  clear borderMask  G1 G2 G3
  
end

% rectangular filtering

if strcmpi(mShape, 'rectangle')
  
  
  [ padVal ] = BH_multi_padVal(size(mWindow), mSize);
  prePad = padVal(1,:);
  postPad= padVal(2,:);
  % Still doesn't really address an off-centered mask, but for this particular
  % case, this works with the new definition of the origin.
  
  shiftOrigin = ~mod(mSize,2);
  
  prePad = prePad + shiftOrigin ;
  postPad = postPad - shiftOrigin;
  
  % First trim the window down in case the area is small.
  preTrim = abs(prePad) .* (prePad < 0);
  postTrim = abs(postPad) .* (postPad < 0);
  prePad = prePad .* (prePad >= 0);
  postPad = postPad .* (postPad >= 0);
  
  
  
  
  mWindow = mWindow(preTrim(1) + 1: end - postTrim(1), ...
    preTrim(2) + 1: end - postTrim(2), ...
    preTrim(3) + 1: end - postTrim(3));
  
  mWindow = BH_padZeros3d(mWindow,prePad,postPad,METHOD,'singleTaper');
  
end
clear X Y Z x y z

if strcmpi(mShape, 'BINARY')
  % Select a pretty aggressive cutoff, then dilate by at 10A or at least 2
  % pixels.
  % % %   binaryMask = binaryMask - mean(binaryMask(:));
  % % %   binaryMask = binaryMask ./ rms(binaryMask(:));
  pixelSize = SIZE
  rectMask = BH_mask3d_cpu('rectangle',size(binaryMask),(size(binaryMask)./2 - 7),[0,0,0]);
  binaryVol = BH_bandLimitCenterNormalize_cpu(binaryMask.*rectMask, ...
    BH_bandpass3d(size(binaryMask),0,0,24,METHOD,pixelSize),...
    (rectMask > .01),...
    [0,0,0;0,0,0], 'double');
  binaryVol = real(ifftn(binaryVol));
  
  %  figure, imshow3D(binaryMask)
  binaryMask = (binaryVol > 3.0);%(log(kurtosis(binaryMask(:)./3))));
  
  a = single((binaryMask));
  
  dilationKernel = (BH_multi_gaussian3d([3,3,3],3.0));
  b = a;
  %  figure, imshow3D(gather(b))
  for i = 1:mRadius + 5
    b = single(~a.*convn(b,dilationKernel,'same') > 0.2)+a;
  end
  a = b;
  
  
  binaryMask = (binaryVol.*(a > 0.5) > 0.3);%(log(kurtosis(binaryMask(:)./3))));
  a = single((binaryMask));
  b = a;
  %  figure, imshow3D(gather(b))
  for i = 1:mRadius
    b = single(~a.*convn(b,dilationKernel,'same') > 0.2)+a;
  end
  a = b;
  
  
  
  %  figure, imshow3D(gather(a))
  taperKernel = (BH_multi_gaussian3d([3,3,3],1.5));
  for i = 1:5
    b = ~(a).*convn(b,taperKernel,'same')+a;
  end
  %  figure, imshow3D(gather(b))
  smoothKernel = (BH_multi_gaussian3d([5,5,5],0.65));
  for i = 1:2
    b = convn(b,smoothKernel,'same')+b;
  end
  % figure, imshow3D(gather(b))
  mWindow = gather(b ./ max(b(:)));
  % figure, imshow3D(gather(mWindow))
  clear a b avgM
  
end

%mWindow(mWindow < .5.*fallOffVal) = 0;
MASK = mWindow;

clear mWindow

clearvars -except MASK
end % end of BH_mask3d function.


function [mShape, mSize, mRadius, mCenter, binaryMask] = ...
  parseVariables( SHAPE, SIZE, RADIUS, CENTER)
%Check for correct input parameters:

% If SHAPE is an image, default to rectangular size of image and apply mask,
% returning the masked image rather than the image itself.

binaryMask = '';
if isnumeric(SHAPE)
  
  binaryMask = SHAPE;
  clear SHAPE
  mShape = 'BINARY';
  mSize = size(binaryMask);
  %mRadius = mSize ./ 2 - 6;
  mCenter = [0,0,0];
  
  % Should be Ang/pixel, used to calculate dilation.
  pixelSize = SIZE;
  mRadius = max(3,floor(10./pixelSize))
else
  
  if strcmpi(SHAPE, 'sphere')
    mShape = SHAPE;
    
    dim = 1;
  elseif strcmpi(SHAPE, 'cylinder')
    mShape = SHAPE;
    
    dim = 2;
  elseif strcmpi(SHAPE, 'rectangle')
    mShape = SHAPE;
    dim = 3;
  else
    error('SHAPE must be sphere, cylinder, or rectangle not %s', SHAPE)
  end
  
  if ~isnumeric(SIZE) || ~(length(SIZE) == 3)
    error('SIZE must be a vector in R3')
  else
    mSize = SIZE;
  end
  
  if ~isnumeric(RADIUS) || ~(length(RADIUS) == 3)
    error('RADIUS must be radius radius [height], dimension %d', dim)
  else
    mRadius = RADIUS;
    % % %       % For cylindircal/spherical reduce the radius by one (instead of 1 at the edge
    % % %   % then the value will be ~.995)
    % % %     if strcmpi(SHAPE, 'sphere')  || strcmpi(SHAPE, 'cylinder')
    % % %       mRadius = mRadius - 1;
    % % %     end
  end
  
  if ~isnumeric(CENTER) || ~(length(CENTER) == 3)
    error('CENTER must be a vector in R3')
  else
    mCenter = CENTER;
    
  end
end




end % end of parseVariables function

