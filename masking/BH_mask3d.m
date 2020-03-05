function [ MASK, volCOM ] = BH_mask3d( SHAPE, SIZE, RADIUS, CENTER, varargin)
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
volCOM = [0,0,0];
volTight = 0;
flgCOM = 0;
flg3d = 1;

  global bh_global_binary_mask_low_pass
  global bh_global_binary_mask_threshold
  global bh_global_vol_est_scaling
  if isempty(bh_global_binary_mask_low_pass)
    bh_global_binary_mask_low_pass = 14;
  end
  if isempty(bh_global_binary_mask_threshold)
    bh_global_binary_mask_threshold = 2.5;
  end
  if isempty(bh_global_vol_est_scaling)
    % The low pass version of the map used for the estimate overestimates
    % the molecular volume at the hydration radius of the underlying atoms.
    % This flag will override the value I've calculated which depends on
    % the masking resolution. TODO when the map resolution is lower than
    % the masking resolution, this will again underestimate the scaling,
    % artificialy *de*pressing the FSC
    bh_global_vol_est_scaling = 0.0;
  end
  
  if (bh_global_vol_est_scaling == 0)
    localParticleScaling = (-2.8e-3) .* bh_global_binary_mask_low_pass.^2 + ...
                            0.14 .*     bh_global_binary_mask_low_pass + 1.5;
  else
    localParticleScaling = bh_global_vol_est_scaling;
  end
    
  
asymmetricRestriction = 0;
if nargin > 5
  if strcmp(varargin{1},'2d')
    flg3d = 0;
  end
  asymmetricRestriction = varargin{2};

elseif nargin > 4
  
  if strcmp(varargin{1},'2d')
    flg3d = 0;
  else
    flgCOM = 1; 
  end
  
end
pixelFallOff = 6;

% Make sure that the begining and end of the taper happens at a predictable
% place after convolution.
convCutLow  = 0.025;

%gaussian = @(x,m,APOSIZE) (exp( -1.*(x-m).^2 ./ (2.*APOSIZE.^2) ));
taper = 0.5+0.5.*cos((((1:pixelFallOff)).*pi)./(length((1:pixelFallOff+1))));


% Check that input is approprate.
[mShape, mSize, mRadius, mCenter, binaryMask,fscMask] = ...
                              parseVariables( SHAPE, SIZE, RADIUS, CENTER, flg3d);


if (asymmetricRestriction ~= 0 && ~strcmpi(mShape,'Cylinder'))
  error('Experimental symmetry restricting mask is only available for cylinders right now.');
end

  METHOD = 'GPU';

  clear mWindow
  if strcmpi(mShape, 'rectangle')
    if (flg3d)
      mWindow(2*mRadius(1)+14-1, 2*mRadius(2)+14-1, 2*mRadius(3)+14-1) = ...
                                                            gpuArray(single(0));
    else
      mWindow(2*mRadius(1)+14-1, 2*mRadius(2)+14-1) = gpuArray(single(0));
    end
    mWindow = mWindow + 1; 
  else
    mWindow = zeros(mSize,'single','gpuArray');
    mWindow = mWindow + 1;
  end
  
  if (flg3d)
    [ gaussKernel ] = gpuArray(BH_multi_gaussian3d(5, 0.5 ));
  else
    [ gaussKernel ] = gpuArray(BH_multi_gaussian2d(5, 0.5 ));
  end



if strcmpi(mShape, 'sphere')

  if (flg3d)
    [ G1,G2,G3,~,~,~ ] = BH_multi_gridCoordinates( mSize, 'Cartesian',METHOD,...
                                                  {'single',...
                                                  [1,0,0;0,1,0;0,0,1],...
                                                  mCenter','forward',1,1}, ...
                                                  0, 1, 0 );
  else
    [ G1,G2,~,~,~,~ ] = BH_multi_gridCoordinates( mSize, 'Cartesian',METHOD, ...
                                                  {'single',...
                                                  [1,0,0;0,1,0;0,0,1],...
                                                  mCenter','forward',1,1}, ...
                                                  0, 1, 0 );
  end

  
    ellipsoid = (G1./mRadius(1)).^2 + (G2./mRadius(2)).^2;
  if (flg3d)
    ellipsoid = ellipsoid + (G3./mRadius(3)).^2;
  end
  
fullMask = (ellipsoid <= 1) ;
mWindow = mWindow .* fullMask; 
for iShell = 1:pixelFallOff

  ellipsoid = (G1./(mRadius(1)+iShell)).^2 + (G2./(mRadius(2)+iShell)).^2;
  if (flg3d)
    ellipsoid = ellipsoid + (G3./(mRadius(3)+iShell)).^2 ;
  end
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
  
  if (flg3d)
    [ G1,G2,G3,~,~,~ ] = BH_multi_gridCoordinates( mSize, 'Cartesian',METHOD,...
                                                  {'single',...
                                                  [1,0,0;0,1,0;0,0,1],...
                                                  mCenter','forward',1,1}, ...
                                                  0, 1, 0 );
    G3 = abs(G3);                                                
  else
    [ G1,G2,~,~,~,~ ] = BH_multi_gridCoordinates( mSize, 'Cartesian',METHOD,...
                                                  {'single',...
                                                  [1,0,0;0,1,0;0,0,1],...
                                                  mCenter','forward',1,1}, ...
                                                  0, 1, 0 );
  end


ellipsoid = (G1./mRadius(1)).^2 + (G2./mRadius(2)).^2;

if (flg3d)
  fullMask = (ellipsoid <= 1) & (G3 <= mRadius(3));
else
  fullMask = (ellipsoid <= 1);
end

if (asymmetricRestriction)
  if (flg3d)
    [ ~,angles,~,~,~,~ ] = BH_multi_gridCoordinates( mSize, 'Cylindrical',METHOD,...
                                                  {'single',...
                                                  [1,0,0;0,1,0;0,0,1],...
                                                  mCenter','forward',1,1}, ...
                                                  0, 1, 0 );
    G3 = abs(G3);                                                
  else
    [ ~,angles,~,~,~,~ ] = BH_multi_gridCoordinates( mSize, 'Cylindrical',METHOD,...
                                                  {'single',...
                                                  [1,0,0;0,1,0;0,0,1],...
                                                  mCenter','forward',1,1}, ...
                                                  0, 1, 0 );
  end
  
  sectorMax = 2*pi/asymmetricRestriction * 1.025;
  angles = (angles > (2*pi-sectorMax/2) | angles < sectorMax/2);
  gc = BH_multi_gaussian3d(-1.*size(angles),1.5);
  mWindow = real(ifftn(fftn(angles.*fullMask).*gc));
  clear gc
  mWindow = mWindow ./ max(angles(:));
  
else
   
  mWindow = mWindow .* fullMask; 

for iShell = 1:pixelFallOff
  ellipsoid = (G1./(mRadius(1)+iShell)).^2 + ...
              (G2./(mRadius(2)+iShell)).^2;
  if (flg3d)
    ellipsoid = ((ellipsoid <= 1).*(G3 <= mRadius(3)+iShell)) - mWindow;
  else
    ellipsoid = ((ellipsoid <= 1) - mWindow);
  end

    mWindow = mWindow + (ellipsoid .* taper(iShell));
 
    mWindow = convn(mWindow, gaussKernel, 'same') ;
    mWindow = mWindow ./ max(mWindow(:));
    % Set the edges at mRadius back to 1 after convolution
    mWindow(fullMask) = 1; 
    mWindow(mWindow < convCutLow)  = 0;
end
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
    
  
  % UPDATE TO USE padZeros3d TO TRIM
  if (flg3d)
    mWindow = mWindow(preTrim(1) + 1: end - postTrim(1), ...
                      preTrim(2) + 1: end - postTrim(2), ...
                      preTrim(3) + 1: end - postTrim(3));
  else
    mWindow = mWindow(preTrim(1) + 1: end - postTrim(1), ...
                      preTrim(2) + 1: end - postTrim(2));
  end

  mWindow = BH_padZeros3d(mWindow,prePad,postPad,METHOD,'singleTaper');

end
clear X Y Z x y z

if strcmpi(mShape, 'BINARY')
  % Select a pretty aggressive cutoff, then dilate by at 10A or at least 2
  % pixels.
% % %   binaryMask = binaryMask - mean(binaryMask(:));
% % %   binaryMask = binaryMask ./ rms(binaryMask(:));
  pixelSize = SIZE;

  if (flg3d)
    rectMask = BH_mask3d('rectangle',size(binaryMask), ...
                                             (size(binaryMask)./2 - 7),[0,0,0]);
  else
    rectMask = BH_mask3d('rectangle',size(binaryMask), ...
                                        (size(binaryMask)./2 - 7),[0,0,0],'2d');
  end



  if (flg3d)
    binaryMask = BH_bandLimitCenterNormalize(medfilt3(gather(binaryMask),[3,3,3]).*rectMask, ...
                                          BH_bandpass3d(size(binaryMask), ...
                                                        0,0, ...
                                                        bh_global_binary_mask_low_pass, ...
                                                        'GPU',pixelSize),...
                                           (rectMask > .01),...
                                           [0,0,0;0,0,0], 'single');
  else
    binaryMask = BH_bandLimitCenterNormalize(medfilt2(gather(binaryMask),[3,3]).*rectMask, ...
                                          BH_bandpass3d([size(binaryMask),1], ...
                                                        0,0,...
                                                        bh_global_binary_mask_low_pass, ...
                                                        'GPU',pixelSize),...
                                           (rectMask > .01),...
                                           [0,0;0,0], 'single');
  end
                                        
  binaryMask = real(ifftn(binaryMask)).*rectMask;
% % %   
% % %   binaryOutside = (binaryVol > 0.5);
% % %   if (fscMask)
% % %     maxThreshold = mean(binaryVol(binaryOutside)) + ...
% % %                  1.5.* std(binaryVol(binaryOutside))
% % %   else
% % %     
% % %     maxThreshold = mean(binaryVol(binaryOutside)) + ...
% % %                  3.* std(binaryVol(binaryOutside))
% % %   end
% % % %   maxThreshold = (kurtosis(binaryVol(:))-3).^0.5;
% % %   binaryMask = (binaryVol > maxThreshold );

%   % Just to test a hunch
%   if (flg3d)
%     binarySmooth = (medfilt3(gather(binaryMask),[3,3,3]));
%   else
%     binarySmooth = (medfilt2(gather(binaryMask),[3,3]));
%   end
  

  maxThreshold = bh_global_binary_mask_threshold.*(std(binaryMask(binaryMask(:)>0)));
  binaryVol = binaryMask > maxThreshold;
  clear binarySmooth

  if (fscMask)
    dilationThresholds = [ 0.9 0.85  0.75 0.7 0.65 0.5 0.35 0.2 0.1   ] ;
  else
    dilationThresholds =  [ 1.0000 0.9   ] ;
  end
  for threshold =  dilationThresholds.*maxThreshold
    threshold 
    maxThreshold
  %  figure, imshow3D(binaryMask) 0.6923 0.6154 
   if threshold >= 0 
    currentMask = single(gpuArray(binaryVol));

    if (flg3d)
      dilationKernel = gpuArray(BH_multi_gaussian3d(3.*[1,1,1],3.0));
    else
      dilationKernel = gpuArray(BH_multi_gaussian2d(3.*[1,1],3.0));
    end
    size(dilationKernel)
    %figure, imshow3D(gather(b))
  
    if (fscMask)
      dilationIter = ceil(threshold.^2./3);
    else
      dilationIter = ceil(threshold.^2./3);
    end
    
    % Grow 
    for i = 1:dilationIter
      currentMask = single(~currentMask.*convn(currentMask,dilationKernel,'same') > 0.00)+currentMask;   
    end   
    

      binaryVol = (binaryMask.*currentMask > threshold);

%     figure, imshow3D(gather(binaryVol))
    

   end
  end
  % It should already be binary here right?
 
  particleVolEstimate = sum(binaryVol(:) == 1);
  % Return volume estimate for visualization/ analysis. Zero out first pixel in FSC.
  volCOM = single(binaryVol);
  
  for i = 1:size(binaryVol,3)
    bw = bwdist(binaryVol(:,:,i));
    binaryVol(:,:,i) = binaryVol(:,:,i) + bw < mRadius / 2 ; %%% Comment out for tight mask demo
  end
   clear bw
%   a = (binaryMask.*binaryVol > threshold);
  % Bringback the edges
  %particleVolEstimate = sum(binaryVol(:));

  if (fscMask)
  %  figure, imshow3D(gather(a)) 
    currentMask = single(gpuArray(binaryVol));
   
    if (flg3d)
      taperKernel = gpuArray(BH_multi_gaussian3d(4.*[1,1,1],1.75));
    else
      taperKernel = gpuArray(BH_multi_gaussian2d(4.*[1,1],1.75));
    end
    
    volCOM = convn(volCOM,taperKernel,'same');
    volCOM = convn(sqrt(volCOM),taperKernel,'same'); clear taperKernel
    
    
    if (flg3d)
      smoothKernel = gpuArray(BH_multi_gaussian3d(mRadius.*[1,1,1],mRadius./2));
    else
      smoothKernel = gpuArray(BH_multi_gaussian2d(mRadius.*[1,1],mRadius./2));      
    end
    for i = 1:2 %%% set to one for tight mask demo
      currentMask = convn(currentMask,smoothKernel,'same');
      %       figure, imshow3D(gather(currentMask))  
    end
    currentMask = currentMask ./ max(currentMask(:));
   
    % We assume the particle envelope only cuts through solvent, which
    % is not always the case in cryoSTAC of extended assemblies. To estimate
    % The power reduction in the signal is (I think) a good idea, because
    % Otherwise the noise reduction term is too strong. 
    % Get the reduction in power in the taper region due to masking
    powerReduction = sum(abs(binaryMask(:)).^2.*(currentMask(:)>0))./ ...
                     sum(abs(binaryMask(:)).^2.* currentMask(:));
    
    maskVolume = sum(currentMask(:)>0);
    particleVolEstimate = particleVolEstimate ./ localParticleScaling;
    particleFraction = particleVolEstimate ./ maskVolume .* powerReduction;
    fprintf('Estimated partVol, %d voxels\nmaskVol %d voxels\npwrReduction %2.3f\npartFract %2.3f\n',...
             particleVolEstimate, maskVolume, powerReduction,particleFraction);
           
    % Should probably use varargout, but for now, returning the center of
    % mass is not done at the same stage as a particle volume estimate.

    volCOM(1) = 1/particleFraction;
   
  % figure, imshow3D(gather(b))  
    mWindow = gather(currentMask);
  else
    mWindow = gather(binaryVol);
  end
    %figure, imshow3D(gather(mWindow))
  clear a b avgM  
    
end

%mWindow(mWindow < .5.*fallOffVal) = 0;
if (flg3d)
  MASK = BH_padZeros3d(mWindow, [0,0,0],[0,0,0],'GPU','singleTaper');
else
  MASK = BH_padZeros3d(mWindow, [0,0],[0,0],'GPU','singleTaper');
end
clear mWindow
if (flgCOM)
  % get the center of mass of the low-pass filtered volume and return.
  [X,Y,Z,~,~,~] = BH_multi_gridCoordinates(size(binaryVol),'Cartesian',...
                                                     'GPU',{'none'},0,1,0);
                                                   
  binaryVol = (binaryVol - min(binaryVol(MASK > 0.01)) ).*MASK; 
  
  volCOM = [sum(sum(sum(binaryVol.*X))), ...
            sum(sum(sum(binaryVol.*Y))), ...
            sum(sum(sum(binaryVol.*Z)))] ./ sum(binaryVol(:));
  clear binaryVol X Y Z binaryMask

end  

  clearvars -except MASK  volCOM                     
end % end of BH_mask3d function.


function [mShape, mSize, mRadius, mCenter, binaryMask,fscMask] = ...
                              parseVariables( SHAPE, SIZE, RADIUS, CENTER, flg3d)
%Check for correct input parameters:

% If SHAPE is an image, default to rectangular size of image and apply mask,
% returning the masked image rather than the image itself.

binaryMask = '';
fscMask = '';

if isnumeric(SHAPE)
  binaryMask = gpuArray(SHAPE);
  clear SHAPE
  mShape = 'BINARY';
  mSize = size(binaryMask);
  %mRadius = mSize ./ 2 - 6;
  if flg3d
    mCenter = [0,0,0];
  else
    mCenter = [0,0];
  end
  
  % Should be Ang/pixel, used to calculate dilation.
  pixelSize = abs(SIZE);
  if SIZE > 0
    % Use a more conservative dilation for FSC masking
    fscMask = 1;
  else
    % Use a more stringent dilation (mask is just as soft) to focus alignment on
    % stronger density
    fscMask = 0;
  end
  
  mRadius = max(3,floor(10./pixelSize));
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

  if ~isnumeric(SIZE) || ~((length(SIZE) == 3) || (length(SIZE) == 2))
  
    error('SIZE must be a vector in R3')
  else
    mSize = SIZE;
  end

  if ~isnumeric(RADIUS) || ~((length(RADIUS) == 3) || (length(SIZE) == 2))
    error('RADIUS must be radius radius [height], dimension %d', dim)
  else
    mRadius = RADIUS;

  end

  if ~isnumeric(CENTER) || ~((length(CENTER) == 3) || (length(SIZE) == 2))
      error('CENTER must be a vector in R3')
  else
      mCenter = CENTER;

  end
end




end % end of parseVariables function

