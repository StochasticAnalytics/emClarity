function [ TRANS_IMAGE ] = BH_resample2d( IMAGE, ANGLES, SHIFTS, ...
                                          CONVENTION, METHOD, DIRECTION, ...
                                          MAG, SIZEOUT, varargin)
%Transform an image in 3d.
%
%   
%   Input Variables:
%   
%   IMAGE = 2d volume, or a string specifing a volume to read in.
%
%   ANGLES = Euler angles defining the desired transformation
%            Can also be a rotation matrix, predetermined, which is useful when
%            symmetrizing.
%
%   SHIFTS = Translational shifts to apply.
%
%   CONVENTION = Convention defining the transformation, e.g. 'Bah'
%       'Protomo', 'Imod', 'Spider' have yet to be updated.
%
%   METHOD = GPU, case sensitive will use gpu, otherwise cpu
%       Formerly this was an option to use interpolation other than linear,
%       however only linear interpolation is used here, and special cases where
%       more precise interpolation might be truly beneficial are dealt with
%       seperately.
%
%   DIRECTION = 'forward' standard basis to particle frame.
%               'inv'     particle frame to standard basis.
%
%
%   Output Variables:
%
%   TRANS_IMAGE = the transformed image. 
%       Independent of the input image being passed in or read in from disk, the
%       output is an image in memory.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Goals & Limitations:
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   TODO
%
%   -Test with GPU flag on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % % only cpu
if strcmp(METHOD, 'GPU')
  useGPU = true;
else
  useGPU = false;
end

doHalfGrid = 0;
extrapVal = 0;
returnComplex = false;
if nargin == 9
  if isa(varargin{1},'fourierTransformer')
    doHalfGrid = 1;
    hgSHIFTS = SHIFTS; SHIFTS = [0,0];
    hgMAG = 1/MAG; MAG = 1;
    
    if length(varargin) > 1
      if isnumeric(extrapVal)
        extrapVal = varargin{2};
      else
        returnComplex = true;
      end
    end
  elseif isnan(varargin{1})
    extrapVal = nan;
  else
    error('did not recognize the varargin');
  end
end

% Check for indication that transformations are from imod, in which case the
% origin of even images needs to be shifted.
if strcmpi(CONVENTION, 'IMOD')
  % Override in gridCoordinates b/c of switch to strictly odd size. 20171201
  shiftOrigin = [1,-0.5,-0.5,-0.5];
  CONVENTION = 'Bah';
else
  shiftOrigin = 1;
end

if ischar(IMAGE)
    % Read in the image
    stackIN = getVolume(MRCImage(IMAGE),[],[],[]);
else
    stackIN = IMAGE; clear IMAGE
end



% For individual resampling, pushing to the gpu is (with the current hardware)
% slow enough that it negates the benefit. If however the supplied image is
% already on the GPU, then this is much faster. Careful not to pull from the GPU
% and then push back.
if (useGPU)
  stackIN  = gpuArray(stackIN);
else
  stackIN = gather(stackIN); % Just in case it is passed as gpuArray
  
end


flgSeq = 0;
% Get transformed coords, if not already a matrix
if numel(ANGLES) == 3
  R = BH_defineMatrix(ANGLES, CONVENTION, DIRECTION);
  transformation = {'single',R,SHIFTS',DIRECTION,1,MAG(1)};
elseif numel(ANGLES) == 4
  % Apply an imod tansformation which includes any mag/stretch in the matrix
  % which is stored as it is output in imod a11 a12 a21 a22 ie column order
  % first reorder, and then take the inverse (not necessarily just the transpose
  % if non-rotational transformations included.
  rImod = [ANGLES(1:2),0; ANGLES(3:4),0; 0,0,1];
  R = inv(rImod);
  transformation = {'single',R,SHIFTS',DIRECTION,1,1};
elseif numel(ANGLES) == 6
  flgSeq = 1;
  % Sequential 2d alignment, extend later if needed (for mapBack)
  R = BH_defineMatrix(ANGLES(1,:), CONVENTION, DIRECTION);
  R2 = BH_defineMatrix(ANGLES(2,:), CONVENTION, DIRECTION);
  
  % The transformation from IMOD is a rotation matrix scaled by the mag, but the
  % shift values are also already scaled.
  transformation = {'sequential',R,SHIFTS(1,:)',DIRECTION,1,MAG(1) ; ...
                    'sequential',R2,SHIFTS(2,:)',DIRECTION,1,MAG(2)};
 
else
  error('ANGLES must be either three eulers or 9 rot matrix')
end

% In case the image is being expanded, pad accordingly (only for real space

if (returnComplex && (size(stackIN) ~= SIZEOUT))
  error('To return complex the input and output sizes should be the same.');
end

[ trimVal ] = BH_multi_padVal( size(stackIN), SIZEOUT );
padLow = trimVal(1,:);
padHigh= trimVal(2,:);

padLow = padLow .* (padLow > 0);
padHigh = padHigh.* (padHigh > 0);

if ~doHalfGrid
  if ( useGPU )
    stackIN = BH_padZeros3d(stackIN,padLow,padHigh,'GPU','single');
  else
    stackIN = BH_padZeros3d(stackIN,padLow,padHigh,'cpu','single');
  end 
end

if doHalfGrid
  % I don't think it matters if the interpolant is scaled to 0.5 or the
  % grid indices. The latter makes getting the values along the origin
  % easier.
  [ Xnew,Ynew,~,x1,y1,~ ] = BH_multi_gridCoordinates( varargin{1}.inputSize, ...
                                                   'Cartesian', ...
                                                    METHOD,transformation,...
                                                    0, 1, 0,{'halfgrid'});



else
  [ Xnew,Ynew,~,x1,y1,~ ] = BH_multi_gridCoordinates( size(stackIN), ...
                                                   'Cartesian', ...
                                                    METHOD,transformation,...
                                                    0, shiftOrigin, 0 );
end


if ~doHalfGrid 
  [ padVal ] = BH_multi_padVal( SIZEOUT, size(stackIN) );
  cutLow = padVal(1,:);
  cutHigh= padVal(2,:);

  cutLow = cutLow .* (cutLow > 0);
  cutHigh = cutHigh.* (cutHigh > 0);
  Xnew = Xnew(cutLow(1) + 1:end - cutHigh(1), ...
              cutLow(2) + 1:end - cutHigh(2));

  Ynew = Ynew(cutLow(1) + 1:end - cutHigh(1), ...
              cutLow(2) + 1:end - cutHigh(2));
end

if (flgSeq)
  Xnew = Xnew{1};
  Ynew = Ynew{1};
end 


% Interpolate and write out the image.
if ( useGPU )
  if (doHalfGrid)
    
    stackIN = varargin{1}.fwdSwap(varargin{1}.fwdFFT(stackIN));
    
    Xnew = Xnew ./ hgMAG;
    Ynew = Ynew ./ hgMAG;
 
    % Values from X < 0 can simply be conj(X > 0)
    hermitianMates = Xnew < 0;   
    % Values coming from X = 0 will not be correct if simply flipped.
    x_border_mask = Xnew < 1 & Xnew > -1;
    % Create an interpolant that also has +/- 2
    [Xborder,Yborder] = ndgrid(gpuArray(single(-2:2)),y1);
    
    % cut out and shift, giving the correct values for x = 0:2
    values_on_origin = stackIN(1:5,:);
    values_on_origin =  circshift(values_on_origin,-3);
    % If it is even sized, the first column is not handled correctly, but
    % is almost always going to be zero anyway b/c its at Nyquist
    values_on_origin(Xborder == -1 & Yborder > 0) = conj(values_on_origin(Xborder == 1 & Yborder < 0 & Yborder >= -y1(end)));
    values_on_origin(Xborder == -2 & Yborder > 0) = conj(values_on_origin(Xborder == 2 & Yborder < 0 & Yborder >= -y1(end)));
    values_on_origin(Xborder == -1 & Yborder < 0 & Yborder >= -y1(end)) = conj(values_on_origin(Xborder == 1 & Yborder > 0 & Yborder <= y1(end)));
    values_on_origin(Xborder == -2 & Yborder < 0 & Yborder >= -y1(end)) = conj(values_on_origin(Xborder == 2 & Yborder > 0 & Yborder <= y1(end)));
    
    % We need to cut these out prior to flipping the rest of the hermitian
    % mates
    Xnew_border = Xnew(x_border_mask);
    Ynew_border = Ynew(x_border_mask);

    % Invert the coordinates, take conjugate after interpolating
    Xnew(hermitianMates) = -1.*Xnew(hermitianMates);
    Ynew(hermitianMates) = -1.*Ynew(hermitianMates);
      
  
    TRANS_IMAGE = interpn(x1,y1,stackIN,Xnew,Ynew,'linear',0);
    TRANS_IMAGE(hermitianMates) = conj(TRANS_IMAGE(hermitianMates));
    % Now go back and replace the values that came from locations near x =
    % 0.
    TRANS_IMAGE(x_border_mask) = interpn(Xborder,Yborder,values_on_origin,Xnew_border, Ynew_border,'linear',0);

    clear values_on_origin Xborder Yborder x_border_mask
    

    if (hgMAG ~= 1 || hgSHIFTS(1) || hgSHIFTS(2))
      isCentered=1;
      TRANS_IMAGE = varargin{1}.shiftStretch(TRANS_IMAGE,hgSHIFTS,hgMAG,isCentered);
    end
    

    if (returnComplex)
     TRANS_IMAGE = varargin{1}.invSwap(TRANS_IMAGE);
                               
    else
     TRANS_IMAGE = BH_padZeros3d( varargin{1}.invFFT(varargin{1}.invSwap(TRANS_IMAGE),2),...
                                trimVal(1,:),trimVal(2,:),'GPU','single');      
    end
    

  else
    

    TRANS_IMAGE = interpn(x1,y1,stackIN,Xnew,Ynew,'linear',extrapVal);

    
  end
else
  TRANS_IMAGE = interpn(x1,y1,stackIN,Xnew,Ynew,'spline',NaN);
  TRANS_IMAGE(isnan(TRANS_IMAGE)) = extrapVal;
end
clear Xnew Ynew Znew x1 y1 z1 stackIN ANGLES SHIFTS CONVENTION METHOD DIRECTION


end % this is the end of resample3D




