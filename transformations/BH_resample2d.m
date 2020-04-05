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
if nargin == 9
  if isa(varargin{1},'fourierTransformer')
    doHalfGrid = 1;
    hgSHIFTS = SHIFTS; SHIFTS = [0,0];
    hgMAG = 1/MAG; MAG = 1;
    
    if length(varargin) > 1
      extrapVal = varargin{2};
      if ~isnumeric(extrapVal)
        error('Extrapolation value passed to resample2d is non numeric!')
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
  [ Xnew,Ynew,~,x1,y1,~ ] = BH_multi_gridCoordinates( varargin{1}.inputSize, ...
                                                   'Cartesian', ...
                                                    METHOD,transformation,...
                                                    1, 1, 0,{'halfgrid'});



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
    hermitianMates = Xnew < 0;
    Xnew(hermitianMates) = -1.*Xnew(hermitianMates);
    Ynew(hermitianMates) = -1.*Ynew(hermitianMates);


    stackIN = varargin{1}.fwdSwap(varargin{1}.fwdFFT(stackIN));
    TRANS_IMAGE = interpn(x1,y1,stackIN,Xnew.*(1/hgMAG),Ynew.*(1/hgMAG),'linear',0);
    
    TRANS_IMAGE(hermitianMates) = conj(TRANS_IMAGE(hermitianMates));

    if (hgMAG ~= 1 || hgSHIFTS(1) || hgSHIFTS(2))
      isCentered=1;
      TRANS_IMAGE = varargin{1}.shiftStretch(TRANS_IMAGE,hgSHIFTS,hgMAG,isCentered);
    end
    

     TRANS_IMAGE = BH_padZeros3d( varargin{1}.invFFT(varargin{1}.invSwap(TRANS_IMAGE),2),...
                                trimVal(1,:),trimVal(2,:),'GPU','single');

  else
    

    TRANS_IMAGE = interpn(x1,y1,stackIN,Xnew,Ynew,'linear',extrapVal);

    
  end
else
  TRANS_IMAGE = interpn(x1,y1,stackIN,Xnew,Ynew,'spline',NaN);
  TRANS_IMAGE(isnan(TRANS_IMAGE)) = extrapVal;
end
clear Xnew Ynew Znew x1 y1 z1 stackIN ANGLES SHIFTS CONVENTION METHOD DIRECTION


end % this is the end of resample3D




