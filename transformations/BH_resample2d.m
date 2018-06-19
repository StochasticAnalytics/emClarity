
function [ TRANS_IMAGE ] = BH_resample2d( IMAGE, ANGLES, SHIFTS, ...
                                          CONVENTION, METHOD, DIRECTION, ...
                                          MAG, SIZEOUT)
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

% % % % Allow for anisotropic stretching 
% % % if numel(MAG) == 3
% % %   if ~strcmpi(DIRECTION, 'forward')
% % %     error('Just implementing anisotropic stretch for forward %s.\n','xform');
% % %   else
% % %     aStretch = 1/MAG(3);
% % %   end
% % % else
% % %   aStretch = 1;
% % % end

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

% In case the image is being expanded, pad accordingly
[ padVal ] = BH_multi_padVal( size(stackIN), SIZEOUT );
padLow = padVal(1,:);
padHigh= padVal(2,:);

padLow = padLow .* (padLow > 0);
padHigh = padHigh.* (padHigh > 0);

if ( useGPU)
  stackIN = BH_padZeros3d(stackIN,padLow,padHigh,'GPU','single');
else
  stackIN = BH_padZeros3d(stackIN,padLow,padHigh,'cpu','single');
end    

[ Xnew,Ynew,~,x1,y1,~ ] = BH_multi_gridCoordinates( size(stackIN), ...
                                                   'Cartesian', ...
                                                    METHOD,transformation,...
                                                    0, shiftOrigin, 0 );


     
[ padVal ] = BH_multi_padVal( SIZEOUT, size(stackIN) );
cutLow = padVal(1,:);
cutHigh= padVal(2,:);

cutLow = cutLow .* (cutLow > 0);
cutHigh = cutHigh.* (cutHigh > 0);
Xnew = Xnew(cutLow(1) + 1:end - cutHigh(1), ...
            cutLow(2) + 1:end - cutHigh(2));
          
Ynew = Ynew(cutLow(1) + 1:end - cutHigh(1), ...
            cutLow(2) + 1:end - cutHigh(2));
          
if (flgSeq)
  Xnew = Xnew{1};
  Ynew = Ynew{1};
end 

% Interpolate and write out the image.
if ( useGPU )
  TRANS_IMAGE = interpn(x1,y1,stackIN,Xnew.*1/MAG,Ynew.*1/MAG,'linear',0);
else
  TRANS_IMAGE = interpn(x1,y1,stackIN,Xnew.*1/MAG,Ynew.*1/MAG,'spline',NaN);
  TRANS_IMAGE(isnan(TRANS_IMAGE)) = mean(TRANS_IMAGE(~isnan(TRANS_IMAGE)));
end
clear Xnew Ynew Znew x1 y1 z1 stackIN ANGLES SHIFTS CONVENTION METHOD DIRECTION


end % this is the end of resample3D




