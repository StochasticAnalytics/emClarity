function [ TRANS_IMAGE, x1, y1, z1 ] = BH_resample3d( IMAGE, ANGLES, SHIFTS, ...
                                          CONVENTION, METHOD, DIRECTION, ...
                                          varargin)
%Transform an image in 3d.
%
%   
%  Input Variables:
%   
%   IMAGE = 3d volume, or a string specifing a volume to read in.
%
%   ANGLES = Single Volume: Euler angles defining the desired transformation
%            Can also be a rotation matrix, predetermined, which is useful when
%            symmetrizing.
%
%            Multi Volumes: Will return a stack of resampled images if ANGLES is
%            a cell. In this case, the cell should be 
%            {nImages,2} = {Angles, shifts}
%
%            then the SHIFTS variable should be a cell{1,1} with the 2x3 matrix
%            specifying the values to trim each rotated volume to, which defines
%            windowSize --> maskSize in most of my programs.
%
%   SHIFTS = Single volume: Translational shifts to apply.
%            Multi volumes: Trim value, windowSize --> maskSize
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

volBinary = 'noVol';
flgMask = 0;

if strcmp(METHOD, 'GPU')
  useGPU = 1;
else
  useGPU = 0;
end

inputVectors = 0;
if nargin == 7
  inputVectors = 1;
end



if ischar(IMAGE)
    % Read in the image
    IMAGE = getVolume(MRCImage(IMAGE));
end

flgComplex = 0;
flgComplexShift = 0;
if ~isreal(IMAGE)
  flgComplex = 1;
  if length(varargin) < 2 && any(SHIFTS)
    % need to calculate the shift vectors. Assuming thatt the fft is
    % centered.
    preCalc = 0;
    flgComplexShift = 1;
    [dU,dV,dW] = BH_multi_gridCoordinates(size(IMAGE),'Cartesian',METHOD, ...
                                                                {'none'},1,1,0);
  elseif any(SHIFTS)
    preCalc = 1;
    flgComplexShift = 1;
  else
    preCalc = 0;
    flgComplexShift = 0;
  end     
  
  % Set the direction of the shift by the sign of the phase change. 
  if strcmpi(DIRECTION,'forward')
    phaseDir = -1;
     
  elseif strcmpi(DIRECTION, 'inv')
    phaseDir = 1;
  else
    error('DIRECTION in complex resample is what?');
  end
  
  % Shift prior to rotating if inv, and then set shifts to zero so that the
  % rotational grids are not shifted at all.
  phaseShifts = SHIFTS;
  SHIFTS = 0 .* SHIFTS;
  if strcmpi(DIRECTION,'inv')
    if ( flgComplexShift && preCalc)
      % Extra lines rather than re-assigning varargin to dU,dV,dW.
      IMAGE = IMAGE .* exp((phaseDir*2i*pi).*(varargin{2}{1}.*phaseShifts(1) + ...
                                              varargin{2}{2}.*phaseShifts(2) + ...
                                              varargin{2}{3}.*phaseShifts(3)));

    elseif ( flgComplexShift )
      IMAGE = IMAGE .* exp((phaseDir*2i*pi).*(dU.*phaseShifts(1) + ...
                                              dV.*phaseShifts(2) + ...
                                              dW.*phaseShifts(3)));
                                                
    end
  end
end
% For individual resampling, pushing to the gpu is (with the current hardware)
% slow enough that it negates the benefit. If however the supplied image is
% already on the GPU, then this is much faster. Careful not to pull from the GPU
% and then push back.
if (useGPU)
  IMAGE  = gpuArray(IMAGE);
end

flgSymmetry = 0;
symmetry = 1;
mag = 1;
if isa(CONVENTION, 'cell')
  interpMethod = CONVENTION{3};  
  symmetry = CONVENTION{2};
  if (symmetry > 1)
    flgSymmetry = true;
  end
  
  if length(CONVENTION) >= 4
    % Note that unlike BH_reScale3D, the image size is NOT changed here, which
    % could result in clipping.
    mag = CONVENTION{4};
  end
  
 
  if length(CONVENTION) == 5
    % a binary mask to limit the region searched
    volBinary = CONVENTION{5};
    flgMask = 1;
% % %     This is a wasted pre-allocation. Just zero outside interp mask in
% the output vol.
% % %     if (useGPU)
% % %       imOUT = zeros(size(IMAGE),'single','gpuArray');
% % %     else
% % %       imOUT = zeros(size(IMAGE),'single');
% % %     end
% % %     if (flgComplex)
% % %       imOUT = complex(imOUT);
% % %     end
  end
  
  CONVENTION = CONVENTION{1};


end


if strcmpi(DIRECTION,'forward') && (flgMask) && flgComplexShift
  % These are shifted after rotation but prior to reshaping, so
  % linearize the shift vectors; It probably makes sense to just move
  % the shift to be after reshaping.
  dU = dU(volBinary);
  dV = dV(volBinary);
  dW = dW(volBinary);
end

% Get transformed coords, if not already a matrix
if numel(ANGLES) == 3
  R = BH_defineMatrix(ANGLES, CONVENTION, DIRECTION);
elseif numel(ANGLES) == 9
  R = reshape(ANGLES, 3,3);
else
  error('ANGLES must be either three eulers or 9 rot matrix')
end

if (inputVectors)
  if (flgSymmetry)  

    [ Xnew,Ynew,Znew,x1,y1,z1 ] = BH_multi_gridCoordinates( size(IMAGE), 'Cartesian', METHOD, ...
                                                    {'single',R(:),SHIFTS',DIRECTION,symmetry,mag,volBinary},...
                                                    0, 1, 0, varargin{1} );
  else
    [ Xnew,Ynew,Znew,x1,y1,z1 ] = BH_multi_gridCoordinates( size(IMAGE), 'Cartesian', METHOD, ...
                                                    {'single',R(:),SHIFTS',DIRECTION,1,mag,volBinary},...
                                                    0, 1, 0, varargin{1});
  end  
else
  if (flgSymmetry)  

    [ Xnew,Ynew,Znew,x1,y1,z1 ] = BH_multi_gridCoordinates( size(IMAGE), 'Cartesian', METHOD, ...
                                                    {'single',R(:),SHIFTS',DIRECTION,symmetry,mag,volBinary},...
                                                    0, 1, 0 );
  else
    [ Xnew,Ynew,Znew,x1,y1,z1 ] = BH_multi_gridCoordinates( size(IMAGE), 'Cartesian', METHOD, ...
                                                    {'single',R(:),SHIFTS',DIRECTION,1,mag,volBinary},...
                                                    0, 1, 0 );
  end  
end


%%%%%%%%%%%% First interpolate the real part. If complex then also the
%%%%%%%%%%%% imaginary.
if (flgSymmetry)

  if strcmpi(interpMethod, 'spline')
     fgrid = griddedInterpolant({x1,y1,z1},real(IMAGE), 'spline', 'none');
     TRANS_IMAGE = fgrid(Xnew{1}, Ynew{1},Znew{1});
     
     TRANS_IMAGE(isnan(TRANS_IMAGE)) = 0;
  else
     TRANS_IMAGE = interpn(x1,y1,z1,real(IMAGE),Xnew{1},Ynew{1},Znew{1},'linear',0);
     
  end
else
  
  TRANS_IMAGE = interpn(x1,y1,z1,real(IMAGE),Xnew,Ynew,Znew,'linear',0);
end

if (flgSymmetry)
  symInc = 360/symmetry;
  
  for iSym = 2:symmetry

%     Rsym = R * BH_defineMatrix([0,0,iSym*symInc], CONVENTION, DIRECTION);
%     [ Xnew,Ynew,Znew,x1,y1,z1 ] = BH_multi_gridCoordinates( size(stackIN), ...
%                                                       'Cartesian', METHOD, ...
%                                                {Rsym(:);SHIFTS';DIRECTION},...
%                                                    0, 1, 0 );      
    if strcmpi(interpMethod, 'spline')       
       symVol =  fgrid(Xnew{iSym}, Ynew{iSym},Znew{iSym});
       
       symVol(isnan(symVol)) = 0;
       TRANS_IMAGE = TRANS_IMAGE + symVol;
    else
         
       TRANS_IMAGE = TRANS_IMAGE + ...
                          interpn(x1,y1,z1,real(IMAGE),Xnew{iSym},Ynew{iSym},Znew{iSym},'linear',0);
      
    end 
  end
  
TRANS_IMAGE = TRANS_IMAGE ./ symmetry; clear firstVol symVol fgrid 
end

%%%%%%%%%%%% Imaginary
% If clearing the interpolants sequentially for symmetry helps with real
% space OOM on large volumes, make sure to update this too.
if ( flgComplex )
  if (flgSymmetry)
    % if symmetry applied, return a cell, with the first being the asymmetric, and
    % second being the symmetrized volume
    if strcmpi(interpMethod, 'spline')
       fgrid = griddedInterpolant({x1,y1,z1},imag(IMAGE), 'spline', 'none');
       TRANS_IMAGE_Imag = fgrid(Xnew{1}, Ynew{1},Znew{1});
       TRANS_IMAGE_Imag(isnan(TRANS_IMAGE_Imag)) = 0;
    else
       TRANS_IMAGE_Imag = interpn(x1,y1,z1,imag(IMAGE),Xnew{1},Ynew{1},Znew{1},'linear',0);

    end
  else

    TRANS_IMAGE_Imag = interpn(x1,y1,z1,imag(IMAGE),Xnew,Ynew,Znew,'linear',0);
  end

  if (flgSymmetry)
    symInc = 360/symmetry;

    for iSym = 2:symmetry

  %     Rsym = R * BH_defineMatrix([0,0,iSym*symInc], CONVENTION, DIRECTION);
  %     [ Xnew,Ynew,Znew,x1,y1,z1 ] = BH_multi_gridCoordinates( size(stackIN), ...
  %                                                       'Cartesian', METHOD, ...
  %                                                {Rsym(:);SHIFTS';DIRECTION},...
  %                                                    0, 1, 0 );      
      if strcmpi(interpMethod, 'spline')       
         symVol =  fgrid(Xnew{iSym}, Ynew{iSym},Znew{iSym});
         symVol(isnan(symVol)) = 0;
         TRANS_IMAGE_Imag = TRANS_IMAGE_Imag + symVol;
      else

         TRANS_IMAGE_Imag = TRANS_IMAGE_Imag + ...
                            interpn(x1,y1,z1,imag(IMAGE),Xnew{iSym},Ynew{iSym},Znew{iSym},'linear',0);
      end 
    end

  TRANS_IMAGE_Imag = TRANS_IMAGE_Imag ./ symmetry; clear firstVol symVol fgrid 
  end
 
  %for iSymImg = 1:1+flgSymmetry

    TRANS_IMAGE = complex(TRANS_IMAGE,TRANS_IMAGE_Imag);
    TRANS_IMAGE_Imag = [];
  %end
 
  % A forward transform is shifted after rotation.
  if strcmpi(DIRECTION,'forward')
    for iSymImg = 1:1+flgSymmetry
      if ( flgComplexShift && preCalc)
        % Extra lines rather than re-assigning varargin to dU,dV,dW.
        TRANS_IMAGE = TRANS_IMAGE .* ...
                               exp((phaseDir*2i*pi).*(varargin{2}{1}.*phaseShifts(1) + ...
                                                varargin{2}{2}.*phaseShifts(2) + ...
                                                varargin{2}{3}.*phaseShifts(3)));

      elseif ( flgComplexShift )


        TRANS_IMAGE = TRANS_IMAGE .* ...
                               exp((phaseDir*2i*pi).*(dU.*phaseShifts(1) + ...
                                                dV.*phaseShifts(2) + ...
                                                dW.*phaseShifts(3)));

      end
    end
  end
end

% In alignRaw any ctf correction has been applied so we aren't worried about
% delocalized information, and so interpolate only within the minimum mask
% needed and not in the padded area. In this case, return only a single image
% and not a cell, and let it be the symmetric one if that is what is requested.
if (flgMask)

    TRANS_IMAGE(~volBinary) = 0;

end
dU = []; dV = []; dW = []; varargin = []; IMAGE = [];
Xnew = []; Ynew = []; Znew = []; volBinary = []; imOUT = [];




end % this is the end of resample3D




