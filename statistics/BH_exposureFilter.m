function [ exposureFilter ] = BH_exposureFilter( SIZE, TILT_GEOMETRY, METHOD,...
                                                 SAMPLING, SHIFT_ORIGIN, ...
                                                 varargin)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% Default behavior is to leave the origin alone (reciprocal space.)
if nargin == 3
  SHIFT_ORIGIN = 0;
end

% If making only a projection and not a stack, overwrite the projections index
% in the stack.
if size(TILT_GEOMETRY,1) == 1
  TILT_GEOMETRY(1) = 1;
end
if nargin < 6
  % The A-D are the defaults as published in the optimal exposure paper
  % optW controls the exponential fall-off when the current exposure
  % exceeds the optimal exposure. In unblur, this is "inf" as everything
  % beyond this point is set to zero. A small number still gives a steep
  % fall off, but with some taper.
  expA = 0.24499;
  expB =-1.66490;
  expC = 2.81410;
  optD = 2.51284;
  optW = 0.0;
else
  expA = varargin{1}(1);
  expB = varargin{1}(2);
  expC = varargin{1}(3);
  optD = varargin{1}(4);
  optW = varargin{1}(5);
end
nPrjs = size(TILT_GEOMETRY,1);

% Assumed constant across tilts
WAVELENGTH = TILT_GEOMETRY(1,18);

% Optimal exposure filter from Grant,Grigorieff 2015 eLife
% Assuming either 200 or 300 KV
if WAVELENGTH > 2.1*10^-12
  % 200KV which needs to be scaled since Grant's formula is fit for 300KV
  kvScale = 0.8;
else
  kvScale = 1.0;
end

% The only time it wouldn't be square is in template matching where the chunks
% aren't strictly cubic (a little less accurate but faster.)
if numel(SIZE) == 1
  SIZE = [SIZE,SIZE,nPrjs];
elseif numel(SIZE) == 2
  SIZE = [SIZE,nPrjs];
end

% This is from an "odd" behavior where zeros fails due to size being non-numeric.
% iT turns out SIZE cannot be a gpu array for some reason. Note in development docs
% and add check before removing.
%fprintf('size is numeric --> %d , size is %d %d %d\n',isnumeric(SIZE),SIZE);

clear exposureFilter


if strcmp(METHOD,'GPU')
  exposureFilter = zeros(SIZE, 'single', 'gpuArray');
elseif strcmp(METHOD, 'cpu')
  exposureFilter = zeros(SIZE, 'single');
else
    error('METHOD must be GPU or %s\n','cpu');
end


[criticalDose,~,~,~,~,~] =  BH_multi_gridCoordinates( ...
                                              SIZE(1:2),'Cartesian',...
                                              METHOD,{'none'},1,SHIFT_ORIGIN,1);

% Assuming the pixelSize is constant across projections, the only change is
% the cummulative dose, so precompute everything
pixelSize = TILT_GEOMETRY(1,16).*SAMPLING;
criticalDose =  kvScale.*(expA.* (criticalDose./pixelSize.*10^-10).^expB +expC);
optimalDose = (optD.*criticalDose);

% Precompute some values

% criticalDose = exp(-0.5.*criticalDose.^-1);
                  
                    
                       
for iPrj = 1:nPrjs 

  
  CUMeDOSE = TILT_GEOMETRY(iPrj,11);

  % Frequency where the dose exceeds the optimal dose (2.51*critical) are
  % set to zero in the unblur code .. use a gaussian falloff.
  
  
%%%%% This creates a more agressive dose filter perpendicular to the tilt axis.
%%%%% With these parameters I saw no change, better or worse which is a little
%%%%% odd. Maybe check it out in the future.
% % %   sX = SIZE(1);
% % %   sY = floor(SIZE(2).*(2-cosd(TILT_GEOMETRY(iPrj,4))));
% % %   [radialGrid,~,~,~,~,~] =  BH_multi_gridCoordinates( ...
% % %                                                 [sX,sY],'Cartesian',...
% % %                                                 METHOD,{'none'},1,1,1);
% % %   oX = ceil((sX+1)./2);    
% % %   oY = ceil((sY+1)./2);
% % %   radialGrid = radialGrid ./ (radialGrid(oX,oY+oX-2).*2);
% % %   trimVal = BH_multi_padVal([sX,sY],SIZE(1:2));
% % %   radialGrid = BH_padZeros3d(radialGrid,trimVal(1,:),trimVal(2,:),METHOD,'single');
% % %   if ~(SHIFT_ORIGIN)
% % %     radialGrid = ifftshift(radialGrid);
% % %   end
% exp((-0.5*CUMeDOSE)./criticalDose).* ...


    optimalMask = ( (optimalDose>=CUMeDOSE) + exp(-optW*(CUMeDOSE-optimalDose)).*(optimalDose<CUMeDOSE) );
    optimalMask(~isfinite(optimalMask)) = 1;
    exposureFilter(:,:,TILT_GEOMETRY(iPrj,1)) = exp((-0.5*CUMeDOSE)./criticalDose).*optimalMask;
                                                

% % %     exposureFilter(:,:,TILT_GEOMETRY(iPrj,1)) = exp((-0.5*CUMeDOSE)./criticalDose).* ...
% % %                                                 (optimalDose>CUMeDOSE); 
   
end


clear optimalMask criticalDose optimalDose

end
