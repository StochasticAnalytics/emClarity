function [ WEDGE_MASK, padValue] = BH_weightMask3d(SIZE, ORIENTATION, METHOD, ...
  particleRadius, flgIncCtf, ...
  SYMMETRY, samplingRate)
%Create a missing wedge mask.
%
%   Input variables:
%
%   SIZE = size of the padded fft of the image
%
%   ORIENTATION = missing wedge orientation
%       [tiltangles, sample orientation], [-a, a, e1, e2, e3]
%       As elsewhere the sample orientation are the ZXZ euler angles that
%       transform coordinates from the microscope to the sample frame.
%
%   METHOD = case sensitive 'GPU', cpu otherwise.
%
%   particleRadius = Thickness of particle in pixels in beam direction, used to
%               calculate a zone of influence for the mask.
%
%   Output variables:
%
%   WEDGEMASK = missing wedge mask ready to be applied to 3d fft
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Goals & Limitations:
%
%	As a test case, use an asymmetric wedge -50,70 in order to visualize
%	any ambiguities in angles.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   TODO:
%		- test gpu option for function and return value. (template search requires
%		return of gpu)
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Save the original size for binary wedge calc
binaryWedgeSize = abs(SIZE);

if isvector(ORIENTATION)
  flgCTF = 0;
  tiltColumn=1;
  tiltAzimuth=90;
  
else
  % Consider amplitude modulation by CTF (NOT CTF envelope though)
  flgCTF = 1;
  tiltColumn=4;
  tiltAzimuth=ORIENTATION(1,6);
  
end

flgTiltWeight = 1;
flgRadial = 1;
flgSmooth = 1;

if flgCTF && (flgIncCtf == 3)
  % Don't include exposure weight, just CTF
  flgNoExposure = 1;
elseif flgCTF && (flgIncCtf == 4)
  % Calc full CTF
  flgNoExposure = 0;
elseif flgCTF && (flgIncCtf == 2)
  % Calc just the weights, which are still different than a binary wedge, but
  % don't include CTF or exposure.
  flgCTF = 0;
end

if strcmp(METHOD, 'GPU')
  useGPU = true;
else
  useGPU = false;
end

%SIZE = [512,512,512];
if all(SIZE > 0)
  outputScaling = SIZE(1)/512;
  SIZE = [512,512,512]
  %   outputScaling = SIZE(1)/256;
  %   SIZE = [256,256,256]
else
  %optional override for template matching which may not always be cubic, wich
  %trades a little accuracy in the ctf mask for speed.
  SIZE = abs(SIZE)
  outputScaling = 1
end

if strcmpi(METHOD, 'applyMask')
  useGPU = true;
  filterMask = false;
  METHOD = 'GPU';
elseif strcmpi(METHOD, 'applyMaskcpu')
  filterMask = true;
  useGPU = false;
  METHOD = 'cpu';
  
elseif strcmp(METHOD,'binaryWedgeGPU') || strcmp(METHOD,'binaryWedgeCpu')
  if strcmp(METHOD,'binaryWedgeGPU')
    METHOD = 'GPU';
  elseif strcmp(METHOD,'binaryWedgeCpu')
    METHOD = 'cpu';
  else
    error('binary wedge uncertain if cpu or gpu.')
  end
  
  % Move binary wedge to a separate function
  maxTilt = 90-abs(min(ORIENTATION(:,tiltColumn)));
  minTilt = 90-abs(max(ORIENTATION(:,tiltColumn)));
  tiltAxis= BH_defineMatrix(1.*[90-tiltAzimuth,0,0], 'Bah', 'forwardVector');
  [bX,~,bZ,~,~,~] = BH_multi_gridCoordinates(binaryWedgeSize,'Cartesian',METHOD,{'single',tiltAxis,[0,0,0]','forwardVector',1,1},0,1,0);
  theta = atan2d(bX,bZ); clear bX bZ
  
  WEDGE_MASK = single(~( (180-minTilt<=theta | theta<=maxTilt-180) | ...
    (-1.*minTilt<=theta & theta<=maxTilt) ));
  
  % Set the value at the origin = 0.2
  origMask = ceil((size(WEDGE_MASK)+1)./2);
  WEDGE_MASK(origMask(1)-2:origMask(1)+2,:,origMask(3)-2:origMask(3)+2) = 1;
  
  
  clear theta maxTilt minTilt tiltAxis
  
  [ gaussKernel ] = BH_multi_gaussian3d(9, 1.5 );
  if strcmpi(METHOD,'GPU')
    gaussKernel = gpuArray(gaussKernel);
  end
  
  
  
  WEDGE_MASK = convn(WEDGE_MASK, gaussKernel, 'same');
  %
  %   WEDGE_MASK = convn(WEDGE_MASK, gaussKernel, 'same');
  
  %   [ WEDGE_MASK ] = BH_multi_randomizeTaper(WEDGE_MASK);
  
  rad = fftshift(BH_bandpass3d(size(WEDGE_MASK),0,0,0,METHOD,'nyquist'));
  WEDGE_MASK = WEDGE_MASK .* rad;
  WEDGE_MASK = (WEDGE_MASK - min(WEDGE_MASK(:)));
  WEDGE_MASK = WEDGE_MASK./(max(WEDGE_MASK(:)));  clear rad
  padValue = [0,0,0;0,0,0];
  return
else
  filterMask = false;
  
end







padValue = [0,0,0;0,0,0];
% Grid for the binary mask
% min / max tilt angle, inverted for reciprocal space



% Minimum size to pad in each dimension so that rotation will not cause
% extrapolation.

% % % padValue = ceil(0.5.*(sqrt(SIZE*SIZE')-SIZE));
% % % paddedSize = SIZE + 2.*padValue

% % % % 1/Diameter = thickness of projection in reciprocal space
% % % zoneOfInfluence = floor(SIZE(3)/(2.*THICKNESS))
% % % [~,~,Z] = ndgrid(zeros(paddedSize(1),1), ...
% % %                  zeros(paddedSize(2),1), ...
% % %                  window(@hamming,2.*zoneOfInfluence+1));

% Use the shape transform to estimate the extent of the projections influence
% Even though we end up multiplying, because the projection is only a single
% pixel thick to start with, this is equivalent to convolution, st this is
% essentially accounting for creating a finite backprojection volume by
% multiplication in real space with a box, hence this sinc based convolution in
% reciprocal space.lt

blob = BH_mask3d('sphere',SIZE, outputScaling.*[1,1,1].*particleRadius.*2,[0,0,0]);
if ~(useGPU)
  blob = gather(blob);
end
BLOB = abs(fftn(blob)); clear blob
BLOB = fftshift(BLOB ./ max(BLOB(:)));
oB = ceil((SIZE+1)./2);
blobKernel = BLOB(oB(1)-6:oB(1)+6,...
  oB(2)-6:oB(2)+6,...
  oB(3)-6:oB(3)+6);

BLOB = squeeze(gather(BLOB(oB(1), oB(2), :)))';

zWeight = repmat(BLOB, SIZE(2),1,SIZE(1));
zWeight = (permute(zWeight,[3,1,2]));


% For re-weighting, calc the r-weighting applied in tomo reconstruction, note
% this goes from 1 --> nX/2 and not 0 --> 0.5

if (flgRadial)
  
  [ rWeight ] = calc_rWeight( SIZE, 'single', METHOD);
  
  % % % rWeight = ((abs([-1*floor((SIZE(1))/2):0,1:floor((SIZE(1)-1)/2)])'));
  % % % rOrig = ceil((SIZE(1)+1)./2);
  % % % % % % % imod tilt zero freq = 0.2 * first non zero component
  % % % rWeight(rOrig ) = 1;
  % % %
  % % %
  % % % % rWeight = rWeight + 1;
  % % %
  % % %
  % % % % rWeight = (rWeight ./ max(rWeight)).^0.5;
  % % %  rWeight = rWeight + 1./rWeight.^2;
  % % % % resample2d only handles scaling right now, so pad to z=3
  % % % rWeight = repmat((rWeight), 1, SIZE(2),3);
  % % %
  % % % rWeight = BH_resample3d( rWeight, [90-tiltAzimuth,0,0], ...
  % % %                          [0,0,0],'Bah','GPU','forwardVector');
  % % % rWeight = rWeight(:,:,2);
  
else
  rWeight = 1;
end

[mtf,~,~,~,~,~] = BH_multi_gridCoordinates(SIZE(1:2), ...
  'Cartesian','cpu',...
  {'none'},1,1,1);

if (useGPU)
  [radialGrid,~,~,~,~,~] = BH_multi_gridCoordinates(SIZE(1:2), ...
    'Cartesian','GPU',...
    {'none'},1,0,1);
else
  [radialGrid,~,~,~,~,~] = BH_multi_gridCoordinates(SIZE(1:2), ...
    'Cartesian','cpu',...
    {'none'},1,0,1);
  
end




pixelSize = ORIENTATION(1,16).*ORIENTATION(1,14).*samplingRate; % scaled pixel size

radialGrid = radialGrid./(pixelSize.*10.^10);

%       exposureFilter = ones([SIZE(1:2),size(ORIENTATION,1)],'single');
if (flgCTF)
  if (flgNoExposure)
    exposureFilter = zeros(1,1,size(ORIENTATION,1)) + 1;
  else
    [ exposureFilter ] = BH_exposureFilter( SIZE(1:2), ORIENTATION, 'cpu',samplingRate,0 );
  end
  
  % Calc the downweighting due to CTF --> here assuming phases were flipped by
  % multiplying by the CTF past the first zero
  iPrj = 1;
  defocus = [ORIENTATION(iPrj,15) - ORIENTATION(iPrj,12), ...
    ORIENTATION(iPrj,15) + ORIENTATION(iPrj,12), ...
    ORIENTATION(iPrj,13)];
  Cs = ORIENTATION(iPrj,17);
  wavelength = ORIENTATION(iPrj,18);
  ampContrast = ORIENTATION(iPrj,19);
  
  % assuming ampContrast = 0.1, using -0.15 results in a weight with
  % (0.1^0.15)^2~ 0.5 at zero freqency. Allows some recovery of low freq without
  % creating too severe a blur
  [Hqz, HqzUnMod] = BH_ctfCalc(pixelSize,Cs,wavelength,defocus,SIZE(1:2),ampContrast,-1,-1);
  
  Hqz = single(abs(Hqz.*HqzUnMod));
  
  %  Hqz = conv2(Hqz,fspecial('gaussian',[5,5],1.0),'same');
  %    SAVE_IMG(MRCImage(gather(single(Hqz))), 'tmp.mrc');
  
else
  Hqz = zeros(size(radialGrid),'single')+1;
  exposureFilter = zeros(1,1,size(ORIENTATION,1)) + 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





% Make grid where 0.5 = Nyquist
if (useGPU)
  %centralSection = zeros(SIZE,'single', 'gpuArray');
  %centralSection = gpuArray(centralSection);
  wedgeMask = zeros(SIZE,'single', 'gpuArray');
else
  %centralSection = zeros(SIZE,'single');
  wedgeMask = zeros(SIZE,'single');
end

% centralSection(:,:,ceil((SIZE(3)+1)/2) - zoneOfInfluence: ...
%                    ceil((SIZE(3)+1)/2) + zoneOfInfluence) = Z;

% % % centralSection = centralSection .* zWeight;
% % % clear zWeight


nPrjs = size(ORIENTATION,1);
if nPrjs == 1
  ORIENTATION = ORIENTATION';
  nPrjs = size(ORIENTATION,1);
end

for iPrj = 1:nPrjs
  
  
  %exposureFilter = exp(-1.*ORIENTATION(iPrj,11).*(fftshift(radialGrid)).^2);
  
  % Applying a re-weighting to the ctfTiles - remove this
  % Optimal exposure filter from Grant,Grigorieff 2015 eLife
  % Assuming either 200 or 300 KV
  % % % WAVELENGTH = ORIENTATION(iPrj,18);
  % % % if WAVELENGTH > 2.1*10^-12
  % % %   kvScale = 0.8;
  % % % else
  % % %   kvScale = 1.0;
  % % % end
  % % % CUMeDOSE = ORIENTATION(iPrj,11);
  % % % exposureFilter = exp(-0.5*CUMeDOSE .* (kvScale.*0.245.*fftshift(radialGrid) .^ -1.665 + 2.81).^-1);
  
  %expF = 1;
  if  (flgCTF)
    expF = exposureFilter(:,:,ORIENTATION(iPrj,1));
  else
    expF = 1;
  end
  %expF = (expF./sqrt((sum(sum(sum(abs(expF).^2))))./numel(expF)));
  if (flgRadial)
    centralSection =  repmat(rWeight.* fftshift(expF.* ...
      Hqz ),1,1, SIZE(3)) .* ...
      zWeight;
    %       centralSection =  repmat(rWeight.*  ...
    %                              Hqz .* mtf,1,1, SIZE(3)) .* ...
    %                              zWeight;
  else
    centralSection =  repmat(fftshift(exposureFilter(:,:,ORIENTATION(iPrj,1)).*Hqz ),1,1, SIZE(3)) .* ...
      zWeight;
  end
  
  
  
  symInc = 360/SYMMETRY;
  for iSym = 1:SYMMETRY
    % Since the angles here are used bring the rotated projection back to standard
    % basis, they will also rotate the x,y plane to the plane perpendicular to the
    % projection direction in reciprocal space.
    if any(size(ORIENTATION) == 1)
      R = BH_defineMatrix([90, ORIENTATION(iPrj), -90], 'Bah', 'inv');
    else
      %       rSample = BH_defineMatrix((ORIENTATION(iPrj,8:10)), 'Bah', 'inv');
      %       rElevation = BH_defineMatrix([90, ORIENTATION(iPrj, 7), -90], 'Bah', 'inv');
      rTilt = BH_defineMatrix([1.*ORIENTATION(iPrj,6),1.*ORIENTATION(iPrj,4),-1*ORIENTATION(iPrj,6)],'Bah','inv');
      
      % Switch to imod means that the inPlane rotation has been applied to the
      % projections already
      %%%rInPlane = BH_defineMatrix(1.*[ORIENTATION(iPrj,5),0,0], 'Bah', 'inv');
      %       rInPlane = eye(3);
      
      
      %       R = rElevation*rTilt*rSample*rInPlane;
      R = rTilt;
    end
    
    R = R*BH_defineMatrix([(1-iSym)*symInc,0,0],'Bah','forward');
    %
    if (flgTiltWeight)
      %      tiltDiff = 1-cosd(ORIENTATION(iPrj,tiltColumn)).^1;
      %      tiltGrad = 1-(tiltDiff.*exp(-15.*fftshift(radialGrid).^2));
      %      tiltWeight = repmat(tiltGrad,1,1,size(centralSection,3));
      iAng = ORIENTATION(iPrj,tiltColumn);
      tiltWeight = ((exp(-10.*mtf.^(0.5+cosd(iAng).^2.5))+(0.6))./(1.6)).^sind(abs(iAng));
      tiltWeight = repmat(tiltWeight,1,1,size(centralSection,3));
      
      %       tiltWeight = cosd(ORIENTATION(iPrj,tiltColumn));
    else
      tiltWeight = 1;
    end
    
    
    
    wedgeMask  = wedgeMask + ...
      BH_resample3d(centralSection.*tiltWeight, R, ...
      [0,0,0], 'Bah','GPU','inv');
    
    
    
    % % %                   BH_resample3d(centralSection.*tiltWeight, R, [0,0,0], 'Bah','GPU','inv');
  end
end


% wedgeMask = wedgeMask .* mtf;
clear centralSection rWeight zWeight Hqz radialGrid mtf
if (flgSmooth)
  
  [ gaussKernel ] = gpuArray(BH_multi_gaussian3d(5, 0.75 ));
  
  wedgeMask = convn(wedgeMask, gaussKernel, 'same');
  % % % % %   wedgeMask = convn(wedgeMask, gaussKernel, 'same');
end



% wedgeMask = wedgeMask(padValue(1)+1 : end - padValue(1), ...
%                 padValue(2)+1 : end - padValue(2), ...
%                 padValue(3)+1 : end - padValue(3));


% wedgeMask = single(wedgeMask ./ max(wedgeMask(:)));
% wedgeMask = convn(wedgeMask, gaussKernel, 'same');
%
% wedgeMask = wedgeMask ./ max(wedgeMask(:));
%  [rad,~,~,~,~,~] = BH_multi_gridCoordinates(SIZE,'Cylindrical', ...
%                              'GPU',{BH_defineMatrix([0,90,90-tiltAzimuth],'Bah','forwardVector'),[0,0,0]','forward'},1,1,0);
% % %
% %rad = (rad < 0.5);
%  wedgeMask = wedgeMask .* rad;



WEDGE_MASK = BH_reScale3d( gather(single(wedgeMask ./ max(wedgeMask(:)))), ...
  '', sprintf('%f',outputScaling), METHOD);

rad = fftshift(BH_bandpass3d(size(WEDGE_MASK),0,0,0,'GPU','nyquist'));
WEDGE_MASK = WEDGE_MASK .* rad;
WEDGE_MASK = (WEDGE_MASK - min(WEDGE_MASK(:)));
WEDGE_MASK = WEDGE_MASK./(max(WEDGE_MASK(:)));

%    SAVE_IMG(MRCImage(gather(single(wedgeMask))),'wm.mrc');% .* (rad <0.5) ;
%     error('sdfs')
%  BINARY_WEDGE = BINARY_WEDGE .* (rad <0.5) ;
%  BINARY_WEDGE = convn(single(BINARY_WEDGE), BH_multi_gaussian3d(5, 1.0 ),'same');
% % %   BINARY_WEDGE = BINARY_WEDGE ./ max(BINARY_WEDGE(:)) .*rad;
%   clear wedgeMask centralSection radialGrid rad zWeight


clearvars -except WEDGE_MASK padValues

end

function  [ rWeight ] = calc_rWeight( SIZE, PRECISION, METHOD)

rWeight = ((abs([-1*floor((SIZE(1))/2):0,1:floor((SIZE(1)-1)/2)])'));
if strcmp(METHOD,'GPU')
  rWeight = gpuArray(rWeight);
end
rOrig = ceil((SIZE(1)+1)./2);
% % % % imod tilt zero freq = 0.2 * first non zero component
rWeight(rOrig ) = 0.2;
[rCut] = find(rWeight == floor(0.45*SIZE(1)));
pixelFallOff = rCut(1) ;
taperLow = 0.5+0.5.*cos((((1:pixelFallOff)).*pi)./(length((1:pixelFallOff+1))));

pixelFallOff = SIZE(1)-rCut(2)+1 ;
taperTop = 0.5+0.5.*cos((((1:pixelFallOff)).*pi)./(length((1:pixelFallOff+1))));
rWeight(1:rCut(1)) = rWeight(1:rCut(1)).*flip(taperLow)';
rWeight(rCut(2):end) = rWeight(rCut(2):end).*taperTop';


%rWeight = rWeight + 1./rWeight.^2;
% resample2d only handles scaling right now, so pad to z=3
rWeight = repmat((rWeight), 1, SIZE(2),1);
if strcmpi(PRECISION,'single')
  rWeight = single(rWeight);
else
  % This should be the default.
  rWeight = double(rWeight);
end
end
