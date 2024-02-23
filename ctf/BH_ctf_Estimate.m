function [  ] = BH_ctf_Estimate(varargin)

global bh_global_do_2d_fourier_interp;
!mkdir -p aliStacks
modLocal = false;
if length(varargin) > 3
  error('Too many input arguments');
else
  % PARAMETER_FILE, STACK_BASENAME, gpuIDX
  % Assumes that all required files are in ./fixedStacks
  PARAMETER_FILE = varargin{1};
  STACK_BASENAME = varargin{2};
  
  if length(varargin) == 3
    anglesSkipped = EMC_str2double(varargin{3});
    modLocal = true;
  else
    anglesSkipped = 0;
  end
  
  stackNameIN = sprintf('fixedStacks/%s.fixed', STACK_BASENAME);
  tltFile = sprintf('fixedStacks/%s.tlt', STACK_BASENAME);
  stackNameOUT = sprintf('%s_ali1', STACK_BASENAME);
  mapBack = sprintf('fixedStacks/%s', STACK_BASENAME);
  collectionORDER = sprintf('fixedStacks/%s.order',STACK_BASENAME);
  
end
emc = BH_parseParameterFile(PARAMETER_FILE);

gpuIDX = BH_multi_checkGPU(-1);
gDev = gpuDevice(gpuIDX);

% for trouble shooting downstream, default falsw
flgSkip = false;
% For trouble shooting downstream, default true
resample_stack = true; 
% slightly dampen lower resolution information that may overwhelm the CCC calc,
% but don't risk too much noise amplification. 1 = fit just the amplitude (not
% PS) 0.5 = take sqrt prior to normalizing.
cccScale = 1;


flgReScale = true;
mapBackPrfx= mapBack;


try
  rawTLT = load(tltFile);
catch
  error('The tilt angle file %s was not found.\n', tltFile);
end

if (modLocal)
  localName = sprintf('fixedStacks/%s.local',STACK_BASENAME);
  if exist(localName,'file')
    fprintf('\nModifying the local alignments\n');
    BH_trimIMODLocal(localName,length(rawTLT),anglesSkipped)
  end
end
skipFitting = 0;
try
  PHASE_PLATE_SHIFT = emc.('PHASE_PLATE_SHIFT').*pi
catch
  PHASE_PLATE_SHIFT = [0,0]
end
if sum(PHASE_PLATE_SHIFT)
  skipFitting = 1;
end

flgStandardOrdeDoCalc = 1;
try
  flgCosineDose = emc.('oneOverCosineDose');
  startingAngle = emc.('startingAngle');
  startingDirection = emc.('startingDirection');
  doseSymmetricIncrement = emc.('doseSymmetricIncrement');
  doseAtMinTilt = emc.('doseAtMinTilt');
  
  flgOldDose = 0;
  tltOrder = calc_dose_scheme(emc,rawTLT,anglesSkipped,PHASE_PLATE_SHIFT);
  
  
catch
  fprintf('\nFalling back on old dose specification through a *.order file\n\n');
  fprintf('Parameters flgCosineDose=(0/1 bool), \nstartingAngle=, \nstartingDirection=[pos/neg],\ndoseSymmetricIncrement=[0, or # tilts per sweep],\n doseAtMinTilt are needed for the new method.\n');
  pause(2);
  tltOrder = load(collectionORDER);
  flgOldDose = 1;
  if size(tltOrder,2) == 6
    flgStandardOrdeDoCalc = 0;
    flgSkip = 1;
  else
  end
end


scalePixelsBy = 1;

useGPU = 1;
METHOD = 'GPU';


WAVELENGTH = 10^-12*1226.39/sqrt(emc.VOLTAGE + 0.97845*10^-6*emc.VOLTAGE^2) ;

CUM_e_DOSE = 0;
% test astigmatism vals
flgAstigmatism = 1;
if (flgAstigmatism ~=1 && flgAstigmatism ~= 0)
  error('flgAstigmatism should be 0 or 1');
end
% Add sanity checks after further testing.
maxAstig = 2000e-10;
astigStep = 150e-10;
coarseAngStep = deg2rad(10);
fineAngStep = deg2rad(0.5);



eraseSigma = 3;

eraseRadius = ceil(1.2.*(emc.('beadDiameter')./emc.pixel_size_si.*0.5));
flgImodErase = 0



% Assuming that the first CTF zero is always less than this value
FIXED_FIRSTZERO =  emc.pixel_size_si / (70*10^-10) ;
highCutoff = emc.pixel_size_si/emc.('defCutOff');
% I still use the def for underfocus < 0 as this places the origin at the
% focal plan in the microscope rather than on the specimen. Which makes
% more sense to me.
defEST = emc.('defEstimate').*10^6
defWIN =  emc.('defWindow').*10^6
tiltRange = [-1];

backGroundBuffer = 0.9985;

try
  deltaZTolerance = emc.('deltaZTolerance');
catch
  deltaZTolerance = 100e-9;
end

try
  zShift = abs(emc.('zShift'));
catch
  zShift = 150e-9;
end

if abs(zShift) > 100e-7
  error('make sure your zShift values are of reasonable amounts (50-200nm)');
end

try
  maxNumberOfTiles = emc.('ctfMaxNumberOfTiles');
catch
  
  maxNumberOfTiles = 10000;
  
end

% Starting at +/- 100nm
deltaZTolerance = deltaZTolerance / emc.pixel_size_si;
% Use to check for proper gradient.
zShift = zShift / emc.pixel_size_si;

% Tile size & overlap
try
  tileSize = emc.('ctfTileSize');
catch
  tileSize = floor(680e-10 / emc.pixel_size_si);
end
tileOverlap = 2;

tileSize = tileSize + mod(tileSize,2);
% tileSize = max(tileSize, 384);
if (tileSize > 512)
  tileOverlap = tileOverlap * 2;
end
fprintf('Using a tile size of %d\n',tileSize);

overlap = floor(tileSize ./ tileOverlap);

% Size to padTile to should be even, large, and preferably a power of 2
try
  paddedSize = emc.('paddedSize');
catch
  paddedSize = 768;
end

padVAL = BH_multi_padVal([tileSize,tileSize], [paddedSize,paddedSize]);



if exist(stackNameIN, 'file')
  
  if ( flgOldDose )
    TLT = zeros(length(rawTLT),23);
    TLT(:,1) = 1:length(rawTLT);
    TLT(:,23) = 1:length(rawTLT);
    TLT(:,4) = rawTLT; clear rawTLT
    nSkipped = 0;
    
    if ~(flgStandardOrdeDoCalc)
      TLT(:,1) = tltOrder(:,1);
      TLT(:,12)  = (tltOrder(:,4)-tltOrder(:,5))./2 .* 10^-10;
      TLT(:,13)  =  tltOrder(:,6) .* (pi / 180);
      TLT(:,15)  = -1.*(tltOrder(:,4)+tltOrder(:,5))./2 .* 10^-10;
      TLT(:,14) = tltOrder(:,3);
      sorted_dose = sortrows(tltOrder,2);
      cummul_dose = cumsum(sorted_dose(:,3));
      TLT(sorted_dose(:,1),11) = cummul_dose;
    end
  else
    nSkipped = length(rawTLT) - size(tltOrder,1); % used to adjust header
    TLT = zeros(size(tltOrder,1),23);
    TLT(:,1)  = tltOrder(:,1);
    TLT(:,23) = tltOrder(:,5); % This was only the appropriate tilts are read in.
    TLT(:,4)  = tltOrder(:,2); clear rawTLT
    TLT(:,11) = tltOrder(:,3);
    TLT(:,14) = 1;
    TLT(:,19) = tltOrder(:,4);
  end
  
  
  
  [pathName,fileName,extension] = fileparts(stackNameIN);
  if isempty(pathName)
    pathName = '.';
  end
else
  fprintf('ignoring %s, because the file is not found.\n', stackNameIN)
  
end



% Make ctf directory to store diagnostic images
system(sprintf('mkdir -p %s/ctf', pathName));


iMrcObj = MRCImage(stackNameIN,0);

% The pixel size should be previously set correctly, but if it is not, then we
% must maintain whatever is there in case beads are to be erased. The model
% used for this process depends on the pixel size in the header when it was
% created in IMod alignment.

iHeader = getHeader(iMrcObj);

iPixelHeader = [iHeader.cellDimensionX/iHeader.nX .* scalePixelsBy, ...
  iHeader.cellDimensionY/iHeader.nY .* scalePixelsBy, ...
  iHeader.cellDimensionZ/iHeader.nZ];
% Reduce the Z dimension after pixel size is calculated
iHeader.nZ = iHeader.nZ - nSkipped;

iOriginHeader= [iHeader.xOrigin , ...
  iHeader.yOrigin , ...
  iHeader.zOrigin ] ./ scalePixelsBy;

d1 = iHeader.nX ; d2 = iHeader.nY ; d3 = iHeader.nZ;




% Copy with column for defocus = input to CTF correct
% saved as <filename>_ctf.tlt


flgReOrderMapBack = 0;

TLT(:,2:3) = repmat([0.00,0.00],size(TLT,1),1);
TLT(:,5:10) = repmat([0,90.0,1.0,0.0,0.0,1.0],size(TLT,1),1);
% Defocus will go at 15 - 12 and 13 currently unused.
TLT(:,16:18) = repmat([emc.pixel_size_si,emc.Cs,WAVELENGTH],size(TLT,1),1);
TLT(:,19) = TLT(:,19) + emc.AMPCONT;

oddSize = [d1,d2,d3] - (1-mod([d1,d2,d3],2));
TLT(:,20:22) = repmat(oddSize,size(TLT,1),1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Part of the switch to listing dose not bfactor to use the optimal exposure
% filter by Grant/Grigorieff -
nPrjs = size(TLT,1);

if ( flgOldDose && flgStandardOrdeDoCalc )
  if CUM_e_DOSE < 0
    % This is a dose per tilt at 0degress to be scaled by the cosine of the
    % tilt angle
    flgCosineDose = 1;
    exposure = abs(CUM_e_DOSE);
  else
    flgCosineDose = 0;
    exposure = CUM_e_DOSE./nPrjs;
  end
  
  totalExposure = 0;
  % If the fit tilt angles have moved alot, you may end up with duplicates,
  alreadyPicked = zeros(nPrjs,1,'single','gpuArray');
  largeVect = alreadyPicked + 1000;
  for iExposure = 1:nPrjs
    % find the projection angle most closley matching (
    
    [~,iTilt] = min((alreadyPicked.*largeVect)+(abs(TLT(:,4)-tltOrder(iExposure))));
    alreadyPicked(iTilt) = 1;
    if flgCosineDose == 0
      totalExposure = totalExposure + exposure;
    else
      totalExposure = totalExposure + exposure/cosd(TLT(iTilt,4));
    end
    TLT(iTilt,11) = totalExposure;
    TLT(iTilt,14) = exposure;
    
  end
  
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ( resample_stack)

fprintf('Combining tranformations\n\n');
% Load in the mapBack alignment
mbEST = load(sprintf('%s.xf',mapBackPrfx));
mbTLT = load(sprintf('%s.tlt',mapBackPrfx));

outputStackName = sprintf('aliStacks/%s%s',stackNameOUT,extension);

try
  erase_beads_after_ctf = emc.('erase_beads_after_ctf');
catch
  erase_beads_after_ctf = false;
end

if (erase_beads_after_ctf)
  flgEraseBeads = 0;
else
  if exist(sprintf('%s.erase',mapBackPrfx),'file')
    flgEraseBeads = 1;
  else
    flgEraseBeads = 0;
    fprintf('\nDid not find the gold bead file (%s) for erasing, will skip\n\n',sprintf('%s.erase',mapBackPrfx));
  end
end

sizeCropped = [d1,d2,d3]-(1-mod([d1,d2,d3],2));
sizeCropped(3) = d3;

STACK = zeros(sizeCropped,'single');
samplingMaskStack = zeros(sizeCropped,'single');



if (flgReOrderMapBack)
  TLT = sortrows(TLT,1);
end

% if any([d1,d2] > 4096)
%   shiftMETHOD = 'cpu';
%   fprintf('transforming on cpu b/c > 4096\n')
% else
shiftMETHOD = 'GPU';
% end

% Redefine d3 incase views are ignored
d3 = size(TLT,1);

osX = 1-mod(d1,2); osY = 1-mod(d2,2);




for i = 1:d3
  %  fprintf('Transforming prj %d in fourier space oversampled by 2x physical Nyquist\n',i);
  
  
  % Stored in row order as output by imod, st transpose is needed. Inversion
  % of the xform is handled in resample2d.
  origXF = [1,0;0,1];
  
  newXF = reshape(mbEST(TLT(i,23),1:4),2,2)';
  
  
  dXYZ  = [(newXF*TLT(i,2:3)')' + mbEST(TLT(i,23),5:6) , 0];
  TLT(i,2:3) = dXYZ(1:2);
  dXYZ = dXYZ ./ scalePixelsBy;
  
  combinedXF = reshape((newXF*origXF)',1,4);
  TLT(i,7:10) = combinedXF;
  
  
  
  
  
  sizeODD = [d1,d2]-[osX,osY];
  
  % If it is even sized, shift up one pixel so that the origin is in the middle
  % of the odd output here we can just read it in this way, unlike super res.
  
  iProjection = ...
    single(getVolume(iMrcObj,[1+osX,d1],[1+osY,d2],TLT(i,23),'keep'));
  
  iProjection = real(ifftn(fftn(iProjection).* BH_bandpass3d(1.*[d1-osX,d2-osY,1],0,0,0,'GPU','nyquistHigh')));
  
  largeOutliersMean= mean(iProjection(:));
  largeOutliersSTD = std(iProjection(:));
  largeOutliersIDX = (iProjection < largeOutliersMean - 6*largeOutliersSTD | ...
    iProjection > largeOutliersMean + 6*largeOutliersSTD);
  iProjection(largeOutliersIDX) = (3*largeOutliersSTD).*randn([gather(sum(largeOutliersIDX(:))),1],'single');
  
  largeOutliersIDX = [];
  
  
  % Padding to avoid interpolation artifacts. For K3 images this can push
  % a 2080 close to or over the limit, so it is been reduced to 1/4 (from
  % 1) i.e. the image is paded to 1.25 x unless useFourierInterp is set >
  % 1;
  sizeSQ = floor(([1,1]+bh_global_do_2d_fourier_interp*0.25).*max(sizeODD));
  %         sizeSQ = floor(([1,1]).*max(sizeODD));
  
  padVal  = BH_multi_padVal(sizeODD,sizeSQ);
  trimVal = BH_multi_padVal(sizeSQ,sizeCropped(1:2));
  
  
  iProjection = iProjection - mean(iProjection(:));
  

  iProjection = BH_padZeros3d(iProjection,padVal(1,:),padVal(2,:), ...
    shiftMETHOD,'singleTaper');
  
  if (i == 1 && bh_global_do_2d_fourier_interp)
    bhF = fourierTransformer(iProjection,'OddSizeOversampled');
  end
  
  
  
  % Do the phase shift after rotating - need to invert the scaling since
  % we are in reciprocal space
  [imodMAG, imodStretch, imodSkewAngle, imodRot] = ...
    BH_decomposeIMODxf(combinedXF);
  
  
  
  if (bh_global_do_2d_fourier_interp)
    %   combinedInverted = BH_defineMatrix([imodRot,0,0],'Bah','forward').*(1/imodMAG);
    combinedInverted = BH_defineMatrix([imodRot,0,0],'Bah','forward');
    combinedInverted = combinedInverted([1,2,4,5]);
    
    iProjection = BH_resample2d(iProjection,combinedInverted,dXYZ(1:2), ...
    'Bah','GPU','forward',imodMAG,size(iProjection),bhF);
  else
    combinedInverted = BH_defineMatrix([imodRot,0,0],'Bah','forward').*(imodMAG);
    combinedInverted = combinedInverted([1,2,4,5]);
    iProjection = BH_resample2d(iProjection,combinedInverted,dXYZ(1:2),'Bah','GPU','forward',1.0,size(iProjection));
  end
  
  iSamplingMask = BH_resample2d(ones(sizeCropped(1:2),'single','gpuArray'),combinedXF,dXYZ(1:2),'Bah','GPU','forward',1.0,sizeCropped(1:2),NaN);
  
  iSamplingMask(isnan(iSamplingMask(:))) = 0;
  samplingMaskStack(:,:,i) = (gather(real(iSamplingMask)));
  iSamplingMask = [];
  
  % % % % %    iProjection = real(fftshift(ifftn(ifftshift(iProjection))));
  STACK(:,:,i)  = gather(real(BH_padZeros3d(iProjection, ...
    trimVal(1,:),trimVal(2,:),...
    shiftMETHOD,'single')));
  
  
end

if ( flgEraseBeads )
  STACK = BH_eraseBeads(STACK,eraseRadius, fileName, scalePixelsBy,0,sortrows(TLT,1));
end

[ STACK ] = BH_multi_loadAndMaskStack(STACK,TLT,'',100,emc.pixel_size_si*10^10,samplingMaskStack);


SAVE_IMG(MRCImage(STACK),outputStackName,iPixelHeader,iOriginHeader);
SAVE_IMG(MRCImage(samplingMaskStack),sprintf('%s.samplingMask',outputStackName),iPixelHeader,iOriginHeader);

end

if ~(flgSkip)
  
  gpuDevice(gpuIDX)
  [d1,d2,d3] = size(STACK);
  if (emc.pixel_size_si*10^10 < 0)
    
    flgCrop = 1;
    [croppedIMG,pixelOUT] = cropIMG(STACK(:,:,1),emc.pixel_size_si*10^10);
    [d1C,d2C] = size(croppedIMG); clear croppedIMG
    tltForExp = TLT;
    tltForExp(:,16) = pixelOUT*10^-10;
    pixelOUT
    % Redefining things down hear is a stupid thing to do. Fix this if you
    % keep the optino for cropping.
    FIXED_FIRSTZERO =  pixelOUT / 70 ;
    highCutoff = (pixelOUT*10^-10)/emc.('defCutOff');
    
  else
    
    flgCrop = 0;
    pixelOUT = emc.pixel_size_si*10^10;
    d1C = d1;
    d2C = d2;
    tltForExp = TLT;
    
  end
  
  d3 = size(STACK,3); 
  
  [radialForCTF,phi,~,~,~,~] = ...
    BH_multi_gridCoordinates([paddedSize,paddedSize,1],'Cylindrical','GPU', ...
    {'none'},1,0,0);
  
  radialForCTF = {radialForCTF./(pixelOUT.*10^-10),0,phi}; clear phi
  
  flgExpFilter = 0;
  
  inc = (0.5 - FIXED_FIRSTZERO) / (paddedSize/2);
  freqVector = [inc+FIXED_FIRSTZERO:inc:0.5 ];
  
  tic
  nT =  1;
  nT2 = 0;
  nT3 = 0;
  
  halfX = emc_get_origin_index(paddedSize);
  % % % % psTile = zeros([(paddedSize).*[1,1],3],'single','gpuArray');
  psTile = zeros([halfX,paddedSize,3],'single','gpuArray');
  
  bhF2 = fourierTransformer(randn(paddedSize,paddedSize,'single','gpuArray'));
  
  for k = 1:d3
    if (skipFitting)
      break
    end
    
    tiltIDX = TLT(k,1);
    % Center the pixel coordinates
    iEvalMask = BH_multi_gridCoordinates([d1C,1,1],'Cartesian','GPU',{'none'},0,1,0);
    
    % Convert to the z-height in the projection
    iEvalMask = iEvalMask.*(-1.*tand(TLT(tiltIDX,4)));
    
    iEvalPos = iEvalMask;
    iEvalNeg = iEvalMask;
    
    % Shift by any amount wanted
    iEvalPos = iEvalPos - zShift;
    iEvalNeg = iEvalNeg + zShift;
    
    % Select region limited by tolerance
    iEvalPos = ( iEvalPos > gpuArray(-deltaZTolerance) & iEvalPos < gpuArray(deltaZTolerance));
    iEvalNeg = ( iEvalNeg > gpuArray(-deltaZTolerance) & iEvalNeg < gpuArray(deltaZTolerance));
    iEvalMask = ( iEvalMask > gpuArray(-deltaZTolerance) & iEvalMask < gpuArray(deltaZTolerance));
    
    
    tmpTile = zeros([halfX,paddedSize,3],'single','gpuArray');
    
    % % % %    tmpTile = zeros([paddedSize.*[1,1],3],'single','gpuArray');
    
    if flgCrop
      [iProjection,~] = cropIMG(gpuArray(STACK(:,:,TLT(k,1))),emc.pixel_size_si*10^10);
    else
      iProjection = (gpuArray(STACK(:,:,TLT(k,1))));
    end
    
    iProjection = iProjection - ...
      BH_movingAverage(iProjection,[tileSize,tileSize]);
    
    iProjection = iProjection ./ ...
      BH_movingRMS(iProjection,[tileSize,tileSize]);
    
    for i = 1+tileSize/2:overlap:d1C-tileSize/2
      if min([nT,nT2,nT3])< maxNumberOfTiles && (iEvalMask(i) || iEvalPos(i) || iEvalNeg(i))
        for j = 1+tileSize/2:overlap:d2C-tileSize/2
          
          thisTile = abs(bhF2.fwdFFT(BH_padZeros3d(...
            (iProjection( ...
            i-tileSize/2+1:i+tileSize/2,...
            j-tileSize/2+1:j+tileSize/2)),...
            padVAL(1,:),padVAL(2,:),...
            'GPU','singleTaper')));

          tmpTile(:,:,1) = tmpTile(:,:,1) + thisTile;
          
          if (iEvalMask(i))
            nT = nT+1;
            tmpTile(:,:,1) = tmpTile(:,:,1) + thisTile;
          end
          if ( iEvalPos(i) )
            nT2 = nT2+1;
            tmpTile(:,:,2) = tmpTile(:,:,2) + thisTile;
          end
          if (iEvalNeg(i) )
            nT3 = nT3+1;
            tmpTile(:,:,3) = tmpTile(:,:,3) + thisTile;
          end
          
        end
      end
    end
    fprintf('%d tiles at dZ= 0\t%d tiles at dZ > 0\t%d tiles at dZ < 0, after tilt %d\n',nT,nT2,nT3,k);
    
    % Apply the dose filter to the sum of each projection to save a bunch of
    % multiplicaiton
    psTile = psTile + tmpTile;
    
  end
  clear tmpTile
  toc
  
  rotAvgPowerSpec = zeros([paddedSize,paddedSize,3],'single','gpuArray');
  for iTile = 1:3
    tmp =  bhF2.swapIndexFWD(psTile(:,:,iTile));
    psTile(:,:,iTile) = bhF2.swapIndexFWD(psTile(:,:,iTile));
    rotAvgPowerSpec(:,:,iTile) = BH_multi_makeHermitian(psTile(:,:,iTile),[paddedSize,paddedSize],1);
  end
  
  clear psTile
  
  if ~(skipFitting)
    % % % %   for iTile = 1:3
    % % % %     rotAvgPowerSpec(:,:,iTile) = (fftshift(rotAvgPowerSpec(:,:,iTile)));
    % % % %   end
    AvgPowerSpec = rotAvgPowerSpec;
    
    % TODO make a better rotational averaging funciton
    [rot1, rot2, ~, r1,r2, ~] = BH_multi_gridCoordinates(paddedSize.*[1,1], ...
      'Cartesian','GPU', ...
      {'none'},0,1,0);
    
    for i = 0.5:0.5:360
      R = BH_defineMatrix([i,0,0],'Bah','forward');
      ROT1 = R(1).*rot1 + R(4).*rot2;
      ROT2 = R(2).*rot1 + R(5).*rot2;
      
      for iTile = 1:3
        rotAvgPowerSpec(:,:,iTile) = rotAvgPowerSpec(:,:,iTile) + ...
          interpn(r1,r2,AvgPowerSpec(:,:,iTile),...
          ROT1,ROT2,'linear',0);
      end
    end
    clear ROT1 ROT2
    % Normalize on avgerage #, doseWeighting, and a radial filter to account
    % for rotational averaging.
    rotAvgPowerSpec = rotAvgPowerSpec ./ (720); clear a
    % rotAvgPowerSpec = rotAvgPowerSpec ./ (720.*(sqrt(fftshift(radialForCTF{1}.*(pixelOUT.*10^-10))))); clear a
    
    
    is_a_bummer = ~isfinite(rotAvgPowerSpec);
    if sum(is_a_bummer,'all') > 0.5*numel(rotAvgPowerSpec)
      error('the rotated Avg power spectrum is more than half nan or inf');
    else
      rotAvgPowerSpec(is_a_bummer) = 0;
    end
    
    
    currentDefocusEst = defEST;
    currentDefocusWin = defWIN;
    measuredVsExpected = zeros(2,3);
  end % end of not skip fitting
  
  for iTilt = 1:3
    
    if (skipFitting && resample_stack)
      currentDefocusEst = defEST;
      
      % Add the determined defocus, and write out with mic paramters as well.
      TLT(:,15) = repmat(-1.*defEST,size(TLT,1),1);
      
      [~, idx] = sortrows(abs(TLT(:,4)), -1);
      TLT = TLT(idx,:);
      % number in stack, dx, dy, tilt angle, projection rotation, tilt azimuth, tilt
      % elevation, e1,e2,e3, dose number (order in tilt collection), offsetX, offsetY
      % scaleFactor, defocus, emc.pixel_size_si, CS, Wavelength, Amplitude contrast
      fileID = fopen(sprintf('%s/ctf/%s_ctf.tlt',pathName,stackNameOUT), 'w');
      fprintf(fileID,['%d\t%08.2f\t%08.2f\t%07.3f\t%07.3f\t%07.3f\t%07.7f\t%07.7f\t',...
        '%07.7f\t%07.7f\t%5e\t%5e\t%5e\t%7e\t%5e\t%5e\t%5e\t%5e\t%5e\t',...
        '%d\t%d\t%d\t%8.2f\n'], TLT');
      fclose(fileID);
      fprintf('\n\nUsing the estimated value for defocus and phase shift provided\n\n');
      return;
    end
    
    radialAvg = [rotAvgPowerSpec((paddedSize/2)+1:end,(paddedSize/2)+1,iTilt)]';
    
    
    defRange = [currentDefocusEst-currentDefocusWin,currentDefocusEst+currentDefocusWin];
    
    defInc   = 0.01;
    
    defVal = (defRange(1):defInc:defRange(2))';
    
    cccStorage = zeros(length(defVal),2);
    cccStorage(:,1) = defVal;
    nDF = 1;
    for iDF = defVal'
      DF = iDF*10^-6;
      
      % TODO add a global switch for the damping
      
      if (emc.pixel_size_si < 1*10^-10)
        [ Hqz ] = BH_ctfCalc(radialForCTF,emc.Cs,WAVELENGTH,DF,paddedSize,-emc.AMPCONT,-1.0,-1);
      else
        [ Hqz ] = BH_ctfCalc(radialForCTF,emc.Cs,WAVELENGTH,DF,paddedSize,-emc.AMPCONT,-1.0);
      end
      
      try
        [ bg,  bandpass, rV ] = prepare_spectrum( Hqz, highCutoff, freqVector, radialAvg, 0);
      catch
        %       figure, imshow3D(gather(Hqz));
        %       figure, imshow3D(gather(rotAvgPowerSpec));
        %       highCutoff
        %       figure, plot(freqVector);
        %       figure, plot(radialAvg);
      end
      
      [ iCCC ] = calc_CCC( freqVector, bg,  bandpass, radialAvg, rV, cccScale);
      
      cccStorage(nDF, 2) = gather(iCCC);
      nDF = nDF +1;
    end
    
    [~,maxVal] = max(cccStorage(:,2));
    maxDef = cccStorage(maxVal,1);
    
    if (iTilt == 1)
      % Only save for the "true" defocus at the tilt-axes
      figure('Visible','off'), scatter(cccStorage(:,1), cccStorage(:,2));
      title(sprintf('CCC\nmax = %03.3f micron', maxDef));
      xlabel('defocus (micron)'); ylabel('CCC');
      
      saveas(gcf,sprintf('%s/ctf/%s_ccFIT.pdf',pathName,stackNameOUT), 'pdf')
    end
    DF = maxDef*10^-6;
    
    if (emc.pixel_size_si < 1*10^-10)
      [ Hqz ] = BH_ctfCalc(radialForCTF,emc.Cs,WAVELENGTH,DF,paddedSize,-emc.AMPCONT,-1.0,-1);
    else
      [ Hqz ] = BH_ctfCalc(radialForCTF,emc.Cs,WAVELENGTH,DF,paddedSize,-emc.AMPCONT,-1.0);
    end
    
    [ bg, bandpass, rV ] = prepare_spectrum( Hqz, highCutoff, freqVector, radialAvg, 0);
    
    
    if (iTilt == 1)
      % Only save for the "true" defocus at the tilt-axes
      figure('Visible','off'), plot(freqVector(bandpass),backGroundBuffer.*bg(freqVector(bandpass)),freqVector(bandpass),abs(radialAvg(bandpass)));
      title('Background fitting')
      
      saveas(gcf,sprintf('%s/ctf/%s_bgFit.pdf',pathName,stackNameOUT), 'pdf')
    end
    
    
    clear STACK exposureFilter
    
    pdfOUT = sprintf('%s/ctf/%s_psRadial_%d.pdf',pathName,stackNameOUT,iTilt);
    
    
    bgSubPS = (abs(radialAvg) - bg(freqVector)').*bandpass;
    bgSubPS = bgSubPS ./ max(bgSubPS(:));
    figure('Visible','off'), plot(freqVector(bandpass)./(pixelOUT),bgSubPS(bandpass), freqVector(bandpass)./(pixelOUT),abs(rV(bandpass)).^2./max(abs(rV(bandpass)).^2),'-g');
    title(sprintf('CTF fit\n%03.3f μm ', maxDef));
    xlabel('Spatial Frequency (1/Å)'); ylabel('Relative Power');
    
    saveas(gcf,pdfOUT, 'pdf')
    
    
    
    if (flgAstigmatism)
      
      radialForCTF = {fftshift(radialForCTF{1}),1,fftshift(radialForCTF{3})};
      [radialAstig,~,~,~,~,~] = ...
        BH_multi_gridCoordinates(size(Hqz),'Cartesian',...
        'GPU',{'none'},1,1,1);
      
      % Hqz from max defocus still exisists
      
      [ bg, bandpass, ~ ] = prepare_spectrum( Hqz, highCutoff ,...
        freqVector, AvgPowerSpec(:,:,iTilt), radialAstig);
      
      [ bgSubPS ] = calc_CCC( radialAstig, bg, bandpass, AvgPowerSpec(:,:,iTilt),-9999, cccScale) ;
      
      SAVE_IMG(MRCImage(gather(bgSubPS)), ...
        sprintf('%s/ctf/%s_avgPS-bgSub%s',pathName,fileName,extension));
      
      
      
      
      SAVE_IMG(MRCImage(gather(AvgPowerSpec)), ...
        sprintf('%s/ctf/%s_avgPS%s',pathName,fileName,extension));
      
      
      coarseDefSearch = 0:floor(maxAstig/astigStep);
      coarseAngSearch = -pi/2:coarseAngStep:pi/2;
      
      fineAngSearch = -1*coarseAngStep/2:fineAngStep:coarseAngStep/2;
      fineDefSearch = -astigStep/2:astigStep/4:astigStep/2;
      
      
      
      initAstigCCC = zeros(length(coarseDefSearch)*length(coarseAngSearch),3, 'gpuArray');
      
      
      n=1;
      
      for iAng = coarseAngSearch
        for iDelDF = coarseDefSearch
          df1 =  maxDef*10^-6 + iDelDF*astigStep;
          df2 =  maxDef*10^-6 - iDelDF*astigStep;
          
          if (emc.pixel_size_si < 1*10^-10)
            
            [ Hqz ] = BH_ctfCalc(radialForCTF,emc.Cs,WAVELENGTH, ...
              [df1,df2,iAng],size(radialForCTF{1}),-emc.AMPCONT,-1.0,-1);
          else
            [ Hqz ] = BH_ctfCalc(radialForCTF,emc.Cs,WAVELENGTH, ...
              [df1,df2,iAng],size(radialForCTF{1}),-emc.AMPCONT,-1.0);
          end
          
          
          [ iCCC ] = calc_CCC( radialAstig, bgSubPS, bandpass, AvgPowerSpec, Hqz, cccScale);
          
          initAstigCCC(n,:) = [iAng,iDelDF*astigStep,iCCC];
          n = n + 1;
        end
      end
      
      nPeaks = 1;
      top3 = zeros(nPeaks,3,'gpuArray');
      
      for iCCC = 1:nPeaks
        [~,c] = max(initAstigCCC(:,3));
        top3(iCCC,:) = initAstigCCC(c,:);
        initAstigCCC = initAstigCCC(initAstigCCC(:,1)~=initAstigCCC(c,1),:);
      end
      
      refineCCC = zeros(length(fineDefSearch)*length(fineAngSearch)*nPeaks,3,'gpuArray');
      
      n=1;
      for iPeak = 1:nPeaks
        mAng = top3(iPeak,1);
        mDef = top3(iPeak,2);
        
        for iAng = fineAngSearch
          for iDelDF = fineDefSearch
            
            df1 =  maxDef*10^-6 + (mDef + iDelDF);
            df2 =  maxDef*10^-6 - (mDef + iDelDF);
            
            % For values very close to zero, the search range may include
            % values |df1| < |df2| which is against convention.
            if abs(df1) >= abs(df2)
              
              if (emc.pixel_size_si < 1*10^-10)
                
                [ Hqz ] = BH_ctfCalc(radialForCTF,emc.Cs,WAVELENGTH, ...
                  [df1,df2,iAng+mAng], ...
                  size(radialForCTF{1}), -emc.AMPCONT,-1.0,-1);
              else
                [ Hqz ] = BH_ctfCalc(radialForCTF,emc.Cs,WAVELENGTH, ...
                  [df1,df2,iAng+mAng], ...
                  size(radialForCTF{1}), -emc.AMPCONT,-1.0);
              end
              
              
              [ iCCC ] = calc_CCC( radialAstig,bgSubPS, bandpass, ...
                AvgPowerSpec, Hqz,cccScale );
              
            else
              iCCC = -9999;
            end
            
            refineCCC(n,:) = [iAng+mAng,mDef + iDelDF,iCCC];
            n = n + 1;
            fprintf('%d / %d fine astigmatism search\n',n,size(refineCCC,1));
          end
        end
      end
      figure, plot(refineCCC(:,1), refineCCC(:,3), 'bo');
      [~,c] = max(refineCCC(:,3));
      save( sprintf('%s/ctf/%s_CCC.mat',pathName,fileName),'refineCCC','initAstigCCC');
      topScore = fopen(sprintf('%s/ctf/%s_astig.txt',pathName,fileName),'w');
      fprintf(topScore,'%7.7e %7.7e %2.7f\n',refineCCC(c,:));
      fclose(topScore);
      
      
      
      
    end
    
    if ( iTilt == 1)
      
      radialForCTF = {fftshift(radialForCTF{1}),1,fftshift(radialForCTF{3})};
      currentDefocusEst = maxDef;
      currentDefocusWin = (defWIN*.25);
      measuredVsExpected(1,:) = [maxDef + zShift*emc.pixel_size_si*10^6, maxDef, maxDef - zShift*emc.pixel_size_si*10^6];
      measuredVsExpected(2,2) = maxDef;
      % Add the determined defocus, and write out with mic paramters as well.
      TLT(:,15) = repmat(maxDef*10^-6,size(TLT,1),1);
      if (flgAstigmatism) && (refineCCC(c,3)~=-9999)
        TLT(:,12) = repmat(gather(refineCCC(c,2)),size(TLT,1),1);
        TLT(:,13) = repmat(gather(refineCCC(c,1)),size(TLT,1),1);
      end
      
      % Turn off astigmatism and restrict search range for handedness check.
      flgAstigmatism = 0;
      
      % Sort descending along the magnitude of the tilt angles because higher tilts take
      % longer on CTF correction. If more processor available than projections,
      % this doesn't affect anything.
      [~, idx] = sortrows(abs(TLT(:,4)), -1);
      TLT = TLT(idx,:);
      % number in stack, dx, dy, tilt angle, projection rotation, tilt azimuth, tilt
      % elevation, e1,e2,e3, dose number (order in tilt collection), offsetX, offsetY
      % scaleFactor, defocus, emc.pixel_size_si, CS, Wavelength, Amplitude contrast
      fileID = fopen(sprintf('%s/ctf/%s_ctf.tlt',pathName,stackNameOUT), 'w');
      fprintf(fileID,['%d\t%08.2f\t%08.2f\t%07.3f\t%07.3f\t%07.3f\t%07.7f\t%07.7f\t',...
        '%07.7f\t%07.7f\t%5e\t%5e\t%5e\t%7e\t%5e\t%5e\t%5e\t%5e\t%5e\t',...
        '%d\t%d\t%d\t%8.2f\n'], TLT');
      fclose(fileID);
      
    elseif iTilt == 2
      measuredVsExpected(2,1) = maxDef;
    elseif iTilt == 3
      measuredVsExpected(2,3) = maxDef;
    end % Stuff we only do on the full determin (tilt1)
    
  end % Loop on handedness check
  if sum(abs(diff(measuredVsExpected,1))) > sum(abs(measuredVsExpected(1,:) - flip(measuredVsExpected(2,:))))
    warnInvertedHand = 1;
  else
    warnInvertedHand = 0;
  end
  
  fprintf('\n******************************************************\n\n');
  fprintf('\nCloser to focus |\tAt focus |\tFarther from focus\n\n');
  fprintf('Expected defocus %3.2f %3.2f %3.2f\n\n', abs(measuredVsExpected(1,:)));
  fprintf('Measured defocus %3.2f %3.2f %3.2f\n\n' ,abs(measuredVsExpected(2,:)));
  if ( warnInvertedHand )
    fprintf('\nIt looks like your handedness may be inverted!!\n');
  else
    fprintf('\nIt looks like your handedness is probably correct.\n');
  end
  fprintf('\n******************************************************\n\n\n');
  
else
  
  if (resample_stack)
    [~, idx] = sortrows(abs(TLT(:,4)), -1);
    TLT = TLT(idx,:);
    % number in stack, dx, dy, tilt angle, projection rotation, tilt azimuth, tilt
    % elevation, e1,e2,e3, dose number (order in tilt collection), offsetX, offsetY
    % scaleFactor, defocus, emc.pixel_size_si, CS, Wavelength, Amplitude contrast
    fileID = fopen(sprintf('%s/ctf/%s_ctf.tlt',pathName,stackNameOUT), 'w');
    fprintf(fileID,['%d\t%08.2f\t%08.2f\t%07.3f\t%07.3f\t%07.3f\t%07.7f\t%07.7f\t',...
      '%07.7f\t%07.7f\t%5e\t%5e\t%5e\t%7e\t%5e\t%5e\t%5e\t%5e\t%5e\t',...
      '%d\t%d\t%d\t%8.2f\n'], TLT');
    fclose(fileID);
  end
end % end flgSkip

% TODO should I restart the parallel pool
BH_ctf_Refine2(varargin{1},varargin{2});


end % end of ctf estimate function

function [ bg, bandpass, rV ] = prepare_spectrum( Hqz, highCutoff, freqVector, radialAvg, dualAxis)


if numel(dualAxis) ~= max(size(dualAxis))
  bandpass2d = dualAxis;
  dualAxis = 1;
end




paddedSize = size(Hqz, 1);
%     if (dualAxis)
%       rV =  Hqz(1+paddedSize/2,1+paddedSize/2:end);
%     else
rV =  Hqz(1,1:paddedSize/2);
%     end


% smooth the spectrum a little to fit the zeros, and more to fit the max
% values
rVmin = convn(rV,[0.0180,0.0891,0.2327,0.3204,0.2327,0.0891,0.0180],'same');

[~,firstAbsMax] = max(abs(rV));
knots = [];
n = 1;
phs = -1;
while n < length(rV)
  
  % should start negative
  
  
  if phs < 0
    kn = find(rVmin(n:end) >= 0, 1, 'first');
    if isempty(kn)
      break
    end
    phs = 1;
    n = n + kn + 2 ;
    knots = [knots; n-3];
  else
    kn = find(rVmin(n:end) <= 0, 1, 'first');
    if isempty(kn)
      break
    end
    phs = -1;
    n = n + kn + 2;
    knots = [knots; n-3];
  end
  
end


% Take 15% of the first peak
halfPastFirstMax = knots(1) - floor(0.15*(knots(1) - firstAbsMax));
%     halfPastFirstMax = knots(2);


if (dualAxis)
  
  
  nCone = 1;
  
  for iCone = 2.5:5:90-2.5
    rImg = BH_resample2d(radialAvg,[0,0,iCone],[0,0,0],'Bah','GPU','inv',1,[paddedSize,paddedSize]);
    rAvg = abs(gather(rImg(1+paddedSize/2,1+paddedSize/2:end)));
    
    
    try
      bg{nCone} = fit(gather(freqVector(knots))', ...
        double(gather(rAvg(knots))'), ...
        'smoothingSpline');
    catch
      bg{nCone} = fit(gather(freqVector(knots(1:end-1))'), ...
        double(gather(rAvg(knots(1:end-1)))'), ...
        'smoothingSpline');
    end
    nCone = nCone +1;
    
  end
  
else
  
  if length(knots) >= 2
    try
      bg = fit(gather(freqVector(knots))', ...
        double(gather(abs(radialAvg(knots)))'), ...
        'smoothingSpline');
      
    catch
      bg = fit(gather(freqVector(knots(1:end-1)))', ...
        double(gather(abs(radialAvg(knots(1:end-1))))'), ...
        'smoothingSpline');
      
    end
    
    
  else
    fprintf('fewer than 2 minima detected\n')
  end
end

bandpass = false(1,length(freqVector));
bandpass(halfPastFirstMax:find(freqVector > highCutoff, 1, 'first')) = true;

if (dualAxis)
  avg = zeros(size(freqVector))';
  for iFit = 1:nCone-1
    avg = avg + bg{iFit}(freqVector);
  end
  
  avg = avg ./ (nCone-1);
  
  avgFit = fit(gather(freqVector)',avg,'cubicSpline');
  bg = avgFit;
end


if (dualAxis)
  
  freqLow = freqVector(find(bandpass,1,'first'));
  freqTop = freqVector(find(bandpass,1,'last'));
  bandpass = (bandpass2d > freqLow & bandpass2d < freqTop);
  
end

end

function [ iCCC ] = calc_CCC( freqVector, bg, bandpass, radialAvg, rV, cccScale)


if (numel(radialAvg) == max(size(radialAvg)))
  % One dimensional case
  bgSubPS = (abs(radialAvg) - bg(freqVector)').*bandpass;
  %bgSubPS = (bgSubPS - mean(bgSubPS(bandpass))) .* bandpass;
  if cccScale == 0.5
    % amplify the contribution of higher frequency information that may
    % otherwise be overwhelmed - particularly in CCD images. Only option is to
    % take the sqrt prior to normalizing. Provides some balance without really
    % risking amplifying too much noise.
    bgSubPS = (bgSubPS - min(bgSubPS(:)) + 1).^0.5;
  elseif cccScale ~= 1
    error('cccScale must be 1 or 0.5')
  end
  
  bgSubPS = bgSubPS ./ max(bgSubPS(:));
  % % % % %     bgSubPS = bgSubPS - mean(bgSubPS(:));
else
  if isnumeric(bg)
    % Use existing (otherwise bg is a cfit object)
    bgSubPS = bg;
  else
    bgSubPS = abs(radialAvg) - reshape(bg(freqVector),size(radialAvg));
    
    bgSubPS = bgSubPS .* bandpass;
    
    S = std2(bgSubPS(bandpass));
    bgSubPS(bgSubPS > S*2.5) = bgSubPS(bgSubPS>S*2.5).* ...
      (rand(gather(sum(bgSubPS(:)>S*2.5)),1)+0.5)./2 ;
    
    %     bgSubPS = bgSubPS - BH_movingAverage(bgSubPS,[8,8]);
    %     bgSubPS = bgSubPS ./ BH_movingRMS(bgSubPS,[8,8]);
    
    % % % % %     bgSubPS = bgSubPS ./ max(abs(bgSubPS(bandpass)));
    % % % % %     bgSubPS = (bgSubPS - mean(bgSubPS(bandpass))).*bandpass;
    bgSubPS = bgSubPS ./ max(abs(bgSubPS(bandpass))).*bandpass;
  end
end




if gather(rV(1)) == -9999
  % return the bgSubPS to save
  iCCC = bgSubPS;
else
  ctfSQ = abs(rV).*bandpass;
  % % % % %     ctfSQ = (ctfSQ- mean(ctfSQ(bandpass)).*bandpass);
  ctfSQ = ctfSQ ./ max(abs(ctfSQ(:)));
  iCCC = sum(sum((bgSubPS .* ctfSQ))) ./ ...
    ( numel(ctfSQ(bandpass)).*...
    std2(ctfSQ(bandpass)).*...
    std2(bgSubPS(bandpass)) );
  
end

end

function [ diagnosticIMG ] = make_diagnosticIMG( Hqz, pixel_size_si, bandpass, bg, IMG)


iImg = 1;
Hqz = fftshift(Hqz);
paddedSize = size(Hqz,1);
[radialGrid,~,~,~,~,~] = ...
  BH_multi_gridCoordinates(size(Hqz),'Cartesian',...
  'GPU',{'none'},1,0,1);

radialGrid = radialGrid ./ pixel_size_si;
lowCut = radialGrid(1, find(bandpass , 1,'first'));
highCut= radialGrid(1, find(bandpass , 1,'last'));

radialGrid = fftshift(radialGrid);
bandpass2d = (radialGrid < highCut & radialGrid > lowCut);
bgSubPS2d = (abs(IMG{iImg}) - reshape(bg(radialGrid),size(Hqz))).*bandpass2d;

diagnosticIMG = zeros(size(Hqz));

Hqz = abs(Hqz).* bandpass2d;
Hqz = 1.0.*Hqz ./ max(Hqz(bandpass2d)).*bandpass2d;
diagnosticIMG(1:(paddedSize/2),:,1) = gather(Hqz(1:(paddedSize/2),:));


bgSubPS2d = bgSubPS2d ./ max(bgSubPS2d(:));
diagnosticIMG((paddedSize/2)+1:(paddedSize/2)*2,:,1) = gather(bgSubPS2d((paddedSize/2)+1:(paddedSize/2)*2,:));

end

%%%%%%%%%%

function [ croppedIMG,pixelOUT ] = cropIMG(IMG,pixelIN)

maxRes = 3.3;
targetNyquist = (0.45*maxRes)/0.5;
[d1,d2] = size(IMG);

[radialGrid] = BH_multi_gridCoordinates( [d1,d2,1],'Cartesian','GPU', ...
  {'none'},1,0,1);

% For now, pixelIN only applies to pre-fixedpattern noise removal, while
% pixelOUT is the desired final cropping.
radialGrid = radialGrid./pixelIN;
rVX = radialGrid(1:ceil((d1+1)/2));
rVY = radialGrid(1:ceil((d2+1)/2));
clear radialGrid

cutX = find(rVX >= 1/targetNyquist, 1,'first');
cutY = find(rVY >= 1/targetNyquist, 1,'first');

% The actual nyquist will deviate from the target since it is trunctated
% to some pixel value
newNyquist = 1/rVX(cutX);
pixelOUT = gather(0.5*newNyquist);




% Prior approach: fftshift, trim (with taper and bandpass), ifftshift
% Trial approach: shift and trim with logical fftMask, ifftshift
fftMask = BH_fftShift([cutX,cutY],[d1,d2],1);
ifftMask = BH_fftShift(0,-2.*[cutX,cutY],1);

croppedIMG = fftn(gpuArray(IMG));
croppedIMG = croppedIMG(fftMask);

croppedIMG = real(ifftn(croppedIMG(ifftMask)));

clear fftMask ifftMask



end

function [ tltOrder ] = calc_dose_scheme(emc,rawTLT,anglesSkipped,PHASE_PLATE_SHIFT)

flgCosineDose = emc.('oneOverCosineDose');
startingAngle = emc.('startingAngle');
startingDirection = emc.('startingDirection');
doseSymmetricIncrement = emc.('doseSymmetricIncrement');
doseAtMinTilt = emc.('doseAtMinTilt');
nPrjs = length(rawTLT);
tltOrder = zeros(nPrjs,5);

nAngle = 2;

if (doseSymmetricIncrement < 0)
  % For doseSymmetricIncrement = 2, 3 deg
  % 0, 3, -3, -6, 6, 9 ...
  doseSymmetricIncrement = abs(doseSymmetricIncrement);
  flgFirstTilt = 0;
else
  % For doseSymmetricIncrement = 2, 3 deg
  % 0, 3, 6, -3, -6, 9 ...
  flgFirstTilt=1;
end

if any(PHASE_PLATE_SHIFT)
  if diff(PHASE_PLATE_SHIFT) < 1e-3
    PHASE_PLATE_SHIFT(2) = PHASE_PLATE_SHIFT(1) + 1e-3;
  end
  extraPhaseShift = [PHASE_PLATE_SHIFT(1):(PHASE_PLATE_SHIFT(2) - PHASE_PLATE_SHIFT(1))./nPrjs:PHASE_PLATE_SHIFT(2)];
else
  extraPhaseShift = zeros(nPrjs,1);
end

totalDose = doseAtMinTilt;

if (anglesSkipped)
  % Get the actual angles from the index
  anglesToSkip = rawTLT(anglesSkipped);
  anglesToKeep = ~ismember(1:nPrjs,anglesSkipped);
else
  anglesToSkip = [];
  anglesToKeep = true(nPrjs,1);
end


% We always start from the first tilt.
[~,firstTilt] = min(abs(rawTLT-startingAngle));
tltOrder(firstTilt,:) = [firstTilt,rawTLT(firstTilt),totalDose,extraPhaseShift(1),0];
% Remove this angle to get those remaining
tmpTLT = rawTLT([1:firstTilt-1,firstTilt+1:end]);
largerAngles = tmpTLT(tmpTLT-startingAngle > 0);
smallerAngles= tmpTLT(tmpTLT-startingAngle < 0);

% Now split into thos that are larger or smaller than the min tilt
if ( startingAngle >= 0 )
  largerAngles = sort(largerAngles,'ascend');
  smallerAngles = sort(smallerAngles,'descend');
else
  largerAngles = sort(largerAngles,'ascend');
  smallerAngles = sort(smallerAngles,'descend');
end

clear tmpTLT

if ( doseSymmetricIncrement )
  % It is assumed that blocks of this many tilts are collected NOT
  % including the first tilt. If the original dose symmetric scheme is
  % requested (negative Increment) then the first tilt IS included, and
  % so we need to subtract one from the counter.
  switchAfterNTilts = doseSymmetricIncrement - (1-flgFirstTilt)
  flgFirstTilt=0;
else
  if strcmpi(startingDirection,'pos')
    switchAfterNTilts = length(largerAngles);
  elseif strcmpi(startingDirection,'neg')
    switchAfterNTilts = length(smallerAngles);
  else
    error('flgDose symmetric is 0 and starting direction must be pos or neg');
  end
end


%    while (~isempty(largerAngles) || ~isempty(smallerAngles)) && nAngle <= nPrjs
for iPrj = 1:nPrjs
  if iPrj == firstTilt
    continue;
  end
  
  if strcmpi(startingDirection,'pos')
    try
      nextTilt = largerAngles(1);
      if length(largerAngles) > 1
        largerAngles = largerAngles(2:end);
      else
        largerAngles = [];
      end
    catch
      nextTilt = smallerAngles(1);
      if length(smallerAngles) > 1
        smallerAngles = smallerAngles(2:end);
      else
        smallerAngles = [];
      end
    end
    
  else
    try
      nextTilt = smallerAngles(1);
      if length(smallerAngles) > 1
        smallerAngles = smallerAngles(2:end);
      else
        smallerAngles = [];
      end
    catch
      nextTilt = largerAngles(1);
      if length(largerAngles) > 1
        largerAngles = largerAngles(2:end);
      else
        largerAngles = [];
      end
    end
  end
  
  [~,iTilt] = min(abs(rawTLT-nextTilt));
  
  switchAfterNTilts = switchAfterNTilts -1;
  
  if (switchAfterNTilts == 0)
    if strcmpi(startingDirection,'pos')
      startingDirection = 'neg';
    elseif strcmpi(startingDirection,'neg')
      startingDirection = 'pos';
    end
    if ( doseSymmetricIncrement )
      switchAfterNTilts = doseSymmetricIncrement;
    end
  end
  
  if (flgCosineDose)
    totalDose = totalDose + (1/cosd(rawTLT(iTilt)))*doseAtMinTilt;
  else
    totalDose = totalDose + doseAtMinTilt;
  end
  
  
  % The dose is incremented but don't add to the list.
  
  if ~ismember(rawTLT(iTilt),anglesToSkip)
    tltOrder(iTilt,:) = [iTilt,rawTLT(iTilt),totalDose,extraPhaseShift(nAngle),-1];
  else
    tltOrder(iTilt,:) = [iTilt,rawTLT(iTilt),-1,-1,-1];
  end
  nAngle = nAngle + 1;
  
end



tltOrder = tltOrder(anglesToKeep,:);
tltOrder(:,5) = tltOrder(:,1);
tltOrder(:,1) = 1:size(tltOrder,1);
tltOrder

size(tltOrder)



% This will be a naive run through that works only if the angles are in
% order.

end

