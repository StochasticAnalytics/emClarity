
function [  ] = BH_ctf_Find4(varargin)

!mkdir -p aliStacks

nGPUs=1;

  % PARAMETER_FILE, STACK_BASENAME, gpuIDX
  % Assumes that all required files are in ./fixedStacks
  PARAMETER_FILE = varargin{1};
  STACK_BASENAME = varargin{2};
  anglesSkipped = 0;
  gpuToUse = -1;
  
  findExtraPhaseShift = 0

if length(varargin) > 2

  anglesSkipped = str2num(varargin{3});

end
if length(varargin) > 3
 
  gpuToUse = str2num(varargin{4});  

end
if length(varargin) > 4

  findExtraPhaseShift = 1

end

  stackNameIN = sprintf('fixedStacks/%s.fixed', STACK_BASENAME);
  tltFile = sprintf('fixedStacks/%s.tlt', STACK_BASENAME);
  stackNameOUT = sprintf('%s_ali1', STACK_BASENAME);
  mapBack = sprintf('fixedStacks/%s', STACK_BASENAME);
  collectionORDER = sprintf('fixedStacks/%s.order',STACK_BASENAME);
  

pBH = BH_parseParameterFile(PARAMETER_FILE);

gpuIDX = BH_multi_checkGPU(gpuToUse);

fprintf('Removing tilts, from 1:\n %2.1f,\n',anglesSkipped);
% 
% gDev = gpuDevice(gpuIDX);

flgResume = 0;
flgSkip = 0;
% slightly dampen lower resolution information that may overwhelm the CCC calc,
% but don't risk too much noise amplification. 1 = fit just the amplitude (not
% PS) 0.5 = take sqrt prior to normalizing.
cccScale = 1;

% Assumes a directory with mapBack data on the same level as the directory
% holding the rawStacks, and that the basename of the stack is also the basename
% for the mapBack transformations.
flgReScale = 0;
flgMapBack = str2num(mapBack);
if isempty(flgMapBack)
  % The string specifies the basename for mapback files .xf .mag .tlt
  % and so is empty after str2num
  flgReScale = true;
  mapBackPrfx= mapBack;
elseif (flgMapBack)
  error('mapBack should be false (a 0) or a string with mapBack/basename.[xf,mag,tlt]\n')
end


try
  rawTLT = load(tltFile);
catch
  error('The tilt angle file %s was not found.\n', tltFile);
end

flgStandardOrdeDoCalc = 1;
try
  flgCosineDose = pBH.('oneOverCosineDose');
  startingAngle = pBH.('startingAngle');
  startingDirection = pBH.('startingDirection');
  doseSymmetricIncrement = pBH.('doseSymmetricIncrement');
  doseAtMinTilt = pBH.('doseAtMinTilt');

  flgOldDose = 0;
  tltOrder = calc_dose_scheme(pBH,rawTLT,anglesSkipped);
  
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



PIXEL_SIZE = pBH.('PIXEL_SIZE');
Cs = pBH.('Cs');
VOLTAGE = pBH.('VOLTAGE');
AMPCONT = pBH.('AMPCONT')
SuperResolution = pBH.('SuperResolution');

if (SuperResolution)
  if SuperResolution == 1
    % Standard scenario crop to physical nyquist
    scalePixelsBy = 2;    
  elseif SuperResolution > 10^10*PIXEL_SIZE
    % Crop to the given pixels size
    error('Scaling to arbitrary pixel size is not working\n');
    % Need to factor in the trunctation to integer pixel size.
%     scalePixelsBy = SuperResolution/(10^10*PIXEL_SIZE);
  else
    error('SuperResolution must be 0 (off) 1 (crop to physical Nyquist) or a pixel Size larger than current\n');
  end
  PIXEL_SIZE = scalePixelsBy.* PIXEL_SIZE;
else
  scalePixelsBy = 1;
end

if 10^10*PIXEL_SIZE < 1.2
  fprintf('PixelSize is less than 1.2 Ang so we have to use the cpu\n');
  useGPU = 0;
  METHOD = 'cpu';
else
  useGPU = 1;
  METHOD = 'GPU';
end

% Sanity check
if (PIXEL_SIZE > 20e-10 || PIXEL_SIZE < 0)
  error('pixel size should be [0,20e-10]');
elseif (Cs > 10e-3 || Cs < 1e-3)
  error('Cs should be[1e-3,10e-3]');
elseif(VOLTAGE > 1000e3 || VOLTAGE < 20e3)
  error ('VOLTAGE should be [20e3,1000e3]');
elseif (AMPCONT < 0.025 || AMPCONT > 0.25)
  error('AMPCONT should be [0.025,0.25]');
else
  WAVELENGTH = 10^-12*1226.39/sqrt(VOLTAGE + 0.97845*10^-6*VOLTAGE^2) ;
end
 
CUM_e_DOSE = pBH.('CUM_e_DOSE');
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

eraseRadius = ceil(1.2.*(pBH.('beadDiameter')./PIXEL_SIZE.*0.5));
flgImodErase = 0

  
  
% Assuming that the first CTF zero is always less than this value 
FIXED_FIRSTZERO =  PIXEL_SIZE / (70*10^-10) ;
highCutoff = PIXEL_SIZE/pBH.('defCutOff');
% I still use the def for underfocus < 0 as this places the origin at the
% focal plan in the microscope rather than on the specimen. Which makes
% more sense to me.

defEST = pBH.('defEstimate').*10^10

tiltRange = [-1];

backGroundBuffer = 0.9985;

try
  deltaZTolerance = pBH.('deltaZTolerance');
catch
  deltaZTolerance = 50e-9;
end

try 
  zShift = abs(pBH.('zShift'));
catch
  zShift = 150e-9;
end

if abs(zShift) > 100e-7
  error('make sure your zShift values are of reasonable amounts (50-200nm)');
end

% Starting at +/- 100nm
deltaZTolerance = deltaZTolerance / PIXEL_SIZE;
% Use to check for proper gradient.
zShift = zShift / PIXEL_SIZE;

% Tile size & overlap
try 
  tileSize = pBH.('ctfTileSize');
catch
  tileSize = floor(680e-10 / PIXEL_SIZE);
end
tileOverlap = 12;

tileSize = tileSize + mod(tileSize,2);
fprintf('Using a tile size of %d\n',tileSize);
overlap = floor(tileSize ./ tileOverlap);

% Size to padTile to should be even, large, and preferably a power of 2
try
  paddedSize = pBH.('paddedSize');
catch
  paddedSize = 512;
end
padVAL = BH_multi_padVal([tileSize,tileSize], [paddedSize,paddedSize]);



if exist(stackNameIN, 'file')
  
  if size(tltOrder,1) ~= length(rawTLT)
    error('The length of the collectionOrder %d and tilt Geometry %d are different\n', ...
                                        length(tltOrder),length(tltInfo));
  end 
  
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
    included = (tltOrder(:,3) ~= -9999);
    nSkipped = (sum(~included)); % used to adjust header
    TLT = zeros(sum(included),23);
    TLT(:,1) = 1:sum(included);
    TLT(:,23) = tltOrder(included,1); % This was only the appropriate tilts are read in.
    TLT(:,4) = tltOrder(included,2); clear rawTLT
    TLT(:,11) = tltOrder(included,3);
    TLT(:,14) = 1;
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

% if ~(flgResume)



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
  
if (SuperResolution)
  halfMask = fftshift(BH_bandpass3d(1.*[d1,d2,1],0,0,4,METHOD,1));
end

% Copy with column for defocus = input to CTF correct
% saved as <filename>_ctf.tlt


flgReOrderMapBack = 0;

TLT(:,2:3) = repmat([0.00,0.00],size(TLT,1),1);
TLT(:,5:10) = repmat([0,90.0,1.0,0.0,0.0,1.0],size(TLT,1),1);
% Defocus will go at 15 - 12 and 13 currently unused.
TLT(:,16:19) = repmat([PIXEL_SIZE,Cs,WAVELENGTH,AMPCONT],size(TLT,1),1);
TLT(:,20:22) = repmat([d1,d2,d3],size(TLT,1),1);


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



fprintf('Combining tranformations\n\n');
% Load in the mapBack alignment
mbEST = load(sprintf('%s.xf',mapBackPrfx));
mbTLT = load(sprintf('%s.tlt',mapBackPrfx));

outputStackName = sprintf('aliStacks/%s%s',stackNameOUT,extension)
if exist(sprintf('%s.erase',mapBackPrfx),'file')
  flgEraseBeads = 1;
else
  flgEraseBeads = 0;
  fprintf('\nDid not find the gold bead file (%s) for erasing, will skip\n\n',sprintf('%s.erase',mapBackPrfx));
end

% % % % % 

if (SuperResolution)
  % Forcing output to odd size.
  sizeCropped = floor([d1,d2,d3]./2)-(1-mod(floor([d1,d2,d3]./2),2));
else
  sizeCropped = [d1,d2,d3]-(1-mod([d1,d2,d3],2));
end
sizeCropped(3) = d3; 

STACK = zeros(sizeCropped,'single');



if (flgReOrderMapBack)
  TLT = sortrows(TLT,1);
end

% if any([d1,d2] > 4096)
%   shiftMETHOD = 'cpu';
%   fprintf('transforming on cpu b/c > 4096\n')
% else
  shiftMETHOD = 'GPU';
% end

for i = 1:d3
  fprintf('Transforming prj %d in fourier space oversampled by 2x physical Nyquist\n',i);




  % Stored in row order as output by imod, st transpose is needed. Inversion
  % of the xform is handled in resample2d.
  origXF = [1,0;0,1];
   newXF = reshape(mbEST(i,1:4),2,2)';


  dXYZ  = [(newXF*TLT(i,2:3)')' + mbEST(i,5:6) , 0];
  TLT(i,2:3) = dXYZ(1:2);
  dXYZ = dXYZ ./ scalePixelsBy;

  combinedXF = reshape((newXF*origXF)',1,4);
  TLT(i,7:10) = combinedXF;


  osX = 1-mod(d1,2); osY = 1-mod(d2,2);

  % Pad the projection prior to xforming in Fourier space.
  if (SuperResolution)

     iProjection = single(getVolume(iMrcObj,-1,-1,TLT(i,23),'keep'));

     % Remove any very large outliers
     largeOutliersMean= mean(iProjection(:));
     largeOutliersSTD = std(iProjection(:));
     largeOutliersIDX = (iProjection < largeOutliersMean - 10*largeOutliersSTD | ...
                         iProjection > largeOutliersMean + 10*largeOutliersSTD);
     iProjection(largeOutliersIDX) = (3*largeOutliersSTD).*randn([gather(sum(largeOutliersIDX(:))),1],'single','gpuArray');
                 
    % Information beyond the physical nyquist should be removed to limit
    % aliasing of noise prior tto interpolation.
    iProjection = BH_padZeros3d(iProjection,[0,0],[0,0],shiftMETHOD,'singleTaper',largeOutliersMean);
    trimVal = BH_multi_padVal(1.*size(iProjection),sizeCropped(1:2));

    iProjection = real(ifftn(ifftshift(...
                               BH_padZeros3d(halfMask.*fftshift(...
                                             fftn(iProjection)), ...
                                             trimVal(1,:),trimVal(2,:),...
                                             shiftMETHOD,'single'))));  




     sizeODD = size(iProjection)-[osX,osY];
  else
     sizeODD = [d1,d2]-[osX,osY];

    % If it is even sized, shift up one pixel so that the origin is in the middle
    % of the odd output here we can just read it in this way, unlike super res. 
    
      iProjection = ...
                 single(getVolume(iMrcObj,[1+osX,d1],[1+osY,d2],TLT(i,23),'keep'));
     largeOutliersMean= mean(iProjection(:));
     largeOutliersSTD = std(iProjection(:));
     largeOutliersIDX = (iProjection < largeOutliersMean - 6*largeOutliersSTD | ...
                         iProjection > largeOutliersMean + 6*largeOutliersSTD);
     iProjection(largeOutliersIDX) = (3*largeOutliersSTD).*randn([gather(sum(largeOutliersIDX(:))),1],'single');

             

  end

     % Because the rotation/scaling and translation are done separately,
     % we must use a square transform; otherwise, a rotation angle dependent
     % anisotropic distortion (like mag distortion) is introduced.                                         
    sizeSQ = floor(2.*[1,1].*max(sizeODD));
    padVal  = BH_multi_padVal(sizeODD,sizeSQ);
     trimVal = BH_multi_padVal(sizeSQ,sizeCropped(1:2));
      % Only need to do this for the first pass
  if ( i == 1 )
    
    if strcmpi(shiftMETHOD, 'GPU')
      [ fftMask ] = BH_fftShift(0,sizeSQ,1); 
      [ ifftMask ] = BH_fftShift(0,-1.*sizeSQ,1); 
    else
      [ fftMask ] = BH_fftShift(0,sizeSQ,0); 
      [ ifftMask ] = BH_fftShift(0,-1.*sizeSQ,0); 
    end
    % Calculate grids in reciprocal pixels including 2pi for phase shifting
    [ dU, dV ] = BH_multi_gridCoordinates(sizeSQ,'Cartesian',shiftMETHOD, ...
                                                    {'none'},1,1,0);
    dU = dU .* (-2i*pi);
    dV = dV .* (-2i*pi);

    % Remake with potential new size, and adjusted cutoff
    halfMask2 = fftshift(BH_bandpass3d([sizeSQ,1],0,0,2,shiftMETHOD,1));
  end

  iProjection = iProjection - mean(iProjection(:));
  

  if ( SuperResolution )
    iProjection = BH_padZeros3d(iProjection(1+osX:end,1+osY:end), ...
                            padVal(1,:),padVal(2,:),shiftMETHOD,'singleTaper');
  else
    iProjection = BH_padZeros3d(iProjection,padVal(1,:),padVal(2,:), ...
                                                    shiftMETHOD,'singleTaper');
  end



  
  iProjection = fftn(iProjection(ifftMask));
  iProjection = iProjection(fftMask);
  

% % % % %    iProjection = fftshift(fftn(ifftshift(iProjection)));



 % Do the phase shift after rotating - need to invert the scaling since
 % we are in reciprocal space
 [imodMAG, imodStretch, imodSkewAngle, imodRot] = ...
                                       BH_decomposeIMODxf(combinedXF);
 % Assuming stretch and skew are not fit, leave defined for possible
 % later consideration.

 ddXddY = [0,0];%[combinedXF(1:2),0;combinedXF(3:4),0;0,0,imodMAG] * [-1.5,0.5,0]';
 combinedInverted = BH_defineMatrix([imodRot,0,0],'Bah','forward').*(1/imodMAG);

 combinedInverted = combinedInverted([1,2,4,5]);


 iProjectionR = (BH_resample2d(real(iProjection).*halfMask2, combinedInverted,[0,0,0],...
                           'Bah',shiftMETHOD,'forward',1.0,size(iProjection))); 
 iProjectionI = (BH_resample2d(imag(iProjection).*halfMask2, combinedInverted,[0,0,0],...
                           'Bah',shiftMETHOD,'forward',1.0,size(iProjection)));  


                         
 iProjection = imodMAG.^2.*exp(dU.*(dXYZ(1)+ddXddY(1)) + dV.*(dXYZ(2)+ddXddY(2))) .* ...
               complex(iProjectionR,iProjectionI);


 iProjectionR = []; iProjectionI = [];

 iProjection = real(ifftn(iProjection(ifftMask)));
 iProjection = iProjection(fftMask);
 

% % % % %    iProjection = real(fftshift(ifftn(ifftshift(iProjection))));
 STACK(:,:,i)  = gather(real(BH_padZeros3d(iProjection, ...
                                       trimVal(1,:),trimVal(2,:),...
                                       shiftMETHOD,'single')));


end 

if ( flgEraseBeads )
    system(sprintf('imodtrans -i fixedStacks/%s.fixed fixedStacks/%s.erase fixedStacks/%s.tmpMod',fileName,fileName,fileName));
    system(sprintf('model2point -float fixedStacks/%s.tmpMod  fixedStacks/%s.erase2',fileName,fileName));
    system(sprintf('rm fixedStacks/%s.tmpMod',fileName));
    %beadList = importdata(sprintf('fixedStacks/%s.erase2',fileName));
    beadList = load(sprintf('fixedStacks/%s.erase2',fileName));
   
    beadList(:,1:2) = beadList(:,1:2) ./ scalePixelsBy;
    STACK = BH_eraseBeads(STACK,eraseRadius, beadList);
end 

[ STACK ] = BH_multi_loadAndMaskStack(STACK,TLT,'',100,PIXEL_SIZE*10^10);

SAVE_IMG(MRCImage(STACK),outputStackName,iPixelHeader,iOriginHeader);

clear dU dV fftMask ifftMask halfMask2 iProjection largeOutliersIDX
% % % % % 


if ~(flgSkip)
  

[d1,d2,d3] = size(STACK);
if (PIXEL_SIZE*10^10 < 0)
  
  flgCrop = 1;
  [croppedIMG,pixelOUT] = cropIMG(STACK(:,:,1),PIXEL_SIZE*10^10);
  [d1C,d2C] = size(croppedIMG); clear croppedIMG
  tltForExp = TLT;
  tltForExp(:,16) = pixelOUT*10^-10;
  pixelOUT
  % Redefining things down hear is a stupid thing to do. Fix this if you
  % keep the optino for cropping.
  FIXED_FIRSTZERO =  pixelOUT / 70 ;
  highCutoff = (pixelOUT*10^-10)/pBH.('defCutOff');
  
else
  
  flgCrop = 0;
  pixelOUT = PIXEL_SIZE*10^10;
  d1C = d1;
  d2C = d2;
  tltForExp = TLT;
  
end

d3 = size(STACK,3)
% Check for extra large (8k) data which will be too big for the gpu.
% Should set this up to be a hybrid where  each slice is on GPU but 
% But then pull to the cpu and store there in stack.
if d1C > 4096 || d2C > 4096 || d3 > 40 
   prjMaskMethod = 'GPU'
 else
  prjMaskMethod = 'GPU'
end



tic

nStrips = 13; % Must be odd
maxWorkers = 7*nGPUs; % TODO get from paramFile
nWorkers = min(nStrips,maxWorkers);
whos
delete(gcp('nocreate'))
EMC_parpool(nWorkers);
parfor iWorker = 1:nWorkers
  gpuDevice(gpuIDX);
end

% For visualization
bandNyq = BH_bandpass3d([paddedSize.*[1,1],1],0,0,0,'GPU','nyquistHigh');
bandLow = fftshift(BH_bandpass3d([paddedSize.*[1,1],1],0,0,2.15,'GPU',1));
tiltBaseName = {[pathName '/ctf/' stackNameOUT],stackNameOUT};

  
% First get ranges for defocus, check quality and handedness.
[~,tiltLow] = min(abs(TLT(:,4)));
[~,tiltpos20 ] = min(abs(TLT(:,4)-30)); 
[~,tiltneg20 ] = min(abs(TLT(:,4)+30)); 
% Make sure we use the tilt with less exposure
if (TLT(tiltpos20,11) < TLT(tiltneg20,11))
  fprintf('Using +20 degree tilt to check handedness\n');
  tilt20 = tiltpos20;
  tiltDirection = 1;
else
  fprintf('Using -20 degree tilt to check handedness\n');
  tilt20 = tiltneg20;
  tiltDirection = -1;
end

system(sprintf('mkdir -p %s',tiltBaseName{1}));
lowResFit = defEST*0.001
lowResFit = min(50, max(16, lowResFit))
% Adjusted again after initial check

highResFit = max(2.5*PIXEL_SIZE*10^10, lowResFit/5)

% For zero-tilt
tiltInfo = [PIXEL_SIZE*10^10,TLT(tiltLow,4),0.6*defEST,1.4*defEST,lowResFit,highResFit, ...
             300.0,2.70,0.07,paddedSize,findExtraPhaseShift];
           
[ fitVals, ~, ~, ~ ] = run_ctffind(STACK(:,:,TLT(tiltLow,1)), true, tiltBaseName,tiltInfo, 1);
highestResToFit = min(fitVals(7),10); % TODO set 10 somehwere else
likelyDefocus   = mean(fitVals(2:3));
likelyAstig     = (fitVals(2) - fitVals(3))/2;
likelyAngle     = fitVals(4);
fprintf('For zero angle tilt:\nhighest res fit %2.2f\n',highestResToFit);
fprintf('nominal def %6.0f \nAstig %6.0f Ang at %2.2f deg\n', ...
        likelyDefocus, likelyAstig, likelyAngle);
astigRestraint = sqrt(abs(likelyAstig))*10;
fprintf('Using a restraint on astigmatism of %6.2f\n',astigRestraint);

fprintf('Adjusting lowResFit from lowResFit %2.2f to',lowResFit);
lowResFit = likelyDefocus*0.001;
lowResFit = min(50, max(16, lowResFit));
fprintf(' %2.2f\n',lowResFit);


tiltInfo        = repmat(tiltInfo,nStrips,1);
tiltInfo(:,2)   = repmat(TLT(tilt20,4),nStrips,1);
tiltInfo(:,3:4) = repmat([0.85,1.15].*likelyDefocus,nStrips,1); % TODO add range check
tiltInfo(:,5:6)   = repmat([lowResFit,highestResToFit],nStrips,1);


     

% Chech handedness first
searchVect = 1:d3;
searchVect = searchVect(searchVect ~= tilt20);
searchVect = [tilt20, searchVect];
checkHand = true;

for iPrj = searchVect
  

  tiltInfo(:,2)   = repmat(TLT(iPrj,4),nStrips,1);

  tiltSign = sign(TLT(iPrj,4));
  
  iEvalMask = zeros(d1C,nStrips,'single','gpuArray');
  % Center the pixel coordinates
  iEvalMask(:,1) = BH_multi_gridCoordinates([d1C,1,1],'Cartesian','GPU',{'none'},0,1,0);
  iVectX_selected = zeros(nStrips,1);
  % Convert to the z-height in the projection
  iEvalMask(:,1) = iEvalMask(:,1).*(-1.*tand(TLT(iPrj,4)));
  

  % Get the max shift
  rMaxZ = abs(iEvalMask(1,1)) ./ (nStrips-2);
  fprintf('rMaxZ increment is %3.3f for iPrj %2.2f\n',rMaxZ,TLT(iPrj,4));
  
  for iCopy = 2:nStrips
    iEvalMask(:,iCopy) = iEvalMask(:,1);
  end
  

  
  for iMask = 1:nStrips
    
    % Change in Z-height (positive then closer to focust)

    iDefShift_Pix = -1*tiltSign*rMaxZ*(iMask-(ceil(nStrips/2)));
    iDefShift_Ang = gather(iDefShift_Pix .* tiltInfo(1,1));
    iVectX_selected(iMask) = -iDefShift_Ang;

    % Adjust the search range more narrowly about the expected defocus
    % ctffind uses underfocus > 0 so subtract the z-shift to move closer to focus
    tiltInfo(iMask,3:4) = [0.95,1.05].*(likelyDefocus - iDefShift_Ang);

    % End on 1 as others start as zeros

    iEvalMask(:,iMask) = iEvalMask(:,iMask) - iDefShift_Pix;

    iEvalMask(:,iMask) = (iEvalMask(:,iMask) > gpuArray(-deltaZTolerance/2) &...
                          iEvalMask(:,iMask) < gpuArray(deltaZTolerance/2));
  end
  

  
  par_tmpTile = cell(nWorkers,1);
  par_nT = cell(nWorkers,1);
  for iWorker = 1:nWorkers
    
    par_tmpTile{iWorker} = zeros([paddedSize.*[1,1],nStrips],'single','gpuArray');
    par_nT{iWorker} = zeros(nStrips,1,'single','gpuArray');

  end
  
%     iProjection = double(gpuArray(STACK(:,:,TLT(k,1))));
  if flgCrop
    [iProjection,~] = cropIMG(gpuArray(STACK(:,:,TLT(iPrj,1))),PIXEL_SIZE*10^10);
  else
    iProjection = (gpuArray(STACK(:,:,TLT(iPrj,1))));
  end
    
  iProjection = iProjection - ...
                       BH_movingAverage(iProjection,[tileSize,tileSize]./3);
                  
  iProjection = iProjection ./ ...
                       BH_movingRMS(iProjection,[tileSize,tileSize]./3);



 
  parVect = [1+tileSize/2:overlap:d1C-tileSize/2];
  parfor iWorker = 1:nWorkers
    for iW = iWorker:nWorkers:length(parVect)
      i = parVect(iW);
%       for i = 1+tileSize/2:overlap:d1C-tileSize/2
        iEvalVector = gather(find(iEvalMask(i,:) > 0));
        if iEvalVector
          for j = 1+tileSize/2:overlap:d2C-tileSize/2  

          %if evalMask(i,j,TLT(k,1))


             thisTile = abs(fftn( ...
                                     BH_padZeros3d(...
                                     (iProjection( ...
                                            i-tileSize/2+1:i+tileSize/2,...
                                            j-tileSize/2+1:j+tileSize/2)),...
                                                      padVAL(1,:),padVAL(2,:),...
                                                      'GPU','singleTaper')));

              for iTile = 1:length(iEvalVector)

               par_tmpTile{iWorker}(:,:,iEvalVector(iTile)) = par_tmpTile{iWorker}(:,:,iEvalVector(iTile)) + thisTile;
               par_nT{iWorker}(iEvalVector(iTile)) = par_nT{iWorker}(iEvalVector(iTile)) + 1;
              end

          end
        end
    end
  end
  
  tmpTile = zeros([paddedSize.*[1,1],nStrips],'single','gpuArray');
  nT = zeros(nStrips,1,'single','gpuArray');
  for iWorker = 1:nWorkers
    tmpTile = tmpTile + par_tmpTile{iWorker};
    nT = nT + par_nT{iWorker};
  end
  
  clear par_tmpTile par_nT
  
  fprintf('for iPrj %d, found tiles',iPrj);
  fprintf(' %d\n',int16(nT));
  

  for iStrip = 1:nStrips
%     if iStrip==1 ; figure, imshow3D(gather(tmpTile(:,:,iStrip)));end
    tmpTile(:,:,iStrip) = fftshift(tmpTile(:,:,iStrip).*bandNyq)./nT(iStrip);
%     if iStrip==1 ; figure, imshow3D(gather(tmpTile(:,:,iStrip)));end
    tmpTile(:,:,iStrip) = (tmpTile(:,:,iStrip) - BH_movingAverage(tmpTile(:,:,iStrip),6.*[1,1],2));
    tmpTile(:,:,iStrip) = bandLow.*(tmpTile(:,:,iStrip) ./ BH_movingRMS(tmpTile(:,:,iStrip),6.*[1,1],2));
%    if iStrip==1 ; figure, imshow3D(gather(tmpTile(:,:,iStrip))); figure, imshow3D(gather(BH_movingAverage(tmpTile(:,:,iStrip),0.1.*paddedSize.*[1,1],2)));end
%    return
  end

 

  [ fitVals, defFit, astigMagFit, astigAngFit ] = run_ctffind(tmpTile, false, tiltBaseName,tiltInfo, nStrips, astigRestraint, iVectX_selected); 
  
  TLT(iPrj,5:6)   = [astigMagFit.p1,astigMagFit.p2.*10^-10];
  TLT(iPrj,12:13) = [astigAngFit.p1,astigAngFit.p2.*(pi/180)];
  TLT(iPrj,15)    = [defFit.p2*10^-10];

  % % add a check on both the average score and its std
  % % consider refitting with a weighted avg of the best scores if std is
  % large
  % % consider tracking mean score tilt to tilt to watch for outliers and
  % replace with a neighbor.
  % consider fitting the zero-degree first to find the max res to use and
  % also the ballpark for range to fit.
  % it may also be better to further restrict the search range on tilts by
  % fitting each indivdual rather than a stack
  % a second lower tilt, but big enough to get the gradient would then be
  % useful to confirm such that the restrictions are relevant.

  if (checkHand)


    if ( defFit.p1 < 0)
      % TODO add a flag to override this or to check another tilt first
      fprintf('\n\nIt looks like your handedness is inverted!\n\n');
      defVect = -floor(nStrips_handCheck/2):floor(nStrips_handCheck/2);
      figure('visible','off'), plot(defVect,defFit(defVect),'b',defVect,mean(fitVals(:,2:3),2),'k');
      saveas(gcf,sprintf('%s/%shandednessError_1.pdf',tiltBaseName{:}));
      error('Check the plot for fixedStacks/ctf/handednessError.pdf');
    else
      fprintf('\nIt looks like your handedness is probably correct\n');
    end
    
    checkHand = false;
  end
  

  
end


    % number in stack, dx, dy, tilt angle, projection rotation, tilt azimuth, tilt
    % elevation, e1,e2,e3, dose number (order in tilt collection), offsetX, offsetY
    % scaleFactor, defocus, pixelSize, CS, Wavelength, Amplitude contrast
    fileID = fopen(sprintf('%s/ctf/%s_ctf.tlt',pathName,stackNameOUT), 'w');
    fprintf(fileID,['%d\t%08.2f\t%08.2f\t%07.3f\t%6e\t%6e\t%07.7f\t%07.7f\t',...
             '%07.7f\t%07.7f\t%5e\t%5e\t%5e\t%7e\t%5e\t%5e\t%5e\t%5e\t%5e\t',...
             '%d\t%d\t%d\t%8.2f\n'], TLT');
    fclose(fileID);


    
end % end flgSkip


delete(gcp('nocreate'));

end % end of flag resume (partially killed)





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

function [ tltOrder ] = calc_dose_scheme(pBH,rawTLT,anglesSkipped)

  flgCosineDose = pBH.('oneOverCosineDose');
  startingAngle = pBH.('startingAngle');
  startingDirection = pBH.('startingDirection');
  doseSymmetricIncrement = pBH.('doseSymmetricIncrement');
  doseAtMinTilt = pBH.('doseAtMinTilt');
  
  tltOrder = zeros(length(rawTLT),3);
  
  totalDose = doseAtMinTilt;
  nMax = length(rawTLT)+1;
  
  if (anglesSkipped)
    % Get the actual angles from the index
    anglesToSkip = rawTLT(anglesSkipped);
  else
    anglesToSkip = [];
  end


    % We always start from the first tilt.
   [~,iTilt] = min(abs(rawTLT-startingAngle));
   tltOrder(iTilt,:) = [iTilt,rawTLT(iTilt),totalDose];
   tmpTLT = rawTLT([1:iTilt-1,iTilt+1:end]);
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
     % counting the first tilt. 
     switchAfterNTilts = doseSymmetricIncrement-1;
   else
     if strcmpi(startingDirection,'pos')
       switchAfterNTilts = length(largerAngles);
     elseif strcmpi(startingDirection,'neg')
       switchAfterNTilts = length(smallerAngles);
     else
       error('flgDose symmetric is 0 and starting direction must be pos or neg');
     end
   end

   maxTries = length(largerAngles) + length(smallerAngles) +2;
   nAngle = 2;
   while ~isempty(largerAngles) || ~isempty(smallerAngles)
    
     if strcmpi(startingDirection,'pos')
       try
        nextTilt = largerAngles(1);
       catch
         nextTilt = smallerAngles(1);
       end
       largerAngles = largerAngles(2:end);
     else
       try
        nextTilt = smallerAngles(1);
       catch
         nextTilt = largerAngles(1);
       end
       smallerAngles = smallerAngles(2:end);
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
      tltOrder(iTilt,:) = [iTilt,rawTLT(iTilt),totalDose];
     else
       tltOrder(iTilt,:) = [iTilt,rawTLT(iTilt),-9999];
     end
      nAngle = nAngle + 1;
   
        
     if (maxTries == 0)
       error('max iterations in creating tilt order exceeded, breaking out');
     end
     maxTries = maxTries - 1;
   end
     

     
    tltOrder
    
    
  
  % This will be a naive run through that works only if the angles are in
  % order.
  
end
  
function [fitVals, defFit, astigMagFit, astigAngFit] = run_ctffind(imgStack, doFFT, tiltBaseName, tiltInfo, nStrips, varargin)

  defFit = '';
  astigMagFit = '';
  astigAngFit= '';
  fitVals = zeros(nStrips,7);
  
  pixelSize  = tiltInfo(:,1);
  tiltAngle  = tiltInfo(:,2);
  fit_defLow = tiltInfo(:,3);
  fit_defTop = tiltInfo(:,4);
  fit_resLow = tiltInfo(:,5);
  fit_resTop = tiltInfo(:,6);

  findExtraPhaseShift = tiltInfo(1,11);
  
  % TODO allow
  DEFOCUS_STEP = 100; 
  if nargin > 5
    ASTIGMATISM_RESTRAINT = varargin{1};
  else
    ASTIGMATISM_RESTRAINT = 250.0; % Penalty on Astig, lower = softer restraint, def 100
  end
  if nargin > 6
    defVect = varargin{2};
  else
    defVect  = [-floor(nStrips/2):floor(nStrips/2)]';
  end
  
  parfor iStrip = 1:nStrips
    
    
    SAVE_IMG(MRCImage(gather(imgStack(:,:,iStrip))),sprintf('%s/%s_%2.2f_%d.mrc',tiltBaseName{:},tiltAngle(iStrip),iStrip));

    tFOUT = fopen(sprintf('%s/%s_%2.2f_%d.sh',tiltBaseName{:},tiltAngle(iStrip),iStrip),'w');
    if (doFFT)
      fprintf(tFOUT,'#!/bin/bash\n\nctffind << eof\n');
    else
      fprintf(tFOUT,'#!/bin/bash\n\nctffind --amplitude-spectrum-input  << eof\n');
    end
    fprintf(tFOUT,'%s/%s_%2.2f_%d.mrc\n',tiltBaseName{:},tiltAngle(iStrip),iStrip); % Need a "no" after this line if multi-z stack
    fprintf(tFOUT,'%s/d_%s_%2.2f_%d.mrc\n%f\n',tiltBaseName{:},tiltAngle(iStrip),iStrip,pixelSize(iStrip));
    fprintf(tFOUT,'%f\n%f\n%f\n%f\n',tiltInfo(iStrip,7:10));
    fprintf(tFOUT,'%f\n%f\n%f\n%f\n',fit_resLow(iStrip),fit_resTop(iStrip)./cosd(tiltAngle(iStrip)),fit_defLow(iStrip),fit_defTop(iStrip));
    fprintf(tFOUT,'%f\n%s\n%s\n%s\n%2.2f\n%s\n\n',DEFOCUS_STEP,'no','no','yes',ASTIGMATISM_RESTRAINT);
    if (findExtraPhaseShift)
      minPhase = 0.0;
      maxPhase = pi/2;
      phaseStep= pi/200;
      fprintf(tFOUT,'%s\n%f\n%f\n%f\n%s\n%s\n','yes',minPhase,maxPhase,phaseStep,'yes','yes');
    else
      fprintf(tFOUT,'%s\n%s\n%s\n','no','yes','yes');
    end
    fprintf(tFOUT,'eof\n');
    fclose(tFOUT);
    system(sprintf('chmod a=wrx %s/%s_%2.2f_%d.sh',tiltBaseName{:},tiltAngle(iStrip),iStrip));
    system(sprintf('%s/%s_%2.2f_%d.sh > /dev/null',tiltBaseName{:},tiltAngle(iStrip),iStrip))
%     system(sprintf('%s/%s_%2.2f_%d.sh ',tiltBaseName{:},tiltAngle(iStrip),iStrip))


    system(sprintf('awk ''!/^#/{ print $0 }'' %s/d_%s_%2.2f_%d.txt > %s/d_%s_%2.2f_%d_noHeader.txt', ...
                   tiltBaseName{:},tiltAngle(iStrip),iStrip, ...
                   tiltBaseName{:},tiltAngle(iStrip),iStrip));
                   
  end
          
  for iStrip = 1:nStrips
    fitVals(iStrip,1:7) = load(sprintf('%s/d_%s_%2.2f_%d_noHeader.txt',tiltBaseName{:},tiltAngle(iStrip),iStrip));
  end
  % Use the cross-correlation fits to weight the fitting
  flgUseWeights = 0;
  % Least squares robust method, bisquare, LAR, '' (empty to skip)
  flgRobust = 'bisquare';
  % Check to see if there is variation in astigmatism
  flgAstig = 1;
  if (nStrips > 1)
    
    defMean  = mean(fitVals(:,2:3),2);
    astigMag = (fitVals(:,2)-fitVals(:,3))./2;
    astigAng = fitVals(:,4);
    
    if (flgUseWeights)
      weightVect = [fitVals(:,6)./fitVals(:,7)]
      if isempty(flgRobust)
        defFit = fit(defVect,defMean,'poly1','Weight',weightVect); 
        if (flgAstig)
          astigMagFit = fit(defVect,astigMag,'poly1','Weight',weightVect); 
          astigAngFit = fit(defVect,astigAng,'poly1','Weight',weightVect); 
        end
      else
        defFit = fit(defVect,defMean,'poly1','Weight',weightVect,'Robust',flgRobust);
        if (flgAstig)
          astigMagFit = fit(defVect,astigMag,'poly1','Weight',weightVect,'Robust',flgRobust); 
          astigAngFit = fit(defVect,astigAng,'poly1','Weight',weightVect,'Robust',flgRobust); 
        end
      end
    else
      if isempty(flgRobust)
        defFit = fit(defVect,defMean,'poly1');   
        if (flgAstig)
          astigMagFit = fit(defVect,astigMag,'poly1'); 
          astigAngFit = fit(defVect,astigAng,'poly1'); 
        end
      else
        defFit = fit(defVect,defMean,'poly1','Robust',flgRobust);  
        if (flgAstig)
          astigMagFit = fit(defVect,astigMag,'poly1','Robust',flgRobust); 
          % The fits near +/- 90 can be pretty ambiguous and throw off the overall trend
          % Try all combinations of values between abs(80-90) and pick the one with the best rSquare overall
          % TODO add some indication of changes if they are made.
          posVect = astigAng > 80.0;
          negVect = astigAng <-80.0;
          positive_score = gather(sum(astigAng(posVect))./sum(posVect));
          negative_score = gather(sum(astigAng(negVect))./sum(negVect));
          if (positive_score > negative_score)
            ambigVect = find(posVect);
          else
            ambigVect = find(negVect);
          end
          [astigAngFit,gof] = fit(defVect,astigAng,'poly1','Robust',flgRobust);
          rSQ = gof.rsquare;
          if isempty(ambigVect)
            nTrials = 0;
          else
            nTrials = nchoosek([ambigVect;-1.*ambigVect],length(ambigVect));
          end
          bestVect = 0;
          if (nTrials(1))
            
            for iTrial = 1:size(nTrials,1)
              testVect = astigAng;
              if numel(nTrials > 1)
                iIDX = nTrials(iTrial,:)';
                testVect(abs(iIDX)) = testVect(abs(iIDX)) .* sign(iIDX);
              end
              [testFit,gof] = fit(defVect,testVect,'poly1','Robust',flgRobust);
              if (gof.rsquare > rSQ)
                fprintf('Improved rSq from %f to %f\n',rSQ,gof.rsquare);
                rSQ = gof.rsquare;
                astigAngFit = testFit;
                bestVect = iIDX;
              end
            end
          end 
          if (bestVect(1))
            astigAng(abs(bestVect)) = astigAng(abs(bestVect)) .* sign(bestVect);
          end
        end
      end
    end
        
    
    if (flgAstig && nStrips > 3)
      figure('visible','off'); title(sprintf('Fits for tilt %6.2f',tiltAngle(ceil(nStrips/2))));
      ax1 = subplot(3,1,1);
      plot(ax1,defVect,defMean,defVect,defFit(defVect));
      title(ax1,sprintf('Mean defocus\n%6.2f x + %2.2f',defFit.p1,defFit.p2));
      ylabel(ax1,'Angstrom')

      ax2 = subplot(3,1,2);
      plot(ax2,defVect,astigMag,defVect,astigMagFit(defVect));
      title(ax2,sprintf('Astigmatism magnitude (df1-df2)/2\n%6.2f x + %6.2f',astigMagFit.p1,astigMagFit.p2))
      ylabel(ax2,'Angstrom')
      
      ax3 = subplot(3,1,3);
      plot(ax3,defVect,astigAng,defVect,astigAngFit(defVect));
      title(ax3,sprintf('Astigmatism angle\n%6.2f x + %6.2f',astigAngFit.p1,astigAngFit.p2))
      ylabel(ax3,'Degrees')
         
      saveas(gcf,sprintf('%s/%s_%2.2f.pdf',tiltBaseName{:},tiltAngle(iStrip)));
%       figure, plot(defVect,defMean,defVect,defFit(defVect));
%       figure, plot(defVect,astigMag,defVect,astigMagFit(defVect));
%       figure, plot(defVect,astigAng,defVect,astigAngFit(defVect));
    else
%       figure, plot(defVect,defMean,defVect,defFit(defVect));
    end
  end
end
