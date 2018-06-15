function [  ] = BH_ctf_Estimate(varargin)

!mkdir -p aliStacks

if length(varargin) == 5
  % PARAMETER_FILE, STACK, TILT, NAMEOUT, mapBack, collectionORDER, gpuIDX
  % orig approach, streamlining for publication
  PARAMETER_FILE = varargin{1};
  stackNameIN = varargin{2};
  tltFile = varargin{3};
  stackNameOUT = varargin{4};
  mapBack = varargin{5};
  collectionORDER = varargin{6};

elseif length(varargin) == 2
  % PARAMETER_FILE, STACK_BASENAME, gpuIDX
  % Assumes that all required files are in ./fixedStacks
  PARAMETER_FILE = varargin{1};
  STACK_BASENAME = varargin{2};

  
  stackNameIN = sprintf('fixedStacks/%s.fixed', STACK_BASENAME);
  tltFile = sprintf('fixedStacks/%s.tlt', STACK_BASENAME);
  stackNameOUT = sprintf('%s_ali1', STACK_BASENAME);
  mapBack = sprintf('fixedStacks/%s', STACK_BASENAME);
  collectionORDER = sprintf('fixedStacks/%s.order',STACK_BASENAME);
  
end
pBH = BH_parseParameterFile(PARAMETER_FILE);

gpuIDX = BH_multi_checkGPU(-1);
gDev = gpuDevice(gpuIDX);

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



tltOrder = load(collectionORDER);


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
defEST = -1.*pBH.('defEstimate').*10^6
defWIN =  pBH.('defWindow').*10^6
tiltRange = [-1];

backGroundBuffer = 0.9985;

try
  deltaZTolerance = pBH.('deltaZTolerance');
catch
  deltaZTolerance = 100e-9;
end

try 
  zShift = pBH.('zShift');
catch
  zShift = 0;
end

% Starting at +/- 100nm
deltaZTolerance = deltaZTolerance / PIXEL_SIZE;
% Use to check for proper gradient.
zShift = zShift / PIXEL_SIZE;

% Tile size & overlap
tileOverlap = 4;
tileSize = floor(680e-10 / PIXEL_SIZE);
tileSize = tileSize + mod(tileSize,2);
fprintf('Using a tile size of %d\n',tileSize);
overlap = floor(tileSize ./ tileOverlap);

% Size to padTile to should be even, large, and preferably a power of 2
try
  paddedSize = pBH.('paddedSize');
catch
  paddedSize = 1024;
end
padVAL = BH_multi_padVal([tileSize,tileSize], [paddedSize,paddedSize]);



if exist(tltFile, 'file') 
  rawTLT = load(tltFile);
else
  fprintf('ignoring %s, because the file is not found.\n', tltFile)
end

if exist(stackNameIN, 'file')
  
  if length(tltOrder) ~= length(rawTLT)
    error('The length of the collectionOrder %d and tilt Geometry %d are different\n', ...
                                        length(tltOrder),length(tltInfo));
  end 
  TLT = zeros(length(rawTLT),22);
  TLT(:,1) = 1:length(rawTLT);
  TLT(:,4) = rawTLT; clear rawTLT
  [pathName,fileName,extension] = fileparts(stackNameIN);
  if isempty(pathName)
    pathName = '.';
  end
else
  fprintf('ignoring %s, because the file is not found.\n', stackNameIN)

end

  


% Make ctf directory to store diagnostic images
system(sprintf('mkdir -p %s/ctf', pathName)); 

if ~(flgResume)



  iMrcObj = MRCImage(stackNameIN,0);

  % The pixel size should be previously set correctly, but if it is not, then we
  % must maintain whatever is there in case beads are to be erased. The model
  % used for this process depends on the pixel size in the header when it was
  % created in IMod alignment.

  iHeader = getHeader(iMrcObj);
  iPixelHeader = [iHeader.cellDimensionX/iHeader.nX .* scalePixelsBy, ...
                  iHeader.cellDimensionY/iHeader.nY .* scalePixelsBy, ...
                  iHeader.cellDimensionZ/iHeader.nZ];
                
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
TLT(:,5:14) = repmat([0,90.0,1.0,0.0,0.0,1.0,1,0,0,1],size(TLT,1),1);
% Defocus will go at 15 - 12 and 13 currently unused.
TLT(:,16:19) = repmat([PIXEL_SIZE,Cs,WAVELENGTH,AMPCONT],size(TLT,1),1);
TLT(:,20:22) = repmat([d1,d2,d3],size(TLT,1),1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Part of the switch to listing dose not bfactor to use the optimal exposure
% filter by Grant/Grigorieff - 
  nPrjs = size(TLT,1);


if CUM_e_DOSE < 0
  % This is a dose per tilt at 0degress to be scaled by the cosine of the
  % tilt angle
  flgCosineDose = 1
  exposure = abs(CUM_e_DOSE)
else
  flgCosineDose = 0
  exposure = CUM_e_DOSE./nPrjs
end

totalExposure = 0;
for iExposure = 1:nPrjs
  % find the projection angle most closley matching (
  [~,iTilt] = min(abs(TLT(:,4)-tltOrder(iExposure)));
  if flgCosineDose == 0
    totalExposure = totalExposure + exposure;
  else
    totalExposure = totalExposure + exposure/cosd(TLT(iTilt,4));
  end
  TLT(iTilt,11) = totalExposure;

end
  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~(flgSkip)

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


  % the projection image resampling, set mag to 1.0
  TLT(i,14) = 1.00;


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

     iProjection = single(getVolume(iMrcObj,-1,-1,TLT(i,1)));

    % Information beyond the physical nyquist should be removed to limit
    % aliasing of noise prior tto interpolation.
    iProjection = BH_padZeros3d(iProjection,[0,0],[0,0],shiftMETHOD,'singleTaper');
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
                 single(getVolume(iMrcObj,[1+osX,d1],[1+osY,d2],TLT(i,1)));
             

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
    system(sprintf('imodtrans -i fixedStacks/%s.fixed fixedStacks/%s.erase fixedStacks/%s.tmpMod',fileName,fileName));
    system(sprintf('model2point -float fixedStacks/%s.tmpMod  fixedStacks/%s.erase2',fileName,fileName));
    system(sprintf('rm fixedStacks/%s.tmpMod',fileName));
    %beadList = importdata(sprintf('fixedStacks/%s.erase2',fileName));
    beadList = load(sprintf('fixedStacks/%s.erase2',fileName));
   
    beadList(:,1:2) = beadList(:,1:2) ./ scalePixelsBy;
    STACK = BH_eraseBeads(STACK,eraseRadius, beadList);
end 
SAVE_IMG(MRCImage(STACK),outputStackName,iPixelHeader,iOriginHeader);

else
  STACK = getVolume(MRCImage(sprintf('aliStacks/%s%s',stackNameOUT,extension)));
end

gpuDevice(gpuIDX)
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
%  [ evalMask, ~ ] = BH_multi_projectionMask([d1C,d2C,d3;d1C,d2C,1], tltForExp, ...
%                                       prjMaskMethod, [zShift,deltaZTolerance] );
evalMask = zeros(d1C,d2C,d3,'single');
for iPrj = 1:d3
  tmpTLT = tltForExp(iPrj,:);
  % need to write over the projections position in the stack to not expand beyond 2d
  tmpTLT(1) = 1;
  [ iEvalMask, ~ ] = BH_multi_projectionMask([d1C,d2C,1;d1C,d2C,1], tmpTLT, ...
                                       'GPU', [zShift,deltaZTolerance] ); 
 
  evalMask(:,:,tltForExp(iPrj,1)) = gather(iEvalMask);
end


%evalMask = gather(evalMask);


nTiles = zeros(size(STACK,3),1);

tic
for i = 1+tileSize/2:overlap:d1C-tileSize/2
  for j = 1+tileSize/2:overlap:d2C-tileSize/2
    for k = 1:size(STACK,3)
      if evalMask(i,j,k)
        nTiles(k) = nTiles(k) + 1;
      end
    end
  end
end
toc
sum(nTiles)




[radialForCTF,phi,~,~,~,~] = ...
   BH_multi_gridCoordinates([paddedSize,paddedSize,1],'Cylindrical','GPU', ...
                                                      {'none'},1,0,0);

radialForCTF = {radialForCTF./(pixelOUT.*10^-10),0,phi}; clear phi                                                       

flgExpFilter = 1
%if ( flgExpFilter)
%  [exposureFilter] = BH_exposureFilter(paddedSize.*[1,1], tltForExp,'GPU',1,0);
%  exposureFilter = gather(exposureFilter);
%else
%  exposureFilter = ones([paddedSize.*[1,1],d3],'single');
%end
inc = (0.5 - FIXED_FIRSTZERO) / (paddedSize/2);
freqVector = [inc+FIXED_FIRSTZERO:inc:0.5 ];
clear sumVector radialAvg
sumVector(length(freqVector)) = gpuArray(double(0));
radialAvg(length(freqVector)) = gpuArray(double(0));



tic
nT = 1;

psTile = zeros((paddedSize).*[1,1],'single','gpuArray');


doseWeight = zeros(paddedSize.*[1,1], 'single','gpuArray');
for k = 1:size(STACK,3)
  doseFilter = BH_exposureFilter(paddedSize.*[1,1], tltForExp(k,:),'GPU',1,0);
   tmpTile = zeros(paddedSize.*[1,1],'single','gpuArray');
%     iProjection = double(gpuArray(STACK(:,:,TLT(k,1))));
  if flgCrop
    [iProjection,~] = cropIMG(gpuArray(STACK(:,:,TLT(k,1))),PIXEL_SIZE*10^10);
  else
    iProjection = (gpuArray(STACK(:,:,TLT(k,1))));
  end
    
  iProjection = iProjection - ...
                       BH_movingAverage(iProjection,[tileSize,tileSize]);
                  
  iProjection = iProjection ./ ...
                       BH_movingRMS(iProjection,[tileSize,tileSize]);


    %doseFilter = gpuArray(exposureFilter(:,:,TLT(k,1)));
    % Weight according to the actual contribution
     doseWeight = doseWeight + nTiles(TLT(k,1)).*doseFilter;
 

  for i = 1+tileSize/2:overlap:d1C-tileSize/2
   for j = 1+tileSize/2:overlap:d2C-tileSize/2    
      if evalMask(i,j,TLT(k,1))



         tmpTile = tmpTile + abs(fftn( ...
                                 BH_padZeros3d(...
                                 (iProjection( ...
                                        i-tileSize/2+1:i+tileSize/2,...
                                        j-tileSize/2+1:j+tileSize/2)),...
                                                  padVAL(1,:),padVAL(2,:),...
                                                  'GPU','singleTaper')));



         nT = nT+1;

      end
    end
  end
  % Apply the dose filter to the sum of each projection to save a bunch of
  % multiplicaiton
  psTile = psTile + tmpTile .* doseFilter; 

end
toc


% Normalize the final weighted average of the dose weights
doseWeight = doseWeight ./ sum(nTiles(:));
rotAvgPowerSpec = (fftshift(psTile));
AvgPowerSpec = rotAvgPowerSpec;

% doseWeight = ifftshift(gather(doseWeight./sum(nTiles(:,1))));
doseWeight = fftshift(gather(doseWeight));


	
% doseWeight = ifftshift(gather(doseWeight./sum(nTiles(:,1))));

[rot1, rot2, ~, r1,r2, ~] = BH_multi_gridCoordinates(paddedSize.*[1,1], ...
                                                    'Cartesian','GPU', ...
                                                    {'none'},0,1,0);

for i = 0.5:0.5:360
  R = BH_defineMatrix([i,0,0],'Bah','forward');
  ROT1 = R(1).*rot1 + R(4).*rot2;
  ROT2 = R(2).*rot1 + R(5).*rot2;
  
  rotAvgPowerSpec = rotAvgPowerSpec + interpn(r1,r2,AvgPowerSpec,...
                                              ROT1,ROT2,'linear',0);
end
 clear ROT1 ROT2
% Normalize on avgerage #, doseWeighting, and a radial filter to account
% for rotational averaging.
rotAvgPowerSpec = rotAvgPowerSpec ./ (720.*ifftshift(doseWeight).*(sqrt(fftshift(radialForCTF{1}.*(pixelOUT.*10^-10))))); clear a
% rotAvgPowerSpec = rotAvgPowerSpec ./ (720.*(sqrt(fftshift(radialForCTF{1}.*(pixelOUT.*10^-10))))); clear a

AvgPowerSpec = AvgPowerSpec ./ ifftshift(doseWeight);

checkInfNan = (isnan(rotAvgPowerSpec) | isinf(rotAvgPowerSpec));
rotAvgPowerSpec(checkInfNan) =  mean(mean(rotAvgPowerSpec(~checkInfNan)));

radialAvg = [rotAvgPowerSpec((paddedSize/2)+1:end,(paddedSize/2)+1)]';
radialPS = [AvgPowerSpec((paddedSize/2)+1:end,(paddedSize/2)+1)]';


defRange = [defEST-defWIN,defEST+defWIN]

if defRange(2) > -0.05
  fprintf('\n\nCapping defocus to 50nm from wanted %f. Do you mean to search so close to focus??\n\n',abs(defRange(2)));
  defRange(2) = -0.05;
end
defInc   = [0.01];

defVal = (defRange(1):defInc:defRange(2))';

cccStorage = zeros(length(defVal),2);
cccStorage(:,1) = defVal;
nDF = 1;
for iDF = defVal'
  DF = iDF*10^-6

  [ Hqz ] = BH_ctfCalc(radialForCTF,Cs,WAVELENGTH,DF,paddedSize,-AMPCONT,-1.0);

  [ bg,  bandpass, rV ] = prepare_spectrum( Hqz, highCutoff, freqVector, radialAvg, 0);

  [ iCCC ] = calc_CCC( freqVector, bg,  bandpass, radialAvg, rV, cccScale);

  cccStorage(nDF, 2) = gather(iCCC);
  nDF = nDF +1;
end

[~,maxVal] = max(cccStorage(:,2));
maxDef = cccStorage(maxVal,1)

figure('Visible','off'), scatter(cccStorage(:,1), cccStorage(:,2));
title(sprintf('CCC\nmax = %03.3f micron', maxDef));
xlabel('defocus (micron)'); ylabel('CCC');

  saveas(gcf,sprintf('%s/ctf/%s_ccFIT.pdf',pathName,stackNameOUT), 'pdf')

  DF = maxDef*10^-6;

  [ Hqz ] = BH_ctfCalc(radialForCTF,Cs,WAVELENGTH,DF,paddedSize,-AMPCONT,-1.0);

  [ bg, bandpass, rV ] = prepare_spectrum( Hqz, highCutoff, freqVector, radialAvg, 0);



figure('Visible','off'), plot(freqVector(bandpass),backGroundBuffer.*bg(freqVector(bandpass)),freqVector(bandpass),abs(radialAvg(bandpass)));
title('Background fitting')

  saveas(gcf,sprintf('%s/ctf/%s_bgFit.pdf',pathName,stackNameOUT), 'pdf')


% [ diagnosticIMG ] = make_diagnosticIMG( Hqz, pixelOUT, bandpass, bg, {rotAvgPowerSpec});
% 
% SAVE_IMG(MRCImage(diagnosticIMG), ...
%                        sprintf('%s/ctf/%s_diag%s',pathName,fileName,extension));

                     clear STACK exposureFilter evalMask
                     
 pdfOUT = sprintf('%s/ctf/%s_psRadial.pdf',pathName,stackNameOUT)

  
  bgSubPS = (abs(radialAvg) - bg(freqVector)').*bandpass;
  bgSubPS = bgSubPS ./ max(bgSubPS(:));
  figure('Visible','off'), plot(freqVector(bandpass)./(pixelOUT),bgSubPS(bandpass), freqVector(bandpass)./(pixelOUT),abs(rV(bandpass)).^2./max(abs(rV(bandpass)).^2),'-g');
  title(sprintf('CTF fit\n%03.3f μm ', maxDef));
  xlabel('Spatial Frequency (1/Å)'); ylabel('Relative Power');

  saveas(gcf,pdfOUT, 'pdf')
%   save('resume.mat');
else
%   load('resume.mat');
  
  
end %%% temp to troubleshoot  
  if (flgAstigmatism)
    
    radialForCTF = {fftshift(radialForCTF{1}),1,fftshift(radialForCTF{3})};  
    [radialAstig,~,~,~,~,~] = ...
                         BH_multi_gridCoordinates(size(Hqz),'Cartesian',...
                                                     'GPU',{'none'},1,1,1);
  
    % Hqz from max defocus still exisists
    
    [ bg, bandpass, ~ ] = prepare_spectrum( Hqz, highCutoff ,...
                                            freqVector, AvgPowerSpec, radialAstig);
                                          
    [ bgSubPS ] = calc_CCC( radialAstig, bg, bandpass, AvgPowerSpec,-9999, cccScale) ; 
    
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
        df1 =  maxDef*10^-6 - iDelDF*astigStep;
        df2 =  maxDef*10^-6 + iDelDF*astigStep;

         
        [ Hqz ] = BH_ctfCalc(radialForCTF,Cs,WAVELENGTH, ...
                                [df1,df2,iAng],size(radialForCTF{1}),-AMPCONT,-1.0);

%         [ bg, bandpass, rV ] = prepare_spectrum( Hqz, highCutoff ,...
%                                              freqVector, radialPS, radialAstig); 
            
                
        [ iCCC ] = calc_CCC( radialAstig, bgSubPS, bandpass, AvgPowerSpec, Hqz, cccScale)                                            
       


        initAstigCCC(n,:) = [iAng,iDelDF*astigStep,iCCC]; 
        n = n + 1;                 
        fprintf('%d / %d coarse astigmatism search\n',n,size(initAstigCCC,1));
      end
    end


    nPeaks = 1;
    top3 = zeros(nPeaks,3,'gpuArray');

    for iCCC = 1:nPeaks
      [~,c] = max(initAstigCCC(:,3));
      top3(iCCC,:) = initAstigCCC(c,:);
      initAstigCCC = initAstigCCC(initAstigCCC(:,1)~=initAstigCCC(c,1),:);
    end

    top3
    
    refineCCC = zeros(length(fineDefSearch)*length(fineAngSearch)*nPeaks,3,'gpuArray');

    n=1;
    for iPeak = 1:nPeaks
      mAng = top3(iPeak,1);
      mDef = top3(iPeak,2);

      for iAng = fineAngSearch   
        for iDelDF = fineDefSearch

          df1 =  maxDef*10^-6 - (mDef + iDelDF);
          df2 =  maxDef*10^-6 + (mDef + iDelDF);

          % For values very close to zero, the search range may include
          % values |df1| < |df2| which is against convention.
          if abs(df1) >= abs(df2)
          

            [ Hqz ] = BH_ctfCalc(radialForCTF,Cs,WAVELENGTH, ...
                                [df1,df2,iAng+mAng], ...
                                size(radialForCTF{1}), -AMPCONT,-1.0);
                                           
          
            [ iCCC ] = calc_CCC( radialAstig,bgSubPS, bandpass, ...
                                 AvgPowerSpec, Hqz,cccScale )   
      
          else
            iCCC = -9999
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
  


      
  % Add the determined defocus, and write out with mic paramters as well.
  TLT(:,15) = repmat(maxDef*10^-6,size(TLT,1),1);
  if (flgAstigmatism) && (refineCCC(c,3)~=-9999)
    TLT(:,12) = repmat(gather(refineCCC(c,2)),size(TLT,1),1);
    TLT(:,13) = repmat(gather(refineCCC(c,1)),size(TLT,1),1);
  end

  % Sort descending along the magnitude of the tilt angles because higher tilts take
  % longer on CTF correction. If more processor available than projections,
  % this doesn't affect anything.
  [~, idx] = sortrows(abs(TLT(:,4)), -1);
  TLT = TLT(idx,:);
  % number in stack, dx, dy, tilt angle, projection rotation, tilt azimuth, tilt
  % elevation, e1,e2,e3, dose number (order in tilt collection), offsetX, offsetY
  % scaleFactor, defocus, pixelSize, CS, Wavelength, Amplitude contrast
  fileID = fopen(sprintf('%s/ctf/%s_ctf.tlt',pathName,stackNameOUT), 'w');
  fprintf(fileID,['%d\t%08.2f\t%08.2f\t%07.3f\t%07.3f\t%07.3f\t%07.7f\t%07.7f\t',...
           '%07.7f\t%07.7f\t%5e\t%5e\t%5e\t%7e\t%5e\t%5e\t%5e\t%5e\t%5e\t',...
           '%d\t%d\t%d\n'], TLT');




end

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

function [ diagnosticIMG ] = make_diagnosticIMG( Hqz, pixelSize, bandpass, bg, IMG)

  
  iImg = 1;
  Hqz = fftshift(Hqz);
  paddedSize = size(Hqz,1);
  [radialGrid,~,~,~,~,~] = ...
                         BH_multi_gridCoordinates(size(Hqz),'Cartesian',...
                                                     'GPU',{'none'},1,0,1);
  
  radialGrid = radialGrid ./ pixelSize;
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

