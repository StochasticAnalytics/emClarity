function []  = BH_templateSearch3d_2( PARAMETER_FILE,...
                                        tomoName,tomoNumber,TEMPLATE, ...
                                        SYMMETRY, wedgeType, varargin)
 
                                                       
%3d template matching

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




if length(varargin) == 1
  % Allow for an override of the max number, useful when only a few tomos
  % have a strong feature like carbon that is hard to avoid.
  cmdLineThresh = 0;
  gpuIDX = str2num(varargin{1});
elseif length(varargin) == 2
  cmdLineThresh = str2num(varargin{1});
  gpuIDX = str2num(varargin{2});
end
  tomoNumber = str2num(tomoNumber);


  [ useGPU ] = BH_multi_checkGPU( gpuIDX )



gpuDevice(useGPU);
  
SYMMETRY = str2num(SYMMETRY);
startTime = clock ;

pBH = BH_parseParameterFile(PARAMETER_FILE);
try
  load(sprintf('%s.mat', pBH.('subTomoMeta')), 'subTomoMeta');
  mapBackIter = subTomoMeta.currentTomoCPR
%   clear subTomoMeta
  % Make sure we get a CTF corrected stack
  shouldBeCTF = 1
catch
  mapBackIter = 0;
  shouldBeCTF = -1;
end
samplingRate  = pBH.('Tmp_samplingRate');

try
  tmpDecoy = pBH.('templateDecoy')
catch
  tmpDecoy = 0
end

try
  super_sample = pBH.('super_sample');
  if (super_sample > 0)
    [~,v] = system('cat $IMOD_DIR/VERSION');
    v = split(v,'.');
    if (str2num(v{1}) < 4 || (str2num(v{2}) <= 10 && str2num(v{3}) < 42))
      fprintf('Warning: imod version is too old for supersampling\n');
      super_sample = '';
    else
      super_sample = sprintf(' -SuperSampleFactor %d',super_sample);
    end
  else
    super_sample = '';
  end
catch
  super_sample = '';
  expand_lines = '';
end

if ( cmdLineThresh )
 peakThreshold = cmdLineThresh;
 fprintf('\nOverride peakThreshold from paramfile (%d) with cmd line arg (%d)\n\n',...
         cmdLineThresh, pBH.('Tmp_threshold'));
else
 peakThreshold = pBH.('Tmp_threshold');
end

latticeRadius = pBH.('particleRadius');
try
  targetSize    = pBH.('Tmp_targetSize')
catch
  targetSize = [512,512,512];
end
angleSearch   = pBH.('Tmp_angleSearch');

statsRadius = 1;

convTMPNAME = sprintf('convmap_wedgeType_%d_bin%d',wedgeType,samplingRate)

try 
  eraseMaskType = pBH.('Peak_mType');
catch
  eraseMaskType = 'sphere';
end
try
  eraseMaskRadius = pBH.('Peak_mRadius');
catch
  eraseMaskRadius = 0.75.*latticeRadius;
end


nPreviousSubTomos = 0;

reconScaling = 1;
try
  nPeaks = pBH.('nPeaks');
catch
  nPeaks = 1;
end

ignore_threshold = false;
try
  max_tries = pBH.('max_peaks');
  if max_tries < 0
    max_tries = abs(max_tries);
    ignore_threshold = true;
  end
catch
  max_tries = 10000;
end

pixelSizeFULL = pBH.('PIXEL_SIZE').*10^10;
if pBH.('SuperResolution')
  pixelSizeFULL = pixelSizeFULL * 2;
end

pixelSize = pixelSizeFULL.*samplingRate;

% For testing
try 
  wantedCut = pBH.('lowResCut');
catch
  wantedCut = 28;
end

firstZero = 0;
% Limit to the first zero if we are NOT using the CTF rec
if (shouldBeCTF ~= 1)
TLT = load(sprintf('fixedStacks/ctf/%s_ali%d_ctf.tlt',tomoName,mapBackIter+1));
  def = mean(-1.*TLT(:,15))*10^6; %TODO if you switch to POSITIVEDEFOCUS this will be wrong
  firstZero = -0.2*def^2 +5.2*def +11;

  % Take the lower of firstZero lowResCut or Nyquist
  lowResCut = max(wantedCut, firstZero);
else
  lowResCut = wantedCut;
end

if pixelSize*2 > lowResCut
  fprintf('\nLimiting to Nyquist (%f) instead of user requested lowResCut %f Angstrom\n',pixelSize*2,lowResCut);
  lowResCut = pixelSize*2;
else
  fprintf('\nUsing max (%f) of specified resolution cutoff of %f and first ctf zero %f Angstrom\n',lowResCut, wantedCut, firstZero);
end


mapPath = './cache';
mapName = sprintf('%s_%d_bin%d',tomoName,tomoNumber,samplingRate);
mapExt = '.rec';

sprintf('recon/%s_recon.coords',tomoName)
[ recGeom, ~, ~] = BH_multi_recGeom( sprintf('recon/%s_recon.coords',tomoName) );

reconCoords = recGeom(tomoNumber,:);
clear recGeom


%  [ tomogram ] = BH_multi_loadOrBuild( sprintf('%s_%d',tomoName,tomoNumber),  ...
%                                       reconCoords, mapBackIter, samplingRate,...
%                                       shouldBeCTF*gpuIDX, reconScaling,1,'','ctf'); 

[ tomogram ] = BH_multi_loadOrBuild( sprintf('%s_%d',tomoName,tomoNumber),  ...
                                    reconCoords, mapBackIter, samplingRate,...
                                    shouldBeCTF*gpuIDX, reconScaling,1,'',super_sample); 
                                           

% We'll handle image statistics locally, but first place the global environment
% into a predictible range

  

[template, tempPath, tempName, tempExt] = ...
                              BH_multi_loadOrBin( TEMPLATE, 1, 3 ); 
                            

                            

                            
% The template will be padded later, trim for now to minimum so excess
% iterations can be avoided.
fprintf('size of provided template %d %d %d\n',size(template));
trimTemp = BH_multi_padVal(size(template),ceil(2.0.*max(pBH.('Ali_mRadius')./pixelSizeFULL)));
% template = BH_padZeros3d(template, trimTemp(1,:),trimTemp(2,:),'cpu','singleTaper');
% SAVE_IMG(MRCImage(template),'template_trimmed.mrc');
clear trimTemp
fprintf('size after trim to sqrt(2)*max(lattice radius) %d %d %d\n',size(template));
                            
if isempty(mapPath) ; mapPath = '.' ; end
if isempty(tempPath) ; tempPath = '.' ; end
% Check to see if only tilt angles are supplied, implying a y-axis tilt scheme,
% or otherwise, assume a general geometry as in protomo.
% % % tiltGeometry = load(RAWTLT);
RAWTLT = sprintf('fixedStacks/ctf/%s_ali1_ctf.tlt',tomoName);
tiltGeometry = load(RAWTLT);
% subTomoMeta.('tiltGeometry').(mapName) = tiltGeometry;

% Make sure the template and is an even sized image
template = padarray(template, mod(size(template),2),0, 'post');
template = template - mean(template(:));

templateBIN = BH_reScale3d(template,'',sprintf('%f',1/samplingRate),'cpu');
templateBIN = templateBIN - mean(templateBIN(:));
templateBIN = templateBIN  ./rms(templateBIN(:));


sizeTemp = size(template)
sizeTempBIN = size(templateBIN)



statsRadiusAng = statsRadius.*[2,2,2].*max(latticeRadius);
statsRadius = ceil(statsRadiusAng./pixelSize);
latticeRadius = (0.75 .* latticeRadius) ./ (pixelSize);
latticeRadius = floor(latticeRadius);
latticeRadius = latticeRadius + mod(latticeRadius, 2);

eraseMaskRadius = floor((eraseMaskRadius) ./ (pixelSize));
eraseMaskRadius = eraseMaskRadius + mod(eraseMaskRadius,2);

fprintf('\ntomograms normalized in %f Angstrom cubic window\n',statsRadiusAng(1));

fprintf('\nlatticeRadius = %dx%dx%d pixels\n\n', latticeRadius);
fprintf('\neraseMaskType %s, eraseMaskRadius %dx%dx%d pixels\n',eraseMaskType,eraseMaskRadius);
  % For wedgeMask
particleThickness =  latticeRadius(3);

gpuDevice(useGPU);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize a whole mess of control variables and storage volumes. %
%Out of plane range inc (starts from 1.* inc)
if length(angleSearch) == 5
  helical = angleSearch(5);
else
  helical = 0;
end

[  nInPlane, inPlaneSearch, angleStep, nAngles] ...
                                      = BH_multi_gridSearchAngles(angleSearch)

                                                                        

highThr=sqrt(2).*erfcinv(ceil(peakThreshold.*0.025).*2./(prod(size(tomogram)).*nAngles(1)))

[ OUTPUT ] = BH_multi_iterator( [targetSize; ...
                                 size(tomogram);...
                                 sizeTempBIN; ...
                                 2.*latticeRadius], 'convolution' );


    
tomoPre   = OUTPUT(1,:);
tomoPost  = OUTPUT(2,:);
sizeChunk = OUTPUT(3,:);
validArea = OUTPUT(4,:);
validCalc = OUTPUT(5,:);
nIters    = OUTPUT(6,:);

%[ padVal ] = BH_multi_padVal( sizeTemp, sizeChunk );
%tempPre = padVal(1,:);
%tempPost = padVal(2,:);

[ padBIN ] = BH_multi_padVal( sizeTempBIN, sizeChunk );
[ trimValid ] = BH_multi_padVal(sizeChunk, validArea);

if ( tmpDecoy )
  % This is probably sample dependent. should search a small range and find
  % the maximum rate of change in the ccc
  
  % the -1 searches for the next smallest fast fourier size
  templateBIN = gpuArray(templateBIN);


  
  decoyTest = BH_reScale3d(templateBIN,'',tmpDecoy,'GPU');
  decoyTrim = BH_multi_padVal(size(decoyTest),size(templateBIN));
  decoyTest = fftn(BH_padZeros3d(decoyTest,decoyTrim(1,:),decoyTrim(2,:),'GPU','single'));
  decoyShift = -1.*gather(BH_multi_xcf_Translational(decoyTest,conj(fftn(templateBIN)),'',[3,3,3]));  
  decoyNorm = gather(sum(abs(decoyTest(:)))./sum(abs(fftn(templateBIN(:)))));
  padDecoy = BH_multi_padVal(size(decoyTest),sizeChunk) + decoyTrim;
  clear decoyTest
  templateBIN = gather(templateBIN);
  fprintf('tmpDecoy %f normFactor %f and shift by %2.2f %2.2f %2.2f\n',tmpDecoy,decoyNorm,decoyShift);

end


fprintf('\n-----\nProcessing in chunks\n\n');
fprintf('tomo prepadding  %d %d %d\n', tomoPre);
fprintf('tomo postpadding %d %d %d\n', tomoPost);
fprintf('size to process  %d %d %d\n', sizeChunk);
fprintf('valid Area       %d %d %d\n', validArea);
fprintf('valid Calc       %d %d %d\n', validCalc);
fprintf('# of iterations  %d %d %d\n', nIters);
fprintf('-----\n');

size(tomogram)

[ tomogram ] = BH_padZeros3d(tomogram, tomoPre, tomoPost, ...
                                             'cpu', 'singleTaper',mean(tomogram(:)));
% tomogram = padarray(tomogram,tomoPre,'symmetric','pre');
% tomogram = padarray(tomogram,tomoPost,'symmetric','post');
sizeTomo = size(tomogram);


[ validAreaMask ] = gather(BH_mask3d('rectangle',sizeChunk,validCalc./2,[0,0,0]));
[ vA ] = BH_multi_padVal( validArea, sizeChunk );
% This would need to be changed to take a mask size and not just a radius.
% Currently, this would not produce the correct results for odd size area
% % % fftMask = BH_fftShift(validArea,sizeChunk,0);

% Array for storing chunk results
RESULTS_peak = zeros(sizeTomo, 'single'); 
RESULTS_angle= zeros(sizeTomo, 'single');
if ( tmpDecoy )
  RESULTS_decoy = RESULTS_peak;
end
% Loop over tomogram
% Set this up second


% % % % optimize fft incase a power of two is not used, this will make things run ok.
% % % opt = zeros(sizeChunk, precision,'gpuArray');
% % % fftw('planner','patient');
% % % fftn(opt);
% % % clear opt ans
[ bhF ] = fourierTransformer(randn(sizeChunk, 'single','gpuArray'));


% Temp while testing new dose weighting
TLT = tiltGeometry;
nPrjs = size(TLT,1);


kVal = 0;

% % [ OUTPUT ] = BH_multi_iterator( [sizeTempBIN;kVal.*[1,1,1]], 'extrapolate' );
[ OUTPUT ] = BH_multi_iterator( [sizeChunk;kVal.*[1,1,1]], 'extrapolate' );


% 
% switch wedgeType
%   case 1
%     % make a binary wedge
%     [ wedgeMask ]= BH_weightMask3d(-1.*OUTPUT(1,:), tiltGeometry, ...
%                    'binaryWedgeGPU',particleThickness,...
%                                                  1, 1, samplingRate);
%   case 2
%     % make a non-CTF wedge
%     [ wedgeMask ]= BH_weightMask3d(-1.*OUTPUT(1,:), tiltGeometry, ...
%                    'applyMask',particleThickness,...
%                                                  2, 1, samplingRate);   
%   case 3
%     % make a CTF without exposure weight
%     [ wedgeMask ]= BH_weightMask3d(-1.*OUTPUT(1,:), tiltGeometry, ...
%                    'applyMask',particleThickness,...
%                                                  3, 1, samplingRate);   
%   case 4
%     % make a wedge with full-ctf
%     [ wedgeMask ]= BH_weightMask3d(-1.*OUTPUT(1,:), tiltGeometry, ...
%                    'applyMask',particleThickness,...
%                                                  4, 1, samplingRate);  
%   otherwise
%     error('wedgeType must be 1-4');
% end
% 
% wedgeMask = (ifftshift(wedgeMask));
%                                      
% % Now just using the mask to calculate the power remaining in the template,
% % without actually applying.
% wedgeMask = gather(find(ifftshift(wedgeMask)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Preprocess the tomogram

tomoIDX = 1;
nTomograms = prod(nIters);


tomoStack = zeros([sizeChunk,nTomograms], 'single');
% tomoNonZero = zeros(nTomograms,6,'uint64');

% backgroundVol = zeros(sizeChunk,'single');
tomoCoords= zeros(nTomograms, 3, 'uint16');

% % % [ tomoBandpass ]   = BH_bandpass3d(sizeChunk, 0,maxSizeForHighPass, ...
% % %                                               lowResCut,'cpu', pixelSize );
% In switching to the full 3D-sampling function the high pass is
% already incorporated in the CTF. Still include one for very low
% resolution to deal with gradients in the tomos.
% [ tomoBandpass ]   = BH_bandpass3d(sizeChunk, 1e-3,600, ...
%                                               lowResCut,'cpu', pixelSize );
%                                           
% % if ~(shouldBeCTF)
% %   tomoBandpass = wedgeMask .* tomoBandpass;
% % end
% 
% tomoBandpass = tomoBandpass(1:floor(size(tomoBandpass,1)/2)+1,:,:);
% bhF.bandpass = tomoBandpass; clear tomoBandpass
% clear wedgeMask

try
  doMedFilt = pBH.('Tmp_medianFilter');
  if ~ismember(doMedFilt,[3,5,7])
    error('Tmp_medianFilter can only be 3,5, or 7');
  else
    fprintf('Using median filter, size %d',doMedFilt);
  end
catch
  doMedFilt =0
end

calcStats = 0;
if calcStats
  maskStack = false([sizeChunk,nTomograms]);
  calcMask = 0;
else
  calcMask = 1;
end
firstStats = 1;
flgOOM = 0;

fullX = 0;
fullX2 = 0;
fullnX = 0;
oT  = ceil((validArea./2)+1)
for  iX = 1:nIters(1)
  cutX = 1 + (iX-1).*validArea(1);
  for iY = 1:nIters(2)
    cutY = 1 + (iY-1).*validArea(2);
    for iZ = 1:nIters(3)
      cutZ = 1 + (iZ-1).*validArea(3);

    fprintf('preprocessing tomo_chunk %d/%d col %d/%d row %d/%d plane idx%d\n' , ...
                                  iY,nIters(2),iX,nIters(1),iZ,nIters(3),tomoIDX)


    
    tomoChunk = gpuArray(tomogram(cutX:cutX+sizeChunk(1)-1,...
                                  cutY:cutY+sizeChunk(2)-1,...
                                  cutZ:cutZ+sizeChunk(3)-1));
                                       
    % Make a list of the padded regions of the tomogram to exclude from
    % statistical calculations
    


% % % % %     tomoChunk = real(ifftn(fftn(tomoChunk).*tomoBandpass));
    tomoChunk = bhF.invFFT(bhF.fwdFFT(tomoChunk,0,0,[1e-3,600, ...
                                              lowResCut, pixelSize]),2);
    

    if doMedFilt
      if ( flgOOM )
        tomoChunk = (medfilt3(tomoChunk,doMedFilt.*[1,1,1]));
      else
        tomoChunk = gpuArray(medfilt3(tomoChunk,doMedFilt.*[1,1,1]));
      end
    else
      if ( flgOOM )
        % Leave on CPU
      else        
        tomoChunk = gpuArray(tomoChunk);
        statsRadius = gather(statsRadius);
      end
    end

%     % Handle all mean centering and rms normalization in local window
%     
%     [ averageMask, flgOOM ] = BH_movingAverage(tomoChunk, statsRadius); 
%    
%     if isa(tomoChunk(1),'gpuArray') && flgOOM
%       tomoChunk = gather(tomoChunk);
%     end
% 
%     
%     
%     tomoChunk= tomoChunk - averageMask; clear averageMask 
% 
%     [ rmsMask ] = BH_movingRMS(tomoChunk, statsRadius);
    
  
      
% % % % %     if ( shouldBeCTF == 1 )
% % % % %       tomoStack(:,:,:,tomoIDX) = gather((1.*tomoChunk ./ rmsMask).*validAreaMask);
% % % % %     else 
% % % % %       % Using the non-ctf corrected stack since we limit toA all practical
% % % % %       % defocus (<8um) should be entirely negative, so just flip in real
% % % % %       % space
% % % % %     
% % % % %       tomoStack(:,:,:,tomoIDX) = gather(-1.*(tomoChunk ./ rmsMask).*validAreaMask);
% % % % % %       backgroundVol = backgroundVol + gather(tomoChunk.*maskStack(:,:,:,tomoIDX));
% % % % %     end

    clear rmsMask

    tomoChunk = tomoChunk .* (-1*shouldBeCTF); % This is backwards, but I don't know why
    fullX = fullX + gather(sum(tomoChunk(:)));
    fullX2 = fullX2 + gather(sum(tomoChunk(:).^2));
    fullnX = fullnX + gather(prod(sizeChunk));

    tomoCoords(tomoIDX,:) = [cutX,cutY,cutZ];
    
    tomoStack(:,:,:,tomoIDX) = gather(tomoChunk);

    tomoIDX = tomoIDX + 1;

      
    end % end of loop over Z chunks
  end % end of loop over Y chunks
end % end of loop over X chunks

% Normalize the global variance
globalVariance = (fullX2 - fullX)/fullnX;
fprintf('After local normalization, scaling also the global variance\n');

for iChunk = 1:tomoIDX-1
  tomoStack(:,:,:,iChunk) = tomoStack(:,:,:,iChunk) ./ globalVariance;
end


  
clear tomoWedgeMask averagingMask rmsMask bandpassFilter statBinary validAreaMask tomoChunk

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

kVal = 0;


    

currentGlobalAngle = 1;
ANGLE_LIST = zeros(nAngles(1),3, 'single');
nComplete = 0;
totalTime = 0;
firstLoopOverTomo = true;

    % Center the spectrum by multiplication not swapping
    [dU,dV,dW] = BH_multi_gridCoordinates(size(tomoStack(:,:,:,1)),...
                                         'Cartesian','GPU', ...
                                          {'none'},1,1,0);
                                        
    swapQuadrants = exp((-2i*pi).*(dU.*(floor(size(dU,1)/2)+1) + ...
                                  (dV.*(floor(size(dV,2)/2)+1) + ...          
                                  (dW.*(floor(size(dW,3)/2)+1)))));  
    clear dU dV dW
    
    swapQuadrants = swapQuadrants(1:floor(size(swapQuadrants,1)/2)+1,:,:);

                                                              
for iAngle = 1:size(angleStep,1)
  
  theta = angleStep(iAngle,1);

  % Calculate the increment in phi so that the azimuthal sampling is
  % consistent and equal to the out of plane increment.
  if (helical)
     phi_step = 360;
  else
     phiStep = angleStep(iAngle,3);
  end
  

  numRefIter = angleStep(iAngle,2)*length(inPlaneSearch)+1;
  tempImg = gpuArray(templateBIN); %%%%% NEW switch to bin

  
  interpolationNormFactor = sum(abs(tempImg(:)).^2);


  clear referenceStack tempFilter
  % Calculate all references for each out of plane tilt only once
  referenceStack = zeros([sizeTempBIN,numRefIter], 'single', 'gpuArray');
                                                         
  tomoIDX = 1;
  firstLoopOverAngle = true;
  % Iterate over the tomogram pulling each chunk one at a time.
  for iTomo = 1:nTomograms
    tic;
    iCut =  tomoCoords(tomoIDX,:);
    % reset the angle count and value at the begining of loop
    % inside, while each new outer loop changes the start values.

%       nAngle =  angleIncStart;
    intraLoopAngle = 1;

    % Truth value to initialize temp results matrix each new tomo
    % chunk.
    firstLoopOverChunk = true;
    
    fprintf('working on tilt(%d/%d) tomoChunk(idx%d/%d)\t' ...
                          ,iAngle,size(angleStep,1), tomoIDX,nTomograms);


 

     tomoFou = swapQuadrants.*bhF.fwdFFT(gpuArray(tomoStack(:,:,:,tomoIDX)));

    for iAzimuth = 0:angleStep(iAngle,2)

      if helical == 1
          phi = 90 ;
      else
         phi = phiStep * iAzimuth;
      end

      for iInPlane = inPlaneSearch
        psi = iInPlane;

        %calc references only on first chunk
        if (firstLoopOverAngle)

          ANGLE_LIST(currentGlobalAngle,:) = [phi, theta, psi - phi];
          [phi, theta, psi - phi];
          % Rotate the reference, lowpass and wedge mask, send to gpu
          % Inverse rotation(i.e. rotate particle, angles saved
          % are to rotate frame to particle for extraction.)
% % % % %           tempRot = BH_resample3d(tempImg, [phi, theta, psi - phi], [1,1,1], ...
% % % % %                                   {'Bah', 1,'linear',1,interpMaskGPU},...
% % % % %                                   'GPU','forward');


          tempRot = BH_resample3d(tempImg, [phi, theta, psi - phi], [1,1,1], ...
                                  {'Bah', 1,'linear',1},...
                                  'GPU','forward');



          %%%%%tempFou = BH_bandLimitCenterNormalize(tempRot,tempWedgeMask,'',[tempPre;tempPost],precisionTaper);

          %%%%%tempRot = BH_padZeros3d(real(ifftn(tempFou)),-1.*tempPre,-1.*tempPost,'GPU',precision); 
          
          %%%%%tempRot = gather(BH_reScale3d(tempRot,'',sprintf('%f',1/samplingRate),'GPU'));

         % if (firstLoopOverTomo)
         %   SAVE_IMG(MRCImage(tempRot), sprintf('temp_%s.mrc',convTMPNAME),pixelSize);
         % end
          
         % First correct for any change in power due to
         % rotation/interpolation
         tempRot = tempRot - mean(tempRot(:));
         tempRot = tempRot ./ ((interpolationNormFactor./sum(abs(tempRot(:)).^2)).*rms(tempRot(:)));

%          tempRot = tempRot .* (interpolationNormFactor./sum(abs(tempRot(:)).^2));
         % Then correct for any change in power due to the wedge. These can be combined 
%           normFT = abs(fftn(tempRot).*tempBnd).^2;
% %        
% %          
%            normScore = sum(normFT(:)) ./ sum(normFT(:).*tempWdg(:));
%            clear normFT;
%            tempRot = tempRot .* normScore;
          %clear normScore 
          
          referenceStack(:,:,:,intraLoopAngle) = tempRot;
         
% % % % %           tempFou = fftn(BH_padZeros3d(tempRot,padBIN(1,:),padBIN(2,:),'GPU',precision));
 
%           tempFou = BH_bandLimitCenterNormalize( tempRot, tempBandpass, '', ...
%                                                   padBIN, 'single' );
% % % % %           tempFou = bhF.fwdFFT(BH_padZeros3d(tempRot,padBIN(1,:),padBIN(2,:),'GPU','single',real(mean(tempRot(:)))));

          tempFou = bhF.fwdFFT(BH_padZeros3d(tempRot,padBIN(1,:),padBIN(2,:),'GPU','single'));

          
%               tempFou = fftn(BH_padZeros3d(tempRot,padBIN(1,:),padBIN(2,:),'GPU','single'));

        else

% % %           tempFou = (fftn(BH_padZeros3d( ...
% % %                               referenceStack(:,:,:,intraLoopAngle), ...
% % %                               padBIN(1,:), padBIN(2,:),'GPU', 'single')));
                            
                             
          tempFou = bhF.fwdFFT(BH_padZeros3d( ...
                              referenceStack(:,:,:,intraLoopAngle), ...
                              padBIN(1,:), padBIN(2,:),'GPU', 'single'));
%           tempFou = fftn(BH_padZeros3d( ...
%                               referenceStack(:,:,:,intraLoopAngle), ...
%                               padBIN(1,:), padBIN(2,:),'GPU', 'single'));


        end




%         ccfmap = BH_padZeros3d(fftshift(real(single(...
%                                bhF.invFFT(tomoFou.*conj(tempFou),2)))),...%./(tomoNorm.*tempNorm))))),...
%                                trimValid(1,:),trimValid(2,:),'GPU','single');  
        ccfmap = BH_padZeros3d((real(single(...
                               bhF.invFFT(tomoFou.*conj(tempFou))))),...%./(tomoNorm.*tempNorm))))),...
                               trimValid(1,:),trimValid(2,:),'GPU','single');
%                              

        ccfmap = ccfmap ./ std(ccfmap(:));
        
        if ( tmpDecoy > 0 )
           tempFou = [];
           if (firstLoopOverAngle)

             decoy = BH_padZeros3d(BH_reScale3d(tempRot./decoyNorm,'',tmpDecoy,'GPU',decoyShift),...
                                   padDecoy(1,:),padDecoy(2,:),'GPU','single');
           else
              % Probably just make a second decoy stack to avoid
              % re-interpolating. If it works, then do this.
             decoy = BH_padZeros3d(BH_reScale3d(referenceStack(:,:,:,intraLoopAngle)./decoyNorm,'',tmpDecoy,'GPU',decoyShift),...
                                   padDecoy(1,:),padDecoy(2,:),'GPU','single');                              
           end
           

         
           decoy = BH_padZeros3d(fftshift(real(single( ...
                                 ifftn(tomoFou.*conj(fftn(decoy)))))),..../(decoyNorm.*tomoNorm))))),
                                 trimValid(1,:), ...
                                 trimValid(2,:),'GPU','single');
                               

        elseif ( tmpDecoy < 0 )
          
          % Just use the mirror image of the template, i.e. take the conj
          % (of the conj) so just the padded FFT of the ref.
           decoy = BH_padZeros3d(fftshift(real(single( ...
                                 ifftn(tomoFou.*tempFou)))),..../(decoyNorm.*tomoNorm))))),
                                 trimValid(1,:), ...
                                 trimValid(2,:),'GPU','single');          
        tempFou = [];
        end
        clear tempRot
        % If first loop over tomo, initialize the storage volumes, if
        % first loop over the chunk but not over the tomo, pull storage
        % chunks from storage volume.
        if (firstLoopOverTomo && firstLoopOverChunk)
            %store ccfmap as complex with phase = angle of reference
            magTmp = ccfmap;
            if ( tmpDecoy )
              decoyTmp = decoy;
            end
            angTmp = zeros(size(magTmp), 'single','gpuArray');
            angTmp = angTmp + 1;

            firstLoopOverTomo  = false;
            firstLoopOverChunk = false;
            
            intraLoopAngle = intraLoopAngle + 1;
            currentGlobalAngle = currentGlobalAngle + 1;

        elseif (firstLoopOverChunk)
          % These double cuts are old, and don't really make sense. Make
          % this more consistant with current operations when there is
          % time.
          magTmp =  RESULTS_peak(iCut(1):iCut(1)+sizeChunk(1)-1,...
                                 iCut(2):iCut(2)+sizeChunk(2)-1,...
                                 iCut(3):iCut(3)+sizeChunk(3)-1);
          angTmp = RESULTS_angle(iCut(1):iCut(1)+sizeChunk(1)-1,...
                                 iCut(2):iCut(2)+sizeChunk(2)-1,...
                                 iCut(3):iCut(3)+sizeChunk(3)-1);
          if ( tmpDecoy )
            decoyTmp = RESULTS_decoy(iCut(1):iCut(1)+sizeChunk(1)-1,...
                                     iCut(2):iCut(2)+sizeChunk(2)-1,...
                                     iCut(3):iCut(3)+sizeChunk(3)-1);   
            decoyTmp = gpuArray(decoyTmp(vA(1,1) + 1:end - vA(2,1), ...
                                         vA(1,2) + 1:end - vA(2,2), ...
                                         vA(1,3) + 1:end - vA(2,3))); 
            decoyTmp(decoyTmp < decoy) = decoy(decoyTmp < decoy);                                       
          end

          magTmp = gpuArray(magTmp(vA(1,1) + 1:end - vA(2,1), ...
                                   vA(1,2) + 1:end - vA(2,2), ...
                                   vA(1,3) + 1:end - vA(2,3)));    
          angTmp = gpuArray(angTmp(vA(1,1) + 1:end - vA(2,1), ...
                                   vA(1,2) + 1:end - vA(2,2), ...
                                   vA(1,3) + 1:end - vA(2,3)));
          
          firstLoopOverChunk = false;

          replaceTmp = ( magTmp < ccfmap );
          
          magTmp(replaceTmp) = ccfmap(replaceTmp);
          angTmp(replaceTmp) = currentGlobalAngle;
          
          

          intraLoopAngle = intraLoopAngle + 1;
          currentGlobalAngle = currentGlobalAngle + 1;
          clear replaceTmp

        else
            % update higher values of ccfmap with new reference if applicable.


            replaceTmp = ( magTmp < ccfmap );


            magTmp(replaceTmp) = ccfmap(replaceTmp);
            angTmp(replaceTmp) = currentGlobalAngle;
            if ( tmpDecoy )
              decoyTmp(decoyTmp < decoy) = decoy(decoyTmp < decoy);
            end
            intraLoopAngle = intraLoopAngle + 1;
            currentGlobalAngle = currentGlobalAngle + 1;
            clear replaceTmp
        end
      nComplete = nComplete + 1;
      end
    end

    % After searching all angles on this chunk, but out meaningful
    % portion for storage.

    magStoreTmp =  RESULTS_peak(iCut(1):iCut(1)+sizeChunk(1)-1,...
                                iCut(2):iCut(2)+sizeChunk(2)-1,...
                                iCut(3):iCut(3)+sizeChunk(3)-1);
    angStoreTmp = RESULTS_angle(iCut(1):iCut(1)+sizeChunk(1)-1,...
                                iCut(2):iCut(2)+sizeChunk(2)-1,...
                                iCut(3):iCut(3)+sizeChunk(3)-1);


    magStoreTmp(vA(1,1) + 1:end - vA(2,1), ...
                vA(1,2) + 1:end - vA(2,2), ...
                vA(1,3) + 1:end - vA(2,3)) = gather(magTmp);
    angStoreTmp(vA(1,1) + 1:end - vA(2,1), ...
                vA(1,2) + 1:end - vA(2,2), ...
                vA(1,3) + 1:end - vA(2,3)) = gather(angTmp);


     RESULTS_peak(iCut(1):iCut(1)+sizeChunk(1)-1,...
                  iCut(2):iCut(2)+sizeChunk(2)-1,...
                  iCut(3):iCut(3)+sizeChunk(3)-1) = magStoreTmp;
                
     clear magStoreTmp

    RESULTS_angle(iCut(1):iCut(1)+sizeChunk(1)-1,...
                  iCut(2):iCut(2)+sizeChunk(2)-1,...
                  iCut(3):iCut(3)+sizeChunk(3)-1) = angStoreTmp;
     clear angStoreTmp
    
    if ( tmpDecoy )
    decoyStoreTmp =  RESULTS_decoy(iCut(1):iCut(1)+sizeChunk(1)-1,...
                                iCut(2):iCut(2)+sizeChunk(2)-1,...
                                iCut(3):iCut(3)+sizeChunk(3)-1);
    decoyStoreTmp(vA(1,1) + 1:end - vA(2,1), ...
                vA(1,2) + 1:end - vA(2,2), ...
                vA(1,3) + 1:end - vA(2,3)) = gather(decoyTmp);   
     RESULTS_decoy(iCut(1):iCut(1)+sizeChunk(1)-1,...
                  iCut(2):iCut(2)+sizeChunk(2)-1,...
                  iCut(3):iCut(3)+sizeChunk(3)-1) = decoyStoreTmp;              
      
    end
    tomoTime = toc;
    totalTime = totalTime + toc; timeEstimate = totalTime * (nTomograms*nAngles(1)./(nComplete-1));
    fprintf('elapsed time = %f s  est remain %f s\n', tomoTime, timeEstimate);
    tomoIDX = tomoIDX + 1;
    firstLoopOverAngle = false;
    currentGlobalAngle = currentGlobalAngle - intraLoopAngle + 1;
  end
  
  currentGlobalAngle = currentGlobalAngle + intraLoopAngle -  1;
end
%save('angle_list.txt','angle_list','-ascii');
clear tomoStack
% Cut out the post padding used to iterate over the tomogram
RESULTS_peak = RESULTS_peak(1+tomoPre(1):end-tomoPost(1),...
                            1+tomoPre(2):end-tomoPost(2),...
                            1+tomoPre(3):end-tomoPost(3));
%RESULTS_peak(RESULTS_peak < 0) = 0;                
RESULTS_angle = RESULTS_angle(1+tomoPre(1):end-tomoPost(1),...
                              1+tomoPre(2):end-tomoPost(2),...
                              1+tomoPre(3):end-tomoPost(3));

if ( tmpDecoy )
  RESULTS_decoy = RESULTS_decoy(1+tomoPre(1):end-tomoPost(1),...
                            1+tomoPre(2):end-tomoPost(2),...
                            1+tomoPre(3):end-tomoPost(3));
%   RESULTS_decoy = RESULTS_decoy ./ std(RESULTS_decoy(:));
  RESULTS_decoy(RESULTS_decoy < 1) = 1; 
  
end
gpuDevice(useGPU);
clear bhF


% scale the magnitude of the results to be 0 : 1
szK = latticeRadius;%floor(0.8.*szM);
rmDim = max(max(eraseMaskRadius),max(szK)).*[1,1,1];
mag = RESULTS_peak; clear RESULTS_peak
% Normalize so the difference if using a decoy makes sense. The input decoy
% should have the same power, so I'm not sure why this is needed, but it is
% an easy fix and a problem for future Ben to figure out.
% mag = mag ./ std(mag(:));

system(sprintf('mkdir -p %s',convTMPNAME));
system(sprintf('mv temp_%s.mrc %s',convTMPNAME,convTMPNAME));

resultsOUT = sprintf('./%s/%s_convmap.mrc',convTMPNAME,mapName);
anglesOUT  = sprintf('./%s/%s_angles.mrc',convTMPNAME,mapName);
angleListOUT = sprintf('./%s/%s_angles.list',convTMPNAME,mapName);
SAVE_IMG(MRCImage(mag),resultsOUT);
SAVE_IMG(MRCImage(RESULTS_angle),anglesOUT);
if ( tmpDecoy )
  decoyOUT = sprintf('./%s/%s_decoy.mrc',convTMPNAME,mapName);
  SAVE_IMG(MRCImage((RESULTS_decoy)),decoyOUT);
  diffOUT = sprintf('./%s/%s_convmap-decoy.mrc',convTMPNAME,mapName);
  decoyLogical = mag < RESULTS_decoy;
  mag(decoyLogical) = 0;
  mag(~decoyLogical) = mag(~decoyLogical) - RESULTS_decoy(~decoyLogical); clear RESULTS_decoy
  SAVE_IMG(MRCImage((mag)),diffOUT);
end
angleFILE = fopen(angleListOUT,'w');
fprintf(angleFILE,'%2.2f\t%2.2f\t%2.2f\n', ANGLE_LIST');
fclose(angleFILE);


% mag =  mag - min(mag(:)); mag = mag ./ max(mag(:));

% Zero out one lattice width from the edges to reduce edge effect (but cutting
% out and padding back in.) Also pad by size of removal mask (subtract this from
% coordinates)
mag = mag(szK(1)+1:end - szK(1), ...
          szK(2)+1:end - szK(2), ...
          szK(3)+1:end - szK(3));
mag = BH_padZeros3d(mag,szK+rmDim,szK+rmDim, 'cpu', 'single');
%dev.FreeMemory;
%%%Ang = angle(RESULTS_peak); %clear Results
% negative phase angles mapped back to 0-->pi
%Ang(sign(Ang) < 0) = Ang(sign(Ang)<0) + pi; serotonin_ali1_75_1.mod
%Ang = BH_padZeros3d(round(Ang./angleIncrement),szK,szK,'cpu','single');
Ang = BH_padZeros3d(RESULTS_angle,rmDim,rmDim,'cpu','single');


%mag =  mag - min(mag(:)); mag = mag ./ max(mag(:));

Tmean = mean(mag(( mag ~= 0 )));
Tstd  = std(mag(( mag~=0 )));
threshold = Tmean + peakThreshold*Tstd;
mag((Ang < 0)) = 0;

mag = gpuArray(mag);
sizeTomo = size(mag);


[MAX, coord] = max(mag(:));

peakMat = zeros(peakThreshold,10*nPeaks);

n = 1;

fprintf('rmDim %f szK %f\n',  rmDim,szK);
removalMask = BH_mask3d(eraseMaskType,[2,2,2].*rmDim+1,eraseMaskRadius,[0,0,0]);

maskCutOff = 0.999;
nIncluded = gather(sum(sum(sum(removalMask > maskCutOff))));
nTries = 0;
if strcmpi(eraseMaskType,'rectangle');
  areaPreFactor = 0;
else
  areaPreFactor = (4/3*pi);
end

while nIncluded < areaPreFactor*prod(eraseMaskRadius)
  maskCutOff = 0.99*maskCutOff;
  nIncluded = gather(sum(sum(sum(removalMask > maskCutOff))));
  nTries = nTries + 1;
  if (nTries > 1000)
    error('Did not find an appropriate erase mask');
  end

end

if ignore_threshold
  highThr = 0;
end

this_try = 0;
while n <= 2.*peakThreshold && (this_try < max_tries) && MAX > highThr
this_try = this_try + 1;

%
% Some indicies come back as an error, even when they seem like the
% should be fine. I'm not sure why, and I should think about this
% more, but for now, just set that one index to zero (instead of a
% whole box) and move on with life. It looks like the index that is
% kicking out the error is equal to -1*numberofreferences, which
% might be an issue because that corresonds to the positive upper
% limit of the reference index. Ignoring it still seems to be okay
% but it bothers me not to know.
        
        
[i,j,k] = ind2sub(sizeTomo,coord);
try
  c = gather([i,j,k]);
catch
  print('Ran into some trouble gathering the i,j,k. Breaking out\n');
  break
end

  if Ang(gather(coord)) > 0

    % box for removal and center of mass calc, use a larger box if multiple
    % peaks are being saved.
    bDist = 1+round(log(nPeaks));
    clI  = c(1) - bDist;
    chI  = c(1) + bDist;
    clJ  = c(2) - bDist;
    chJ  = c(2) + bDist;
    clK  = c(3) - bDist;
    chK  = c(3) + bDist;

    magBox = mag(clI:chI,clJ:chJ,clK:chK);
    
    angBox = Ang(clI:chI,clJ:chJ,clK:chK);

    [cmX, cmY, cmZ] = ndgrid(-1*bDist:1*bDist, ...
                             -1*bDist:1*bDist, ...
                             -1*bDist:1*bDist );

    cMass = [ sum(sum(sum(magBox.*cmX))) ; ... 
              sum(sum(sum(magBox.*cmY))) ; ...
              sum(sum(sum(magBox.*cmZ))) ] ./ sum(magBox(:));



%     cenP = [ (c(1)+cMass(1)-1) - sizeTomo(1)./2 ,...
%              (c(2)+cMass(2)-1) - sizeTomo(2)./2 ,...
%              (c(3)+cMass(3)-1) - sizeTomo(3)./2 ];

    % Switching from centered to lower left coordinates and subtracting the
    % padding 
    
    cenP = c + cMass' - rmDim;

    

    % If the most frequent peak is unique use it;
    [peakM, ~, peakC] = mode(angBox(:));
    if length(peakC) == 1 && peakM
      % Need to ensure the mode is none zero which is possible.
      peakMat(n,4:6) = ANGLE_LIST(peakM,:);
      topPeak = peakM;
    else
      % Otherwise use the value at the max for the peak val;
      peakMat(n,4:6) = ANGLE_LIST(Ang(coord),:);
      topPeak = Ang(coord);
    end
    peakMat(n,1:3) = gather(samplingRate.*cenP);
    
    if nPeaks > 1
      oldPeaks = ( angBox == topPeak );
      
      for iPeak = 2:nPeaks
        [peakM, ~, ~] = mode(angBox(~oldPeaks));
        % There could be redundancy, as given by peakC, but just take the
        % first value given by peak M.
        peakMat(n,[1:3]+10*(iPeak-1)) = gather(samplingRate.*cenP);
        peakMat(n,[4:6]+10*(iPeak-1)) = ANGLE_LIST(peakM,:);
        
         oldPeaks = ( angBox == peakM | oldPeaks );
        
      end
      
    end

      

   
    rmMask = BH_resample3d(removalMask,peakMat(n,4:6),[0,0,0],'Bah','GPU','forward');
    % Invert after resampling so that zeros introduced by not extrapolating
    % the corners are swapped to ones, i.e. not removed.
%     rmMask = (1-rmMask);
    
    mag(c(1)-rmDim:c(1)+rmDim,...
        c(2)-rmDim:c(2)+rmDim,...
        c(3)-rmDim:c(3)+rmDim) = ...
    mag(c(1)-rmDim:c(1)+rmDim,...
        c(2)-rmDim:c(2)+rmDim,...
        c(3)-rmDim:c(3)+rmDim) .* (rmMask< maskCutOff);

    peakMat(n,10) = (gather(MAX) - Tmean)./Tstd; % record stds above mean
    n = n + 1;
    
    if ~mod(n,100)
      n
    end

  else
       Ang(gather(coord));
       mag(coord) = 0;
  end


[MAX, coord] = max(mag(:));

end

peakMat = peakMat( ( peakMat(:,1)>0 ),:);

%save('peakMat_post.mat', 'peakMat');

% A temp test, not the correct output just score x y z dx dy dz e1 e2 e3

csv_out = sprintf('./%s/%s.csv',convTMPNAME,mapName);
pos_out = sprintf('./%s/%s.pos',convTMPNAME,mapName);
%fieldOUT = zeros(length(peakMat(:,1)),26);
fileID = fopen(csv_out,'w');
fileID2 = fopen(pos_out,'w');
errID  = fopen(sprintf('./%s/%s.errID',convTMPNAME,mapName));


if SYMMETRY > 1
  symmetry = 0:360/SYMMETRY:359;
  symCell = cell(length(symmetry),1);
  for iSym = 1:length(symmetry)
    symCell{iSym} = BH_defineMatrix([symmetry(iSym),0,0], 'Bah', 'inv');
  end

end
n=1
for i = 1:length(peakMat(:,1))
   if all(peakMat(i,1:3))
        
        if SYMMETRY > 1
          % Generate a uniform distribution over the in-plane
          % randomizations
          iSym = rem( n + SYMMETRY, SYMMETRY)+1;
          r = reshape(BH_defineMatrix(peakMat(i,4:6), 'Bah', 'inv')*...
                      symCell{iSym},1,9);
        else
          r = reshape(BH_defineMatrix(peakMat(i,4:6), 'Bah', 'inv'),1,9);
        end
        fprintf(fileID,['%1.2f %d %d %d %d %d %d %d %d %d %f %f %f %d %d %d ',...
                        '%f %f %f %f %f %f %f %f %f %d '],peakMat(i,10),samplingRate,0, ...
                        i+nPreviousSubTomos,1,1,1,1,1,0,peakMat(i,1:3), ...
                        peakMat(i,4:6),r,1);
         
        if nPeaks > 1
          
          for iPeak = 2:nPeaks
            if SYMMETRY > 1
              % Generate a uniform distribution over the in-plane
              % randomizations
              iSym = rem( n + SYMMETRY, SYMMETRY)+1;
              r = reshape(BH_defineMatrix(peakMat(i,[4:6]+10*(iPeak-1)), 'Bah', 'inv')*...
                          symCell{iSym},1,9);
            else
              r = reshape(BH_defineMatrix(peakMat(i,[4:6]+10*(iPeak-1)), 'Bah', 'inv'),1,9);
            end            
            fprintf(fileID,['%1.2f %d %d %d %d %d %d %d %d %d %f %f %f %d %d %d ',...
                        '%f %f %f %f %f %f %f %f %f %d '],peakMat(i,10),samplingRate,0, ...
                        i+nPreviousSubTomos,1,1,1,1,1,0,peakMat(i,[1:3]+10*(iPeak-1)), ...
                        peakMat(i,[4:6]+10*(iPeak-1)),r,1);          
          end
          
          
        end
                      
         fprintf(fileID,'\n');           
                      
                      
                      
                      
        fprintf(fileID2,'%f %f %f\n',peakMat(i,1:3)./samplingRate);                      

        
        n = n +1;
   end
end

%lastIndex = find(fieldOUT(:,4),1,'last');

fclose(fileID);
fclose(fileID2);

system(sprintf('point2model -number 1 -sphere 3 -scat ./%s/%s.pos ./%s/%s.mod', convTMPNAME,mapName,convTMPNAME, mapName));

fileID = fopen(sprintf('./%s/%s.path',convTMPNAME,mapName),'w');
fprintf(fileID,'%s,%s,%s,%s',mapName,mapPath,mapExt,RAWTLT);
fclose(fileID);
%subTomoMeta.('cycle000').('geometry').(mapName) = fieldOUT;
% subTomoMeta.('mapPath').(mapName) = mapPath;
% subTomoMeta.('mapExt').(mapName) = mapExt;

% if any(ismember(fieldnames(subTomoMeta), 'nSubTomoTotal'))
%   subTomoMeta.('nSubTomoTotal') = subTomoMeta.('nSubTomoTotal') + lastIndex;
% else
%   subTomoMeta.('nSubTomoTotal') = lastIndex;
% end

% preFscSplit = gather(subTomoMeta);
% 
% % Randomly divide the data into half sets.
% [ subTomoMeta ] = BH_fscSplit( preFscSplit );
% subTomoMeta.('currentCycle') = 0;

% save(sprintf('%s.mat', pBH.('subTomoMeta')), 'subTomoMeta');
% save(sprintf('./convmap/%s.mat~', pBH.('subTomoMeta')), 'subTomoMeta');
%save('test.pos','a','-ascii');


fprintf('Total execution time : %f seconds\n', etime(clock, startTime));
   


end % end of templateSearch3d function

