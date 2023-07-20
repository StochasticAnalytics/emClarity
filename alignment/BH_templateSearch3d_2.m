function []  = BH_templateSearch3d_2( PARAMETER_FILE,...
                                        tomoName,tomoNumber,TEMPLATE, ...
                                        SYMMETRY, wedgeType, varargin)
 
                                                       
%3d template matching

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



ctf3dNoSubTomoMeta = true;
if length(varargin) == 1
  % Allow for an override of the max number, useful when only a few tomos
  % have a strong feature like carbon that is hard to avoid.
  gpuIDX = EMC_str2double(varargin{1});
elseif length(varargin) > 1
  error('emClarity templateSearch paramN.m tiltN regionN referenceName symmetry(C1) <optional gpuIDX>');
end

  tomoNumber = EMC_str2double(tomoNumber);


  [ useGPU ] = BH_multi_checkGPU( gpuIDX )



gpuDevice(useGPU);

%  For now just override this as it doesn't do too much (randomizing symmetry mates. And doesn't work with new symmetry ops
% Need to get the symmetry mats from an interpolator object
% SYMMETRY = EMC_str2double(SYMMETRY);
SYMMETRY=1;

startTime = clock ;

pBH = BH_parseParameterFile(PARAMETER_FILE);

if ctf3dNoSubTomoMeta
   mapBackIter = 0;
  shouldBeCTF = 1;
else
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
    if (EMC_str2double(v{1}) < 4 || (EMC_str2double(v{2}) <= 10 && EMC_str2double(v{3}) < 42))
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

 peakThreshold = pBH.('Tmp_threshold');


latticeRadius = pBH.('particleRadius');
try
  targetSize    = pBH.('Tmp_targetSize')
catch
  targetSize = [512,512,512];
end
angleSearch   = pBH.('Tmp_angleSearch');


convTMPNAME = sprintf('convmap_wedgeType_%d_bin%d',wedgeType,samplingRate)

try 
  use_new_grid_search = pBH.('use_new_grid_search');
catch
  use_new_grid_search = true;
end

try
  symmetry = pBH.('symmetry');
catch
  error('You must now specify a symmetry=X parameter, where symmetry E (C1,C2..CX,O,I)');
end

try 
  eraseMaskType = pBH.('Peak_mType');
catch
  eraseMaskType = 'sphere';
end
try
  eraseMaskRadius = pBH.('Peak_mRadius');
catch
  eraseMaskRadius = 1.0.*latticeRadius;
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
catch
  max_tries = 10000;
end

try
  over_ride =  pBH.('Override_threshold_and_return_N_peaks')
  ignore_threshold = true;
  fprintf('Override_threshold_and_return_N_peaks set to true, returning exactly %d peaks\n', over_ride);
  peakThreshold = over_ride;
end

pixelSizeFULL = pBH.('PIXEL_SIZE').*10^10;
if pBH.('SuperResolution')
  pixelSizeFULL = pixelSizeFULL * 2;
end

pixelSize = pixelSizeFULL.*samplingRate;

% For testing
print_warning=false;
try 
  wantedCut = pBH.('lowResCut');
  fprintf('lowResCut is deprecated and will be removed in future versions.\n')
  fprintf('please switch to Tmp_bandpass\n\n');
  bp_vals = [1e-3,600,wantedCut];
  print_warning = true; 
catch
  bp_vals = [1e-3,600,28];
end

try
    bp_vals = pBH.('Tmp_bandpass');
    if numel(bp_vals) ~= 3
        error('Tmp_bandpass is [filter at zero freq, res high-pass cutoff, res low-pass cutoff]');
    end
    if print_warning
        fprintf('WARNING, you specified lowResCut (deprecated) and Tmp_bandpass!\n');
    end
    fprintf('You specified a bandpass with values [%2.2e,%3.2f,%3.2f]\n',bp_vals);
catch
    bp_vals = [1e-3,600,28];
    fprintf('Using default bandpass with values [%2.2e,%3.2f,%3.2f]\n',bp_vals);

end
try
    stats_diameter_fraction = pBH.('diameter_fraction_for_local_stats')
catch
    stats_diameter_fraction = 1
end

sum_of_x =  [];
sum_of_x2 = [];
try 
  rescale_mip = pBH.('rescale_mip');
catch
  rescale_mip = false;
end

% Limit to the first zero if we are NOT using the CTF rec
if (shouldBeCTF ~= 1)
TLT = load(sprintf('fixedStacks/ctf/%s_ali%d_ctf.tlt',tomoName,mapBackIter+1));
  def = mean(-1.*TLT(:,15))*10^6; %TODO if you switch to POSITIVEDEFOCUS this will be wrong
  firstZero = -0.2*def^2 +5.2*def +11;

  % Take the lower of firstZero lowResCut or Nyquist
  bp_vals(3) = max(bp_vals(3), firstZero);
  fprintf('\nUsing max (%f) of specified resolution cutoff of %f and first ctf zero %f Angstrom\n',bp_vals(3), wantedCut, firstZero);

end

if pixelSize*2 >  bp_vals(3)
  fprintf('\nLimiting to Nyquist (%f) instead of user requested low pass cutoff %f Angstrom\n',pixelSize*2,bp_vals(3));
  bp_vals(3) = pixelSize*2;
end


mapPath = './cache';
mapName = sprintf('%s_%d_bin%d',tomoName,tomoNumber,samplingRate);
mapExt = '.rec';

sprintf('recon/%s_recon.coords',tomoName)
[ recGeom, ~, ~] = BH_multi_recGeom( sprintf('recon/%s_recon.coords',tomoName) );

reconCoords = recGeom(tomoNumber,:);
clear recGeom


bp_vals(2) = 2.*max(latticeRadius);
statsRadiusAng = stats_diameter_fraction.*[2,2,2].*max(latticeRadius);
statsRadius = ceil(statsRadiusAng./pixelSize); % Convert to binned pixels
maskRadius  = ceil(0.5.*[1,1,1].*max(latticeRadius)./pixelSize);
latticeRadius = (0.75 .* latticeRadius) ./ (pixelSize);
latticeRadius = floor(latticeRadius);
latticeRadius = latticeRadius + mod(latticeRadius, 2);

eraseMaskRadius = floor((eraseMaskRadius) ./ (pixelSize));
eraseMaskRadius = eraseMaskRadius + mod(eraseMaskRadius,2);

fprintf('EXPERIMENTAL setting the highpass to match the max particle diameter. %3.3f Ang\n\n', bp_vals(2));

fprintf('\ntomograms normalized in %f Angstrom cubic window\n',statsRadiusAng(1));

fprintf('\nlatticeRadius = %dx%dx%d pixels\n\n', latticeRadius);
fprintf('\neraseMaskType %s, eraseMaskRadius %dx%dx%d pixels\n',eraseMaskType,eraseMaskRadius);
  % For wedgeMask
particleThickness =  latticeRadius(3);

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
                            
% Bandpass the template so it is properly normalized
bp_vals
temp_bp = BH_bandpass3d(size(template),bp_vals(1),0.3.*bp_vals(2),bp_vals(3),'GPU',pixelSizeFULL);
template = real(ifftn(fftn(gpuArray(template)).*temp_bp.^2));
clear temp_bp

                            
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


templateBIN = BH_reScale3d(gather(template),'',sprintf('%f',1/samplingRate),'cpu');

templateBIN = templateBIN - mean(templateBIN(:));
templateBIN = templateBIN  ./rms(templateBIN(:));

[templateMask] = gather(EMC_maskReference(gpuArray(templateBIN),pixelSize,{'fsc', true}));


sizeTemp = size(template);
sizeTempBIN = size(templateBIN);




gpuDevice(useGPU);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize a whole mess of control variables and storage volumes. %
%Out of plane range inc (starts from 1.* inc)
rotConvention = 'Bah';
% Check and override the rotational convention to get helical averaging.
% Replaces the former hack of adding a fifth dummy value to the angular search
try
  doHelical = pBH.('doHelical');
catch
  doHelical = 0;
end
if ( doHelical )
  rotConvention = 'Helical'
end

rotConvention

if (use_new_grid_search)

  gridSearch = eulerSearch(symmetry, angleSearch(1),...
        angleSearch(2),angleSearch(3),angleSearch(4), 0, 0, false);
  nAngles = sum(gridSearch.number_of_angles_at_each_theta);
  inPlaneSearch = gridSearch.parameter_map.psi;
  
  
else

  [  nInPlane, inPlaneSearch, angleStep, nAngles] ...
                                      = BH_multi_gridSearchAngles(angleSearch)
end

                                                                       

highThr=sqrt(2).*erfcinv(ceil(peakThreshold.*0.10).*2./(prod(size(tomogram)).*nAngles(1)))

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

RMSFACTOR = sqrt(prod(sizeTempBIN) / prod(sizeChunk));

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

% [ tomogram ] = BH_padZeros3d(tomogram, tomoPre, tomoPost, ...
%                                              'cpu', 'singleTaper',mean(tomogram(:)));
tomogram = padarray(tomogram,tomoPre,'symmetric','pre');
tomogram = padarray(tomogram,tomoPost,'symmetric','post');
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
if (rescale_mip)
  sum_of_x = zeros(sizeTomo, 'single'); 
  sum_of_x2 = zeros(sizeTomo, 'single'); 
end
% Loop over tomogram
% Set this up second


% % % % optimize fft incase a power of two is not used, this will make things run ok.
% % % opt = zeros(sizeChunk, precision,'gpuArray');
% % % fftw('planner','patient');
% % % fftn(opt);
% % % clear opt ans
[ bhF ] = fourierTransformer(randn(sizeChunk, 'single','gpuArray'));

sum_template = mean(templateBIN(:));
sum_templateMask = mean(templateMask(:));
sum_imgMask = prod(sizeChunk);% bhF.halfDimSize * sizeChunk(2) * sizeChunk(3);


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


try
  doMedFilt = pBH.('Tmp_medianFilter');
  if ~ismember(doMedFilt,[3,5,7])
    error('Tmp_medianFilter can only be 3,5, or 7');
  else
    fprintf('Using median filter, size %d',doMedFilt);
  end
catch
  doMedFilt = 0;
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
    
    tomoChunk = tomoChunk - mean(tomoChunk(:));
    tomoChunk = tomoChunk ./ rms(tomoChunk(:));

%       tomoChunk = real(ifftn(fftn(tomoChunk).*tomoBandpass));
    
    tomoChunk = bhF.invFFT(bhF.fwdFFT(tomoChunk,0,0,[bp_vals, pixelSize]),2);
    

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

    [ averageMask, flgOOM ] = BH_movingAverage_2(tomoChunk, statsRadius(1)); 
    rmsMask =  BH_movingAverage_2(tomoChunk.^2, statsRadius(1)); 
    rmsMask = sqrt(rmsMask - averageMask.^2);
    tomoChunk = (tomoChunk - averageMask) ./ rmsMask;
    clear rmsMask averageMask
%     averageMask = gather(averageMask);

%         [ rmsMask ] = gather(BH_movingRMS_3(tomoChunk, statsRadius(1), averageMask)); 
        

%     tomoChunk = tomoChunk - averageMask; 
% tomoChunk = tomoChunk ./ rmsMask;
%     if (save_average_filtered)
%               tomoChunk = gpuArray(tomogram(cutX:cutX+sizeChunk(1)-1,...
%                                   cutY:cutY+sizeChunk(2)-1,...
%                                   cutZ:cutZ+sizeChunk(3)-1));
%       avgFiltRec = zeros(size(tomogram),'single');
%     end        
% statsRadius(1)
%     [ rmsMask ] = BH_movingRMS_2(tomoChunk-averageMask, statsRadius(1));
%     statsRadius(1)
%     tomoChunk = tomoChunk ./ rmsMask;
% figure, imshow3D(BH_padZeros3d(gather(averageMask),'fwd',trimValid,'cpu','single'))
% figure, imshow3D(BH_padZeros3d(gather(rmsMask),'fwd',trimValid,'cpu','single'))
% return
%     clear avgerageMask 
 

    tomoChunk = gather(((-1*shouldBeCTF) .* tomoChunk )).*validAreaMask;   
 
%     tomoChunk = tomoChunk .* (-1*shouldBeCTF); % This is backwards, but I don't know why
    tmp_sum = sum(tomoChunk(validAreaMask > 0.1));
    fullX = fullX + tmp_sum;
    fullX2 = fullX2 + gather(tmp_sum.^2);
    fullnX = fullnX + gather(prod(sizeChunk));
    
    tomoStack(:,:,:,tomoIDX) = tomoChunk;

    tomoCoords(tomoIDX,:) = [cutX,cutY,cutZ];
    tomoIDX = tomoIDX + 1;

      
    end % end of loop over Z chunks
  end % end of loop over Y chunks
end % end of loop over X chunks

% Normalize the global variance
globalVariance = (fullX2/fullnX) - (fullX/fullnX)^2;
%fprintf('After local normalization, scaling also the global variance %3.3e\n',globalVariance);

for iChunk = 1:tomoIDX-1
  tomoStack(:,:,:,iChunk) = tomoStack(:,:,:,iChunk) ./ sqrt(globalVariance);
end


  
clear tomoWedgeMask   bandpassFilter statBinary validAreaMask tomoChunk

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

kVal = 0;


    

currentGlobalAngle = 1;
ANGLE_LIST = zeros(nAngles(1),3, 'single');
nComplete = 0;
totalTime = 0;
firstLoopOverTomo = true;

    % Center the spectrum by multiplication not swapping (this should just
    % be in the fourierTransformer class if it isn't already)
    % swapPhase(obj, inputVol, direction) with fwd should do it
    [dU,dV,dW] = BH_multi_gridCoordinates(size(tomoStack(:,:,:,1)),...
                                         'Cartesian','GPU', ...
                                          {'none'},1,1,0);
                                        
    swapQuadrants = exp((-2i*pi).*(dU.*(floor(size(dU,1)/2)+1) + ...
                                  (dV.*(floor(size(dV,2)/2)+1) + ...          
                                  (dW.*(floor(size(dW,3)/2)+1)))));  
    clear dU dV dW
    
    swapQuadrants = swapQuadrants(1:floor(size(swapQuadrants,1)/2)+1,:,:);

if (use_new_grid_search)
  theta_search = 1:gridSearch.number_of_out_of_plane_angles;
else
  theta_search = 1:size(angleStep,1);
end


for iAngle = theta_search
  
  if (use_new_grid_search)
    theta = gridSearch.parameter_map.theta(iAngle);
    numRefIter = gridSearch.number_of_angles_at_each_theta(iAngle);
  else
    theta = angleStep(iAngle,1);  
    phiStep = angleStep(iAngle,3);
    numRefIter = angleStep(iAngle,2)*length(inPlaneSearch)+1;
  end

  

  tempImg = gpuArray(templateBIN); %%%%% NEW switch to bin

  
  interpolationNormFactor = sum(abs(tempImg(:)).^2);


  clear referenceStack tempFilter
  % Calculate all references for each out of plane tilt only once
  referenceStack = zeros([sizeTempBIN,numRefIter], 'single', 'gpuArray');
                                                         
  tomoIDX = 1;
  firstLoopOverAngle = true;
  
  % Avoid repeated allocations
  tempPAD = zeros(size(tempImg) + padBIN(1,:) + padBIN(2,:),'single','gpuArray');
  tempPADMask = tempPAD;
  
  template_interpolator = '';
  [template_interpolator, ~] = interpolator(tempImg,[0,0,0],[0,0,0], 'Bah', 'forward', 'C1', false);
  
  templateMask_interpolator = '';
  [templateMask_interpolator, ~] = interpolator(gpuArray(templateMask),[0,0,0],[0,0,0], 'Bah', 'forward', 'C1', false);

  


  
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
    
    if (use_new_grid_search)
      fprintf('Working on tilt(%d/%d) tomoChunk(idx%d/%d)\t' ...
                            ,iAngle,gridSearch.number_of_out_of_plane_angles, tomoIDX,nTomograms);      
    else
      fprintf('working on tilt(%d/%d) tomoChunk(idx%d/%d)\t' ...
                            ,iAngle,size(angleStep,1), tomoIDX,nTomograms);
    end


      tomoFou = gpuArray(tomoStack(:,:,:,tomoIDX));
      tomoFou = swapQuadrants.*bhF.fwdFFT(tomoFou);





    if (use_new_grid_search)
      phi_search = gridSearch.parameter_map.phi{iAngle};
    else
      phi_search = 0:angleStep(iAngle,2);
    end
    
    for iAzimuth = phi_search

      if (use_new_grid_search)
        phi = iAzimuth;
      else
        phi = phiStep * iAzimuth;    
      end

      for iInPlane = inPlaneSearch
        psi = iInPlane;

        %calc references only on first chunk
        if (firstLoopOverAngle)

          ANGLE_LIST(currentGlobalAngle,:) = [phi, theta, psi - phi];

        end
                              
        [ tempRot ] = template_interpolator.interp3d(...
                                                   [phi, theta, psi - phi],... 
                                                   [1,1,1],rotConvention,...
                                                  'forward','C1');  

         

          tempPAD = tempPAD .* 0;
          tempPAD(padBIN(1,1)+1: end - padBIN(2,1), ...
                  padBIN(1,2)+1: end - padBIN(2,2), ...
                  padBIN(1,3)+1: end - padBIN(2,3)) = tempRot;
                
          tempFou = conj(bhF.fwdFFT(tempPAD));          


        ccfmap = BH_padZeros3d((real(single(...
                               bhF.invFFT(tomoFou.*tempFou)))),...%./(tomoNorm.*tempNorm))))),...
                               trimValid(1,:),trimValid(2,:),'GPU','single');
%                              

        ccfmap = ccfmap ./ std(ccfmap(:));
        
        if ( tmpDecoy > 0 )
          
           if (firstLoopOverAngle)

             decoy = BH_padZeros3d(BH_reScale3d(tempRot./decoyNorm,'',tmpDecoy,'GPU',decoyShift),...
                                   padDecoy(1,:),padDecoy(2,:),'GPU','single');
           else
              % Probably just make a second decoy stack to avoid
              % re-interpolating. If it works, then do this.
              error('This is temp broken with new interpolator');
%              decoy = BH_padZeros3d(BH_reScale3d(referenceStack(:,:,:,intraLoopAngle)./decoyNorm,'',tmpDecoy,'GPU',decoyShift),...
%                                    padDecoy(1,:),padDecoy(2,:),'GPU','single');                              
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
            angTmp = ones(size(magTmp), 'single','gpuArray');

            if (rescale_mip)
              sum_of_x_tmp = ccfmap; 
              sum_of_x2_tmp = ccfmap.^2; 
            end
            
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
          
          if (rescale_mip)
              
            sum_of_x_tmp = sum_of_x(iCut(1):iCut(1)+sizeChunk(1)-1,...
                                     iCut(2):iCut(2)+sizeChunk(2)-1,...
                                     iCut(3):iCut(3)+sizeChunk(3)-1);   
            sum_of_x_tmp = gpuArray(sum_of_x_tmp(vA(1,1) + 1:end - vA(2,1), ...
                                         vA(1,2) + 1:end - vA(2,2), ...
                                         vA(1,3) + 1:end - vA(2,3))); 
            sum_of_x_tmp = sum_of_x_tmp + ccfmap; 
            
            sum_of_x2_tmp = sum_of_x2(iCut(1):iCut(1)+sizeChunk(1)-1,...
                                     iCut(2):iCut(2)+sizeChunk(2)-1,...
                                     iCut(3):iCut(3)+sizeChunk(3)-1);   
            sum_of_x2_tmp = gpuArray(sum_of_x2_tmp(vA(1,1) + 1:end - vA(2,1), ...
                                         vA(1,2) + 1:end - vA(2,2), ...
                                         vA(1,3) + 1:end - vA(2,3))); 
            sum_of_x2_tmp = sum_of_x2_tmp + ccfmap.^2;              
          end         

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
            
            if (rescale_mip)
              sum_of_x_tmp = sum_of_x_tmp + ccfmap;              
              sum_of_x2_tmp = sum_of_x2_tmp + ccfmap.^2;              
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

    
    % FIXME this double cutting and temporary allocation is ridiculous.
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
    
    if ( rescale_mip )
      sum_of_x_store =  sum_of_x(iCut(1):iCut(1)+sizeChunk(1)-1,...
                                 iCut(2):iCut(2)+sizeChunk(2)-1,...
                                 iCut(3):iCut(3)+sizeChunk(3)-1);
      sum_of_x_store(vA(1,1) + 1:end - vA(2,1), ...
                     vA(1,2) + 1:end - vA(2,2), ...
                     vA(1,3) + 1:end - vA(2,3)) = gather(sum_of_x_tmp);   
     sum_of_x(iCut(1):iCut(1)+sizeChunk(1)-1,...
                  iCut(2):iCut(2)+sizeChunk(2)-1,...
                  iCut(3):iCut(3)+sizeChunk(3)-1) = sum_of_x_store;              
      clear sum_of_x_store
      
      sum_of_x2_store =  sum_of_x2(iCut(1):iCut(1)+sizeChunk(1)-1,...
                                 iCut(2):iCut(2)+sizeChunk(2)-1,...
                                 iCut(3):iCut(3)+sizeChunk(3)-1);
      sum_of_x2_store(vA(1,1) + 1:end - vA(2,1), ...
                     vA(1,2) + 1:end - vA(2,2), ...
                     vA(1,3) + 1:end - vA(2,3)) = gather(sum_of_x2_tmp);   
     sum_of_x2(iCut(1):iCut(1)+sizeChunk(1)-1,...
                  iCut(2):iCut(2)+sizeChunk(2)-1,...
                  iCut(3):iCut(3)+sizeChunk(3)-1) = sum_of_x2_store;              
      clear sum_of_x2_store      
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

if ( rescale_mip )
  sum_of_x = sum_of_x(1+tomoPre(1):end-tomoPost(1),...
                            1+tomoPre(2):end-tomoPost(2),...
                            1+tomoPre(3):end-tomoPost(3)) ./ currentGlobalAngle;
                          
  sum_of_x2 = sum_of_x2(1+tomoPre(1):end-tomoPost(1),...
                            1+tomoPre(2):end-tomoPost(2),...
                            1+tomoPre(3):end-tomoPost(3)) ./ currentGlobalAngle;

   %SAVE_IMG(sum_of_x,'sum_of_x.mrc');
   %SAVE_IMG(sum_of_x2,'sum_of_x2.mrc');
   %SAVE_IMG(RESULTS_peak,'prescaling.mrc');

  RESULTS_peak = RESULTS_peak - sum_of_x;
  sum_of_x  = sqrt(sum_of_x2 - sum_of_x.^2);
  clear sum_of_x2;
%   SAVE_IMG(sum_of_x,'stddev.mrc');
%   mov = mean(sum_of_x(:));
% %   sov = std(sum_of_x(:));
%   sov = 0;
%   sum_of_x(sum_of_x < (mov - 1*sov)) = max(sum_of_x(:));
%   SAVE_IMG(sum_of_x,'stddev_clipped.mrc');

  RESULTS_peak = RESULTS_peak ./ sum_of_x;
  
  RESULTS_peak = RESULTS_peak - mean(RESULTS_peak(:));
  RESULTS_peak = RESULTS_peak ./ rms(RESULTS_peak(:));
  clear sum_of_x;
end
gpuDevice(useGPU);
clear bhF


% scale the magnitude of the results to be 0 : 1
szK = latticeRadius;%floor(0.8.*szM);
rmDim = max(max(eraseMaskRadius),max(szK)).*[1,1,1]+7;
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
rmInt = interpolator(gpuArray(removalMask),[0,0,0],[0,0,0],rotConvention ,'forward','C1');
symOps = interpolator(gpuArray(removalMask),[0,0,0],[0,0,0],rotConvention ,'forward',symmetry);

maskCutOff = 0.98;
nIncluded = gather(sum(sum(sum(removalMask > maskCutOff))));
nTries = 0;
if strcmpi(eraseMaskType,'rectangle')
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


    % Switching from centered to lower left coordinates and subtracting the
    % padding 
    
    cenP = c + cMass' - rmDim;

    

% % %     % If the most frequent peak is unique use it;
% % %     [peakM, ~, peakC] = mode(angBox(:));
% % %     if length(peakC) == 1 && peakM
% % %       % Need to ensure the mode is none zero which is possible.
% % %       peakMat(n,4:6) = ANGLE_LIST(peakM,:);
% % %       topPeak = peakM;
% % %     else
      % Otherwise use the value at the max for the peak val;
      peakMat(n,4:6) = ANGLE_LIST(Ang(coord),:);
      topPeak = Ang(coord);
% % %     end
    peakMat(n,1:3) = gather(samplingRate.*cenP);
    peakMat(n,10) = gather(MAX);

    iSNR = 0;
    if nPeaks > 1
      possible_angles = gather(magBox);
      possible_angles(angBox == topPeak) = 0; 
      nRandom = 2;
      for iPeak = 2:nPeaks
        
        useRandom = false;
        
        if any(possible_angles ~= 0)
          [~, cAng] = max(possible_angles(:));
          topPeak = angBox(cAng);
          iSNR = gather(mean( possible_angles(angBox == topPeak)));
          possible_angles(angBox == topPeak) = 0; 
          Ang(cAng)
          if topPeak <= 0 || Ang(cAng) <= 0
            useRandom = true;
            % If random set the SNR as a fraction of the mean of the
            % previous peaks
            iSNR = mean( peakMat(n,10:10:10*(iPeak-1)) ) .* 0.75;
          else
            iAngles =  ANGLE_LIST(Ang(cAng),:);
          end
        else
          useRandom = true;
        end
        
        if useRandom
          % If we've used up all the possible peaks, just insert a random
          % Incrementally far from the original
          iAngles = [ randn(1) .* (nRandom.^2) + peakMat(n,1), ...
                      randn(1) .* (nRandom.^2) + peakMat(n,2), ...
                      randn(1) .* (nRandom.^2) + peakMat(n,3)];
          if nRandom < 10
            nRandom = nRandom + 1; 
          end
        end
        peakMat(n,[1:3]+10*(iPeak-1)) = gather(samplingRate.*cenP);
        peakMat(n,[4:6]+10*(iPeak-1)) = iAngles;     
        peakMat(n,10+10*(iPeak-1)) = iSNR;
% % %          oldPeaks = ( angBox == peakM | oldPeaks );
        
      end
      
    end


   
%     rmMask = BH_resample3d(removalMask,peakMat(n,4:6),[0,0,0],rotConvention ,'GPU','forward');
    rmMask = rmInt.interp3d(gather(peakMat(n,4:6)),[0,0,0],rotConvention,'forward','C1');

    % Invert after resampling so that zeros introduced by not extrapolating
    % the corners are swapped to ones, i.e. not removed.
%     rmMask = (1-rmMask);
    
    mag(c(1)-rmDim:c(1)+rmDim,...
        c(2)-rmDim:c(2)+rmDim,...
        c(3)-rmDim:c(3)+rmDim) = ...
    mag(c(1)-rmDim:c(1)+rmDim,...
        c(2)-rmDim:c(2)+rmDim,...
        c(3)-rmDim:c(3)+rmDim) .* (rmMask< maskCutOff);

% % %     peakMat(n,10) = (gather(MAX) - Tmean)./Tstd; % record stds above mean
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



n=1
for i = 1:length(peakMat(:,1))
   if all(peakMat(i,1:3))
        
        if SYMMETRY > 1
          % Generate a uniform distribution over the in-plane
          % randomizations
          iSym = rem( n + SYMMETRY, SYMMETRY)+1;
          r = reshape(BH_defineMatrix(peakMat(i,4:6), rotConvention , 'inv') *...
                      symOps.symmetry_matrices{iSym},1,9); 
        else
          r = reshape(BH_defineMatrix(peakMat(i,4:6), rotConvention , 'inv'),1,9);
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
              r = reshape(BH_defineMatrix(peakMat(i,[4:6]+10*(iPeak-1)), rotConvention , 'inv')*...
                          symOps.symmetry_matrices{iSym},1,9);
            else
              r = reshape(BH_defineMatrix(peakMat(i,[4:6]+10*(iPeak-1)), rotConvention , 'inv'),1,9);
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



fprintf('Total execution time : %f seconds\n', etime(clock, startTime));
   


end % end of templateSearch3d function

