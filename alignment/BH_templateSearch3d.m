function [hAvg, hRms, avgRange, rmsRange]  = BH_templateSearch3d( PARAMETER_FILE,...
                                        tomoName,tomoNumber,TEMPLATE, ...
                                        SYMMETRY, wedgeType, varargin)
 
                                                       
%3d template matching
%
%
%   Input variables:
%
%   MAP = Tomogram to be searched for motifs. The final binned size is limited
%         to roughly 5gb on a Tesla K40.
%
%   TEMPLATE = Reference volume containing the motif to search for.
%
%   RAWTLT = File containing the tilt angles as in IMOD convention (underfocus
%            magnitude increases with the X-axis for a positive tilt. These will
%            be used both to create a wedge mask to appropriately distort
%            the
%   
%   SAMPLING = Integer sampling, values of 3 or 4 are generally good for the
%              work I am doing, where the full sampling rate is <= 3A/pixel. For
%              gold standard work, the template matching step must be low-pass
%              filtered quite strongly to avoid biasing the data.
%
% 
%   PEAK_THRESHOLD = minimum multiple of standard deviations above the mean of
%                    cummulative covariance map that a peak needs to be in order
%                    to be counted as valid. By default, the maximum number of
%                    peaks to be accepted is 10,000- but this should only
%                    happen when no real peaks are found (by my thoughts
%                    anyhow!) Generally, it is advisable to test this with other
%                    parameters on your own data, in order to accept just
%                    slightly more (some false positives) than you expect, which
%                    can be later ignored via CCC or classification (or both.)
%
%   ANGLE_SEARCH = [a1 a2 a3 a4 helical] the angular search is a grid searched
%                  designed to exhaustively sample a unit sphere over the range
%                  you specify. Out of plane +/- a1 in steps of size a2, in
%                  plane +/- a3 in steps of a4. These together define a set of
%                  "latitudes" if you will, and the in plane sampling at each
%                  point sampled on that latitude, while the longitudinal
%                  sampling is calculated to be evenly sampled at the same rate
%                  as the latitude. The smaller the out of plane step size, the
%                  more longitudinal sampling points, and also the larger the
%                  absolute out of plane angle, the greater number of steps it
%                  takes to make a full revolution.
%     
%                  helical = 0 or 1 to indicate a helical rather than spherical
%                  grid search. I haven't actually tested this thoroughly so if
%                  set to 1 should be used with some extra caution.
%
%   TARGET_SIZE = [256 256 256] is the dimension, ideally a power of 2, 
%                  that we would like to pad each tomo chunk to, and this 
%                  depends on the available memory. In the future, I should set 
%                  this up as an automatic calculation. Minimal padding is size
%                  chunk + size reference. Should also be smaller than the 
%                  minimum dimenstion of you tomo.
%
%   LATTICE_RADIUS = [x,y,z] unbinned lattice dimension (i.e. particle radius)
%                    that will be excluded after a given peak is picked. In
%                    ANGSTROM
%
%
%   
%   Output variables:
%
%   None - saved to disk, .mat file with the geometry of best match found, as
%          well as the sampled version of the cummulative correlation map, which
%          can be viewed with the positions overlayed to asses the quality,
%          which can depend on not only the reference, but also the bandpass and
%          binning chosen as well as the angular search parameters.
%
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Goals & limitations:
%
%   3D template matching using a GPU. Developed on 4Gb card, consider this
%   a minimum memory requirement.
%
%   Tomogram and reference are read into memory.
%   Chunks of size ~ 384 x 384 x zDIM are processed at a time, which
%   are padded to next power 2 (512), leaving a minimal zero padding
%   buffer for the cross correlation. Ideally the padding would be twice
%   the image size, but this is just too large a memory requirement. 
%
%   The references are calculated on the fly, which takes a significant
%   portion of the computational effort, so larger chunks (fewer repeated
%   calculations) are preferable. See note in TODO section. 
%
%   Missing wedge region and low pass are taken from the tomogram chunk and
%   the rotated reference. They are then set to have mean 0 and std 1 in
%   the image domain, and compared by xcf. It is assumed the template
%   doesn't have a substantial amount of missing information itself.
%   
%   Each reference is given a reference that has an index between (0,pi)
%   and this is stored as the phase in the ccmap (which is inherently
%   real-valued.) Each iteration, the abs(ccmap) is compared to the new
%   ccmap, and updated to include and higher values, such that the final
%   peak search is performed on a "cummulative" ccmap which makes selecting
%   a threshold much easier. The default is mean + 6*std which has worked
%   well so far as judged by leaving only a few false positives (which are
%   easly removed in later classification in subvolume analysis.)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   TODO
%   1) Handle a list of multiple tomograms
%   2) Handle GPU mempory more dynamically
%   3) Handle variable missing wedge geometry
%   4) Determine balance of number of references/ number of chunks, vs.
%   memory requirements for storing the references in memory, and for large
%   grid searches, the benefit of calculating the references ahead of time,
%   even if the must be read into/out of memory onto scratch space. I think
%   for anything but an in plane search, this might be worthwhile.
%   Currently using MRC --> tif, reading in tif (individual images as
%   reading/writing multipage tifs with MATLAB takes forever.) 
%   5) Set-up to take low-pass cut off and pixel size as input
%   6) Test the benefit of filling in the wedge with low amplitude random
%   phase information as a way of decreasing the influence of the wedge.
%   7) Perform a number of test cases for different size tomo/template and
%   different step sizes in order to be able to more accurately estimate the
%   time needed for a given search range. This can also be used to boost the
%   usage on archer.3
%   8) Tomogram padding with numbers other than zero, I can't rember why i did
%   this initially, so ideally it should be figured out or changed.
%
%   - position - 0.5 for protomo display
%   - add field "total number" to track subtomonuber + a check to see if
%   geometry structure already exists
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAP = 'gag2_sirt10.mrc';
% PARAMETER_FILE = 'gagInitParam.m'
% TEMPLATE = 'briggsDup_scaleTrim.mrc';
% RAWTLT = './raw/GagSM4A_7.rawtlt';                                                                                               
% if (nargin ~= 5) 
%   if (nargin ~= 6 && nargin ~= 7)
%     error('args = PARAMETER_FILE,MAP,TEMPLATE, RAWTLT, SYMMETRY, [gpuIDX(1,2)], wedgeType')
%   end/fast_scratch/himesb/testRMS_45/cycle012_testRMS_15_class0_REF_ODD.mrc
% end

precision = 'single';
precisionTaper = 'singleTaper';

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
  mapBackIter = subTomoMeta.currentTomoCPR; clear subTomoMeta
catch
  mapBackIter = 0;
end
samplingRate  = pBH.('Tmp_samplingRate');

try
  tmpDecoy = pBH.('templateDecoy')
catch
  tmpDecoy = 0
end

if ( cmdLineThresh )
 peakThreshold = cmdLineThresh;
 fprintf('\nOverride peakThreshold from paramfile (%d) with cmd line arg (%d)\n\n',...
         cmdLineThresh, pBH.('Tmp_threshold'));
else
 peakThreshold = pBH.('Tmp_threshold');
end

latticeRadius = pBH.('particleRadius');
targetSize    = [512,512,512];%pBH.('Tmp_targetSize');
angleSearch   = pBH.('Tmp_angleSearch');

statsRadius = 1;

convTMPNAME = sprintf('convmap_wedgeType_%d_bin%d',wedgeType,samplingRate)

try 
  eraseMaskType = pBH.('Tmp_eraseMaskType');
catch
  eraseMaskType = 'sphere';
end
try
  eraseMaskRadius = pBH.('Tmp_eraseMaskRadius');
catch
  eraseMaskRadius = 0.75.*latticeRadius;
end

try
  load(sprintf('%s.mat', pBH.('subTomoMeta')), 'subTomoMeta')
  % Make certain this is a new alignment
  if length(fieldnames(subTomoMeta.cycle000)) > 1
    error('It looks like this is not a new project, remove mat file.')
  end
  
  % Make sure each subTomo has a unique idx
  nPreviousSubTomos = 0;
  fieldList = fieldnames(subTomoMeta.cycle000.geometry);
  for iField = 1:size(fieldList,1)
    nPreviousSubTomos = nPreviousSubTomos + ...
        size(subTomoMeta.cycle000.geometry.(sprintf('%s',fieldList{iField})),1);
  end
  subTomoMeta.('nSubTomoTotal') = nPreviousSubTomos;
catch 
  fprintf('\n\nDid not eveRot geometry.\n\n')
  subTomoMeta = struct();
  nPreviousSubTomos = 0;
end
reconScaling = 1;
try
  nPeaks = pBH.('nPeaks');
catch
  nPeaks = 1;
end

pixelSizeFULL = pBH.('PIXEL_SIZE').*10^10;
if pBH.('SuperResolution')
  pixelSizeFULL = pixelSizeFULL * 2;
end

pixelSize = pixelSizeFULL.*samplingRate;

% For testing
try 
  lowResCut = pBH.('lowResCut');
catch
  % Current default, which may be too conservative.
  lowResCut = 40;
end

% % %  % Read and/or bin the tomogram and template
% % % [tomogram,  mapPath, mapName, mapExt] = ...
% % %                                    BH_multi_eveRotOrBin( MAP, samplingRate, 3 );

mapPath = './cache';
mapName = sprintf('%s_%d_bin%d',tomoName,tomoNumber,samplingRate);
mapExt = '.rec';

sprintf('recon/%s_recon.coords',tomoName)
[ recGeom, ~, ~] = BH_multi_recGeom( sprintf('recon/%s_recon.coords',tomoName) );

reconCoords = recGeom(tomoNumber,:);
clear recGeom

    
[ tomogram ] = BH_multi_loadOrBuild( sprintf('%s_%d',tomoName,tomoNumber),  ...
                                     reconCoords, mapBackIter, samplingRate,...
                                     -1*gpuIDX, reconScaling,1); 
                                           
 
% We'll handle image statistics locally, but first place the global environment
% into a predictible range

% Wait until after application of/distortion by tomo's ctf to downsample the template. 
% % % [template, tempPath, tempName, tempExt] = ...
% % %                               BH_multi_loadOrBin( TEMPLATE, samplingRate, 3 );   
[template, tempPath, tempName, tempExt] = ...
                              BH_multi_loadOrBin( TEMPLATE, 1, 3 ); 
                            
% The template will be padded later, trim for now to minimum so excess
% iterations can be avoided.
fprintf('size of provided template %d %d %d\n',size(template));
trimTemp = BH_multi_padVal(size(template),ceil(1.15.*latticeRadius./pixelSizeFULL));
template = BH_padZeros3d(template, trimTemp(1,:),trimTemp(2,:),'cpu','singleTaper');
clear trimTemp
fprintf('size after trim to 1.15x lattice radius %d %d %d\n',size(template));
                            
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


sizeTemp = size(template)
sizeTempBIN = size(templateBIN)

% % borderMask =  gather(BH_mask3d('rectangle',sizeTemp,(sizeTemp./2)-6,[1,1,1])); 

% % % interpMask = BH_mask3d_cpu('sphere',sizeTemp,max(latticeRadius./pixelSizeFULL).*[1,1,1],[0,0,0]);
% % % interpMask = interpMask > 0.01;

% Calculate exclusion area around chosen peaks.
% Convert Ang to pixels
maxSizeForHighPass = 2*max(latticeRadius);

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

[ padVal ] = BH_multi_padVal( sizeTemp, sizeChunk );
tempPre = padVal(1,:);
tempPost = padVal(2,:);

[ padBIN ] = BH_multi_padVal( sizeTempBIN, sizeChunk );
[ trimValid ] = BH_multi_padVal(sizeChunk, validArea);

if ( tmpDecoy )
  % This is probably sample dependent. should search a small range and find
  % the maximum rate of change in the ccc
  
  % the -1 searches for the next smallest fast fourier size
  decoySize = BH_multi_iterator(-1.*floor(tmpDecoy.*sizeChunk),'fourier');
  fprintf('\n\ndecoy size %d %d %d from chunkSize %d %d %d\n\n',decoySize,sizeChunk);  
  padDecoy1 = BH_multi_padVal( sizeTempBIN, decoySize);
  padDecoy2 = BH_multi_padVal( decoySize, sizeChunk );
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
                                                     'cpu', 'singleTaper');
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


% optimize fft incase a power of two is not used, this will make things run ok.
opt = zeros(sizeChunk, precision,'gpuArray');
fftw('planner','patient');
fftn(opt);
clear opt ans



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Preprocess the tomogram

tomoIDX = 1;
nTomograms = prod(nIters);


tomoStack = zeros([sizeChunk,nTomograms], 'single');

% backgroundVol = zeros(sizeChunk,'single');
tomoCoords= zeros(nTomograms, 3, 'uint16');

[ tomoBandpass ]   = BH_bandpass3d(sizeChunk, 0,maxSizeForHighPass, ...
                                              lowResCut,'cpu', pixelSize );
         
try
  doMedFilt = pBH.('Tmp_medianFilter')
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

for statLoop = 1:1+calcStats
for  iX = 1:nIters(1)
  cutX = 1 + (iX-1).*validArea(1);
  for iY = 1:nIters(2)
    cutY = 1 + (iY-1).*validArea(2);
    for iZ = 1:nIters(3)
      cutZ = 1 + (iZ-1).*validArea(3);

    fprintf('preprocessing tomo_chunk %d/%d col %d/%d row %d/%d plane idx%d\n' , ...
                                  iY,nIters(2),iX,nIters(1),iZ,nIters(3),tomoIDX)


    % Cut out chunk and zero pad - double would be more accurate, but for
    % template matching which is fairly crude anyhow, this should be okay, and
    % allows much larger chunks to be processed.
    
    
    tomoChunk = tomogram(cutX:cutX+sizeChunk(1)-1,...
                         cutY:cutY+sizeChunk(2)-1,...
                         cutZ:cutZ+sizeChunk(3)-1);
                                       
              
                                  
    if ( statLoop < 2 )             
      % Calc stats on something smoothed      
      nonZero = tomoChunk(:) ~= 0;
    end
    
    
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

    % Handle all mean centering and rms normalization in local window
    
    [ averageMask, flgOOM ] = BH_movingAverage(tomoChunk, statsRadius); 
    
    if isa(tomoChunk(1),'gpuArray') && flgOOM
      tomoChunk = gather(tomoChunk);
    end
    
    if ( calcStats )
      if ( firstStats )      
        minAvg = min(averageMask(nonZero));
        minAvg = minAvg - 1.25*abs(minAvg);
        maxAvg = max(averageMask(nonZero));
        maxAvg = maxAvg + 1.25*abs(maxAvg);
        avgRange = minAvg:(maxAvg-minAvg)/5000:maxAvg;
        hAvg = hist(averageMask(nonZero),avgRange);
      else
        hAvg = hAvg + hist(averageMask(nonZero),avgRange);
      end
    else
%       maskStack(:,:,:,tomoIDX) = ( averageMask < lowAvg | averageMask > highAvg );
    end
    
    
    tomoChunk= tomoChunk - averageMask; clear averageMask 

    [ rmsMask ] = BH_movingRMS(tomoChunk, statsRadius);
    
    if ( statLoop > 1 )
      maskStack(:,:,:,tomoIDX) = ( maskStack(:,:,:,tomoIDX) | rmsMask > highRms | rmsMask < lowRms);
    end      
      % Using the non-ctf corrected stack since we limit toA all practical
      % defocus (<8um) should be entirely negative, so just flip in real
      % space
   
      tomoChunk = -1.*(tomoChunk ./ rmsMask);
%       backgroundVol = backgroundVol + gather(tomoChunk.*maskStack(:,:,:,tomoIDX));

    
    if ( calcStats )
      if ( firstStats )
        minRms = min(rmsMask(nonZero));
        minRms = minRms - 1.25*abs(minRms);
        maxRms = max(rmsMask(nonZero));
        maxRms = maxRms + 1.25*abs(maxRms);   
        rmsRange = minRms:(maxRms-minRms)/5000:maxRms;
        hRms = hist(rmsMask(nonZero),rmsRange);
      else
        hRms = hRms + hist(rmsMask(nonZero),rmsRange);
      end
      
    end
    clear rmsMask
    firstStats = 0;
    
    if ~( calcStats )
     % Even if single precision is set, pull to CPU and bandpass at double
     % precision, reverting to single.
      if strcmpi(precision, 'double')
        tomoChunk = double(gather(tomoChunk));
        tomoStack(:,:,:,tomoIDX) = single(real(ifftn(fftn(...
                                          validAreaMask.*tomoChunk).* ...
                                          tomoBandpass)));
      else

       tomoStack(:,:,:,tomoIDX) = gather(single(real(ifftn(fftn(...
                                          validAreaMask.*tomoChunk).* ...
                                          tomoBandpass))));
      end
    end
     
    if ~( calcStats )
      tomoCoords(tomoIDX,:) = [cutX,cutY,cutZ];
      tomoIDX = tomoIDX + 1;
    end
  
    end
  end
end

if ( calcStats )
  calcMask = 1;
  calcStats = 0;
  hRms = gather(hRms); hAvg = gather(hAvg);
  rmsRange = gather(rmsRange); avgRange = gather(avgRange);

  % Remove begining and end of range
  hAvg(1) = 0; hAvg(end) = 0;
  hRms(1) = 0; hRms(end) = 0;


  [maxRms,maxRmsCoord] = max(hRms);
  [maxAvg,maxAvgCoord] = max(hAvg);

  rmsWin = find(hRms>0.6*maxRms,1,'last') - find(hRms>0.6*maxRms,1,'first');
  avgWin = find(hAvg>0.3*maxAvg,1,'last') - find(hAvg>0.3*maxAvg,1,'first');

  lowRms = rmsRange(maxRmsCoord-rmsWin);
  highRms= rmsRange(maxRmsCoord+rmsWin);

  lowAvg = avgRange(maxAvgCoord-avgWin);
  highAvg= avgRange(maxAvgCoord+avgWin);

  % figure, plot(rmsRange,hRms);
  hRms(1:maxRmsCoord-rmsWin) = 0;hRms(maxRmsCoord+rmsWin:end) = 0;


  hAvg(1:maxAvgCoord-avgWin) = 0;hAvg(maxAvgCoord+avgWin:end) = 0;

  
  maskStack = gather(maskStack);
end


end


  
clear tomoWedgeMask averagingMask rmsMask bandpassFilter statBinary validAreaMask tomoChunk

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Temp while testing new dose weighting
      TLT = tiltGeometry;
      nPrjs = size(TLT,1);

      
      kVal = 0;
      [ OUTPUT ] = BH_multi_iterator( [sizeChunk;kVal.*[1,1,1]], 'extrapolate' )


     

      % The negative size is to overried the default cubic dimension enforcement  
      %[ wedgeMask ]= BH_weightMask3d(-1.*OUTPUT(1,:), tiltGeometry, 'applyMask',particleThickness,1, 1, samplingRate);



      switch wedgeType
        case 1
          % make a binary wedge
          [ wedgeMask ]= BH_weightMask3d(-1.*OUTPUT(1,:), tiltGeometry, ...
                         'binaryWedgeGPU',particleThickness,...
                                                       1, 1, samplingRate);
        case 2
          % make a non-CTF wedge
          [ wedgeMask ]= BH_weightMask3d(-1.*OUTPUT(1,:), tiltGeometry, ...
                         'applyMask',particleThickness,...
                                                       2, 1, samplingRate);   
        case 3
          % make a CTF without exposure weight
          [ wedgeMask ]= BH_weightMask3d(-1.*OUTPUT(1,:), tiltGeometry, ...
                         'applyMask',particleThickness,...
                                                       3, 1, samplingRate);   
        case 4
          % make a wedge with full-ctf
          [ wedgeMask ]= BH_weightMask3d(-1.*OUTPUT(1,:), tiltGeometry, ...
                         'applyMask',particleThickness,...
                                                       4, 1, samplingRate);  
        otherwise
          error('wedgeType must be 1-4');
      end
      
      wedgeMask = gather(ifftshift(wedgeMask));
                                     
      
        
 
        
    
 
  
  
    

currentGlobalAngle = 1;
ANGLE_LIST = zeros(nAngles(1),3, 'single');
nComplete = 0;
totalTime = 0;
firstLoopOverTomo = true;
for iAngle = 1:size(angleStep,1)
  
  theta = angleStep(iAngle,1);

  % Calculate the increment in phi so that the azimuthal sampling is
  % consistent and equal to the out of plane increment.
  if (helical)
     phi_step = 360;
  else
     phiStep = angleStep(iAngle,3);
  end
  
  gpuDevice(useGPU);

  %numRefIter = nAngles(1);
  numRefIter = angleStep(iAngle,2)*length(inPlaneSearch)+1;
  tempImg = gpuArray(template);

% % %   interpMaskGPU = gpuArray(interpMask);
  

% % %   tempWedgeMask = tempWedgeMask .^ (0.25./tempWedgeMask);    

  % SAVE_IMG(MRCImage(gather(fftshift(tempWedgeMask))), 'tmpWedge.mrc') ;

  [ tempFilter ] = BH_bandpass3d(sizeChunk,0.1, maxSizeForHighPass, ...
                                          lowResCut,'GPU', pixelSizeFULL );

  tempWedgeMask = gpuArray(wedgeMask) .* tempFilter; 

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


     % fftn(double(gpuArray))) ~ 2.5x faster than transfering a double
     % complex
     if strcmpi(precision, 'double')
      tomoFou = fftn(double(gpuArray(tomoStack(:,:,:,tomoIDX))));
     else
      tomoFou = fftn(gpuArray(tomoStack(:,:,:,tomoIDX)));
     end
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



          tempFou = BH_bandLimitCenterNormalize(tempRot,tempWedgeMask,'',[tempPre;tempPost],precisionTaper);

          tempRot = BH_padZeros3d(real(ifftn(tempFou)),-1.*tempPre,-1.*tempPost,'GPU',precision); 
          
          tempRot = gather(BH_reScale3d(tempRot,'',sprintf('%f',1/samplingRate),'GPU'));

          if (firstLoopOverTomo)
            SAVE_IMG(MRCImage(tempRot), sprintf('temp_%s.mrc',convTMPNAME),pixelSize);
          end
          referenceStack(:,:,:,intraLoopAngle) = tempRot;
         
          tempFou = fftn(BH_padZeros3d(tempRot,padBIN(1,:),padBIN(2,:),'GPU',precision));
 
          
          
          
        else

          tempFou = (fftn(BH_padZeros3d( ...
                              referenceStack(:,:,:,intraLoopAngle), ...
                              padBIN(1,:), padBIN(2,:),'GPU', precision)));

        end


% % %         ccfmapFull = fftshift(real(single(ifftn(tomoFou.*conj(tempFou)))));
% % %                                    
% % % 
% % %         ccfmap = ccfmapFull(vA(1,1) + 1:end - vA(2,1), ...
% % %                             vA(1,2) + 1:end - vA(2,2), ...
% % %                             vA(1,3) + 1:end - vA(2,3));
        
        % Even with local normalization, test with all padding and
        % goodness.
%         tomoNorm = ((sqrt(sum(sum(sum(abs(tomoFou).^2)))) ./ numel(tomoFou)));
        tempNorm = ((sqrt(sum(sum(sum(abs(tempFou).^2)))) ./ numel(tempFou)));
        
        ccfmap = BH_padZeros3d(fftshift(real(single(...
                               ifftn(tomoFou.*conj(tempFou))./tempNorm))),...%./(tomoNorm.*tempNorm))))),...
                               trimValid(1,:),trimValid(2,:),'GPU',precision);
        
        if ( tmpDecoy > 0 )
           tempFou = [];
           if (firstLoopOverAngle)
             decoy = BH_padZeros3d(conj(fftn(BH_padZeros3d(tempRot, ...
                                padDecoy1(1,:),padDecoy1(2,:), ...
                                'GPU',precision))), ...
                                padDecoy2(1,:),padDecoy2(2,:),...
                                'GPU',precision,0,1);
           else
             decoy = BH_padZeros3d(conj(fftn(BH_padZeros3d( ...
                                referenceStack(:,:,:,intraLoopAngle), ...
                                padDecoy1(1,:),padDecoy1(2,:), ...
                                'GPU',precision))), ...
                                padDecoy2(1,:),padDecoy2(2,:),...
                                'GPU',precision,0,1);             
           end
           
           

           
           decoy = BH_padZeros3d(fftshift(real(single( ...
                                 ifftn(tomoFou.*decoy)))),..../(decoyNorm.*tomoNorm))))),
                                 trimValid(1,:), ...
                                 trimValid(2,:),'GPU',precision);
        elseif ( tmpDecoy < 0 )
          
          % Just use the mirror image of the template, i.e. take the conj
          % (of the conj) so just the padded FFT of the ref.
           decoy = BH_padZeros3d(fftshift(real(single( ...
                                 ifftn(tomoFou.*tempFou)))),..../(decoyNorm.*tomoNorm))))),
                                 trimValid(1,:), ...
                                 trimValid(2,:),'GPU',precision);          
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
  RESULTS_decoy = RESULTS_decoy ./ std(RESULTS_decoy(:));
  RESULTS_decoy(RESULTS_decoy < 0) = 0; 
  
end
gpuDevice(useGPU);


% scale the magnitude of the results to be 0 : 1
szK = latticeRadius;%floor(0.8.*szM);
rmDim = max(szK).*[1,1,1];
mag = RESULTS_peak; clear RESULTS_peak
% Normalize so the difference if using a decoy makes sense. The input decoy
% should have the same power, so I'm not sure why this is needed, but it is
% an easy fix and a problem for future Ben to figure out.
mag = mag ./ std(mag(:));

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
  mag = mag - RESULTS_decoy; clear RESULTS_decoy
  SAVE_IMG(MRCImage((mag)),diffOUT);
end
angleFILE = fopen(angleListOUT,'w');
fprintf(angleFILE,'%2.2f\t%2.2f\t%2.2f\n', ANGLE_LIST');
fclose(angleFILE);


mag =  mag - min(mag(:)); mag = mag ./ max(mag(:));

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
%%%removalMask = (removalMask < 1);
while  n <= peakThreshold

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
c = gather([i,j,k]);

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
        c(3)-rmDim:c(3)+rmDim) .* ~(rmMask>0.95);

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

