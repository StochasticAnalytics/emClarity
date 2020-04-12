function [  ] = BH_average3d(PARAMETER_FILE, CYCLE, STAGEofALIGNMENT)
%Extract and interpolate a subTomogram from a 3d volume.
%
%   Input variables:
%
%   CYCLE = 0,1,2 etc.
%
%   NEW_CYCLE = 0,1 truth value, use the specified geometry to initialize a new
%               alignment cycle.
%
%   IMAGE = 3d volume, or a string specifing a volume to read in.exit

%
%   classVector = Extract a subset of class averages.
%              0 = Ignore classes, and average all
%
%              2, # classes = Extract specified classes with unique symmetry.
%                                                      [1, 2, 5, 6;
%                                                       3, 1, 6, 1]
%
%   CLASS_NAME = class to draw from
%
%   MOTIF_SIZE = Size of the window to extract. This must be at least 14 pixels
%                larger than the key features to allow for a soft apodization
%                over 7 pixels in each direction.
%
%
%   SAMPLING = Binning factor, assumed to be integer value. Image is first
%              smoothed by an appropriate low-pass filter to reduce aliasing.
%
%   GEOMETRY = A structure with tomogram names as the field names, and geometry
%              information in a 26 columSn array.
%              Additionally, a field called 'source_path' has a value with the
%              absolute path to the location of the tomograms.
%
%              The input is a string 'Geometry_templatematching.mat' for
%              example, and it is expected that the structure is saved as the
%              variable named geometry.
%
%
%   FSC = Randomly divide the data into halves (not just even/odd) for use in
%         Fourier shell correlation calculation.
%
%   OUTPUT_PREFIX = String to prepend to output volumes.
%
%
%   Output variables:
%
%   None = files are written to disk in the current directory.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Goals & Limitations:
%
%   Cut out a subTomogram and transform to the standard basis (microscope
%   reference frame) from an existing 3d volume. 
%
%   Assumed to run on GPU.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   TODO
%     - Error checking for memory limitations     
%     - In testing verify "implicit" gpu arrays are actually gpu arrays
%     - Check binning
%     - Confirm position 7 is where I want to keep FSC value
%     - Update geometry to record % sampling
%
%     - Store binning & Size in class specific. geometry and use as a check for other programs
%     at runtime.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (nargin ~= 3)
  error('args = PARAMETER_FILE, CYCLE, STAGEofALIGNMENT')
end

try
  CUTPADDING = subTomoMeta.('CUTPADDING')
catch
  CUTPADDING=20
end

% Explicit reference to location of variables in main memory, or on the GPU.
cpu = struct();
GPU = struct();

startTime =  clock;
CYCLE = str2num(CYCLE);

% Until people forget this existed
if strcmpi(STAGEofALIGNMENT,'NoAlignment')
  STAGEofALIGNMENT = 'RawAlignment';
end

cycleNumber = sprintf('cycle%0.3u', CYCLE)

pBH = BH_parseParameterFile(PARAMETER_FILE);
load(sprintf('%s.mat', pBH.('subTomoMeta')), 'subTomoMeta');
reconScaling = 1;
try
  nPeaks = pBH.('nPeaks');
catch
  nPeaks = 1;
end


mapBackIter = subTomoMeta.currentTomoCPR;

if (CYCLE)
  try
    flgQualityWeight = pBH.('flgQualityWeight');
  catch
    flgQualityWeight = 4;
  end
else
  fprintf('No quality weighting in the initial cycle after template matching\n');
  flgQualityWeight = 0;
end

% Experimental downweighting of higher frequency info farther from focus.
% Could also consider filtering pre reconstruction
try 
  flgFilterDefocus = pBH.('filterDefocus');
  fprintf('\nFiltering by defocus using exp[-(%d*(argmax(def-1,0,5).*q)^%d)]\n',flgFilterDefocus);
catch
  flgFilterDefocus = 0;
end

try
  flgCutOutVolumes = pBH.('flgCutOutVolumes');
catch
  flgCutOutVolumes = 0;
end

doCut = 0
if (flgCutOutVolumes)
  if isfield(subTomoMeta,'volumesAreCutOut')
    if ~(subTomoMeta.volumesAreCutOut)
      doCut = 1
    end
  else
    doCut = 1
  end
end


flgClassify= pBH.('flgClassify');
%%% For general release, I've disabled class average alignment and
%%% multi-reference alignment, so set the default to OFF. If either of
%%% these features are re-introduced, this will need to be reverted.
if ( flgClassify ); flgClassify = -1 ; end
try
  flgMultiRefAlignment = pBH.('flgMultiRefAlignment');
catch
  flgMultiRefAlignment = 0;
end
flgGold=1;
pixelSize = pBH.('PIXEL_SIZE').*10^10;
if pBH.('SuperResolution')
  pixelSize = pixelSize * 2;
end
nGPUs = pBH.('nGPUs');
% Optionally specify gpu idxs
if numel(nGPUs) == 1
  gpuList = 1:nGPUs;
else
  gpuList = nGPUs;
  nGPUs = length(gpuList);
end



try
  loadTomo = pBH.('loadTomo')
catch
  loadTomo = 0
end

try
  scaleCalcSize = pBH.('scaleCalcSize');
catch
  scaleCalcSize = 1.5;
end

try
  use_v2_SF3D = pBH.('use_v2_SF3D')
catch
  use_v2_SF3D = false
end
global bh_global_ML_compressByFactor;
global bh_global_ML_angleTolerance;
if isempty(bh_global_ML_compressByFactor)
  bh_global_ML_compressByFactor = 2.0;
end
if isempty(bh_global_ML_angleTolerance)
  bh_global_ML_angleTolerance = 5;
end
  
if (nPeaks > 1)
  fprintf('For ML approach:\nUsing a compression factor %3.3f\nUsing an angulare tolerance of %3.3f degrees\n', ...
          bh_global_ML_compressByFactor, bh_global_ML_angleTolerance);
end
% for now only turn on (optionally) in reference generation.

classVector = cell(2,1);
classCoeffs = cell(2,1);
refVector = cell(2,1);
flgFinalAvg = 0;
oddRot = 0;
refIMG = cell(1,1);
refWDG = cell(1,1);
peakMask = [];
peakBinary = [];
peakCOM = [];
peakSearch = [];
saveClassSum = -1;
switch STAGEofALIGNMENT
  case 'RawAlignment'
    
    
    if (flgClassify)
      fieldPrefix = 'Raw'

    else
      fieldPrefix = 'REF'
    end


      classVector{1}  = pBH.(sprintf('%s_classes_odd','Raw'));
      classVector{2}  = pBH.(sprintf('%s_classes_eve','Raw'));

      className    = pBH.(sprintf('%s_className','Raw'));
      samplingRate = pBH.('Ali_samplingRate'); 
      
     if (flgMultiRefAlignment && ~flgClassify)
       className    = pBH.(sprintf('Raw_className'))
       saveClassSum = pBH.(sprintf('Raw_className'))
     elseif (flgMultiRefAlignment && flgClassify)
       fprintf('\n\nMutliRef and Classify enabled.\n');
       fprintf('Only creating the global class average for PCA\n\n.');
       className = 0;
       saveClassSum = 0;
       classVector{1} = [0;1];
       classVector{2} = [0;1];
     end
       
    
    
  case 'ClassAlignment'
    fieldPrefix = 'Cls';

    classVector{1}  = pBH.(sprintf('%s_classes_odd',fieldPrefix));
    classVector{2}  = pBH.(sprintf('%s_classes_eve',fieldPrefix));

    samplingRate = pBH.(sprintf('Cls_samplingRate'));
    classVector{1}  = pBH.(sprintf('%s_classes_odd',fieldPrefix));
    className    = pBH.(sprintf('%s_className',fieldPrefix));

  case 'NoAlignment'
    % Special case for the first cycle
    fieldPrefix = 'Raw';
    
    classVector{1}  = pBH.(sprintf('%s_classes_odd',fieldPrefix));
    classVector{2}  = pBH.(sprintf('%s_classes_eve',fieldPrefix));

    className    = pBH.(sprintf('%s_className',fieldPrefix));
    samplingRate = pBH.('Ali_samplingRate');
    if (flgClassify)
     % samplingRate = pBH.('Pca_samplingRate');
    else
      %samplingRate = pBH.('Raw_samplingRate');
      fieldPrefix = 'REF';
    end
    
  case 'FinalAlignment'
    % Special case for the final cycle.
    % Assuming RawAlignment already run for this cycle and FSC is calculated
    % Goal is to re-extract odd-half, applying the xform found in fscGold
    fieldPrefix = 'Raw'
   
    classVector{1}  = pBH.(sprintf('%s_classes_odd',fieldPrefix));
    classVector{2}  = pBH.(sprintf('%s_classes_eve',fieldPrefix));
    
    className    = pBH.(sprintf('%s_className',fieldPrefix));
    samplingRate = pBH.('Ali_samplingRate');
    if (flgClassify)
      %samplingRate = pBH.('Pca_samplingRate');
    else
      %samplingRate = pBH.('Raw_samplingRate');
      fieldPrefix = 'REF'
    end
    
    flgFinalAvg = 1;
    % Update at some point to handle multiple classes, for now just test on the
    % "global" or whatever requested class
    iRefPrev = 1;
    

    try
      
      aliParams = subTomoMeta.(cycleNumber).('fitFSC').(sprintf('Resample%s%d','REF',iRefPrev));
      oddRot = reshape(aliParams(1,:),3,3);
      % refine the translation per particle.
    catch
      fprintf('\nReverting from %s to Raw in loading fitFSC\n','REF');
     
      aliParams = subTomoMeta.(cycleNumber).('fitFSC').(sprintf('Resample%s%d','Raw',iRefPrev));
      
      oddRot = reshape(aliParams(1,:),3,3);
      % refine the translation per particle.    
    end
    clear iRefPrev
  


  case 'Cluster_cls'
    STAGEofALIGNMENT = 'Cluster';
    ClusterGeomNAME = 'ClusterClsGeom';
    fieldPrefix = 'Cls';

    classVector{1}  = pBH.(sprintf('%s_classes_odd',fieldPrefix));
    classVector{2}  = pBH.(sprintf('%s_classes_eve',fieldPrefix));

    classCoeffs{1} =  pBH.('Pca_coeffs');
    classCoeffs{2} =  pBH.('Pca_coeffs');

    samplingRate = pBH.(sprintf('Cls_samplingRate'));
    className    = pBH.(sprintf('%s_className',fieldPrefix));
    if flgClassify < 0
      flgGold = 0;
    end     
  case 'Cluster_ref'
    STAGEofALIGNMENT = 'Cluster';
    ClusterGeomNAME = 'ClusterRefGeom';
    fieldPrefix = 'Ref';
  
    classVector{1}  = pBH.(sprintf('%s_classes_odd',fieldPrefix));
    classVector{2}  = pBH.(sprintf('%s_classes_eve',fieldPrefix));

    classCoeffs{1} =  pBH.('Pca_coeffs');
    classCoeffs{2} =  pBH.('Pca_coeffs');

    samplingRate = pBH.('Cls_samplingRate');
    className    = pBH.(sprintf('%s_className',fieldPrefix));
    if flgClassify < 0
      flgGold = 0;
    end     
  case 'RefAlignment'
    
    fieldPrefix = 'REF';
      
    refVector{1}  = pBH.(sprintf('%s_references_odd','Ref'));
    refVector{2}  = pBH.(sprintf('%s_references_eve','Ref'));

    samplingRate = pBH.(sprintf('%s_samplingRate','Cls'));
    % note the exchange of class to ref

    for iGold = 1:2
      % get just unique group names, which are subbed in for class id
      [iGroup, groupIDX, ~] = unique(refVector{iGold}(3,:));
      classVector{iGold} = [ iGroup ; refVector{iGold}(2, groupIDX) ];
    end
    classVector{1}
    classVector{2}
    className    = pBH.(sprintf('%s_className','Ref'));
    
  case 'SnrEstimate'
    
    if CYCLE
      fieldPrefix = 'Raw';
    else
      fieldPrefix = 'NoA';
    end
    classVector{1}  = [1:25;ones(1,25)];
    classVector{2}  = [1:25;ones(1,25)];

    className    = 25;
    samplingRate = pBH.(sprintf('%s_samplingRate','Ali'));
    
  otherwise
    error('STAGEofALIGNMENT incorrect')
end

fprintf('StOAlign = %s, fieldPrefix = %s\n', STAGEofALIGNMENT, fieldPrefix);


flgCones     = pBH.('flgCones');

try
  cccCutOff    = pBH.('flgCCCcutoff');
catch
  cccCutOff = 0.0;
end
cutPrecision = 'single'; %pBH.('flgPrecision');
try
  interpOrder = pBH.('interpOrder');
catch
  interpOrder = 1;
end

try 
  flgLimitToOneProcess = pBH.('flgLimitToOneProcess');
catch
  flgLimitToOneProcess = 0;
end

if ( loadTomo )
  limitToOne = loadTomo;
  if (flgLimitToOneProcess)
    limitToOne = min(limitToOne, flgLimitToOneProcess);
  end
elseif interpOrder == 4
  limitToOne = 1;
elseif (flgLimitToOneProcess)
  limitToOne = flgLimitToOneProcess;
else
  limitToOne = pBH.('nCpuCores'); 
  interpOrder = 1;
end





fprintf('Interporder %d, limitToOneProcess %d\n',interpOrder,limitToOne);

if ~(ismember(interpOrder,[1,4]))
  error('interpolationOrder must be 1,,4 - linear,sinc');
end
outputPrefix = sprintf('%s_%s',cycleNumber, pBH.('subTomoMeta'));

pixelSize = pixelSize .* samplingRate;
peakSearch   = floor(0.85.*pBH.('particleRadius')./pixelSize)
peakCOM      = [1,1,1].*3;
    
if strcmp(cutPrecision, 'single')  
  singlePrecision = true
elseif strcmp(cutPrecision, 'double')
  singlePrecision =false
else
  error('cutPrecision must be single or double')
end
  

  
if ~(isfield(subTomoMeta,'currentCycle'))
  subTomoMeta.('currentCycle') = 0;
end
  
% Somewhere I am saving currentCycle as a string. Haven't taken the time to
% track it down, but probably in this function.
if isa(subTomoMeta.('currentCycle'), 'char')
  subTomoMeta.('currentCycle') = str2num(subTomoMeta.('currentCycle'))
end

  if subTomoMeta.('currentCycle') == CYCLE - 1 
    
    cycleRead = sprintf('cycle%0.3u', CYCLE - 1);
    % Save a backup of the cycles total geometry
    save(sprintf('%s_%s_backup.mat',cycleRead,pBH.('subTomoMeta')), ...
                                                                 'subTomoMeta');
  end
  
  if strcmpi(STAGEofALIGNMENT, 'RawAlignment')
    if ( CYCLE )
      cycleRead = sprintf('cycle%0.3u', CYCLE - 1)
    else
      cycleRead = sprintf('cycle%0.3u', CYCLE)
    end
  else
    cycleRead = sprintf('cycle%0.3u', CYCLE )
  end
  
  % leave averages at size appropriate for interpolation when extracting to use
  % in class avg alignment.
  doNotTrim = false;
  eachTomo  = false;
  flgEstSNR = 0;
  switch STAGEofALIGNMENT
    case 'ClassAlignment'
      geometry = subTomoMeta.(cycleRead).ClassAlignment;
%%%    case 'NoAlignment'
%%%      geometry = subTomoMeta.(cycleRead).geometry; 
      % Extract an average from each tomogram for visual inspection, filter to
      % 0.75xtemplateSearch_lowpass
%%%      eachTomo = true;
% %%      lpTomo = pBH.('Tmp_bandpassFilter');
%%%      doNotTrim = true;
        
    case 'RawAlignment' 
      if ( CYCLE )
        geometry = subTomoMeta.(cycleRead).RawAlign;
      else
        geometry = subTomoMeta.(cycleRead).geometry;
        eachTomo = false;%true;
      end
      if ~(flgClassify)
        doNotTrim = true;
      end
      
    case 'FinalAlignment'
      geometry = subTomoMeta.(cycleRead).Avg_geometry;
      if ~(flgClassify)
        doNotTrim = true;
      end      

    case 'Cluster'
      cN = cell(2,1);

      if (flgGold)
        cN{1} = sprintf('%s_%d_%d_nClass_%d_ODD',outputPrefix,classCoeffs{1}(1,1), ...
                                          classCoeffs{1}(1,end), className)

        geometry{1} = subTomoMeta.(cycleRead).ClusterResults.(cN{1});                                
        cN{2} = sprintf('%s_%d_%d_nClass_%d_EVE',outputPrefix,classCoeffs{2}(1,1), ...
                                          classCoeffs{2}(1,end), className)
        geometry{2} = subTomoMeta.(cycleRead).ClusterResults.(cN{2});  

        geometry = BH_mergeClassGeometry(geometry{1}, geometry{2});
      else
         cN{1} = sprintf('%s_%d_%d_nClass_%d_STD',outputPrefix,classCoeffs{1}(1,1), ...
                                          classCoeffs{1}(1,end), className)

        geometry = subTomoMeta.(cycleRead).ClusterResults.(cN{1});        
      end
      
      
      subTomoMeta.(cycleRead).('KmsSampling') = samplingRate;
      doNotTrim = true;
    case 'RefAlignment'
 
      
      geometry = subTomoMeta.(cycleRead).('RefAlignment')
      subTomoMeta.(cycleRead).('Cls_Sampling') = samplingRate;
      doNotTrim = true;
    case 'SnrEstimate'
      flgEstSNR = 1;
      if ( CYCLE )
        geometry = subTomoMeta.(cycleRead).RawAlign;
      else
        geometry = subTomoMeta.(cycleRead).geometry; 
      end
      % randomly dived all currently included subTomos into 10 bins
     [geometry, nTotal, snrBinSize] = BH_randomSubset(geometry,'snr',-1,[1]);
     [geometry, nTotal, snrBinSize] = BH_randomSubset(geometry,'snr',-1,[2]);
   
    otherwise
      error('STAGE_ALIGNMENTS: [Class,Raw]Alignment, not %s', ...
                                                            STAGEofALIGNMENT);
  end

    
    masterTM = subTomoMeta; clear subTomoMeta



if any(classVector{1}(2,:) > 1) || any(classVector{2}(2,:) > 1)
  classSymmetry = true;
else
  classSymmetry = false;
end


% Get the number of tomograms to process.
tomoList = fieldnames(geometry);
nTomograms = length(tomoList);
ctfGroupList = masterTM.('ctfGroupSize');



if (flgClassify)           

[ maskType, maskSize, maskRadius, maskCenter ] = ...
                                  BH_multi_maskCheck(pBH, 'Ali', pixelSize);
  % These are used when 'Cluster' is called, to take the masking parameters
  % from focused PCA/Classification, to produce a montage with reduced
  % Z-dimension & low pass filtering to be used in decision making but not
  % alignment. The X,Y dimensions must not change compared to the full
  % version or else graphical deletion of classes will fail.


  [~, pcaMaskSize, pcaMaskRadius, pcaMaskCenter ] = ...
                                  BH_multi_maskCheck(pBH, 'Cls', pixelSize);
else
  
[ maskType, maskSize, maskRadius, maskCenter ] = ...
                                  BH_multi_maskCheck(pBH, 'Ali', pixelSize);
  
end

[ sizeWindow, sizeCalc, sizeMask, padWindow, padCalc] = ...
                                       BH_multi_validArea( maskSize, maskRadius, scaleCalcSize )
padREF = [0,0,0;0,0,0];                                     

[ nParProcesses, iterList] = BH_multi_parallelJobs(nTomograms,nGPUs, sizeCalc(1),limitToOne);

origMaskSize = sizeMask;
%%%%% Considering removing doNotTrim and making this the default. Temporarily
%%%%% override here.
% % % %doNotTrim = true; 
% This should be moved into BH_multi_validAra                                    
% if (doNotTrim)
%   
if any( (sizeCalc - sizeWindow) < 0 )
  sizeCalc = BH_multi_iterator( sizeWindow, 'fourier' )
end
 sizeMask = sizeWindow;
 % sizeCALC = sizeWindow;
% Find the next largest size for fft
%[ sizeCalc ] = BH_multi_iterator( sizeCalc, 'fourier' )


         
padCalc = BH_multi_padVal(sizeWindow, sizeCalc);           
padWindow = [0,0,0;0,0,0]

if (flgFinalAvg)
  % Load in the even (low-passed) ref and use this only for a (non-iterative)
  % translational alignment.
  
  % Assuming the mask/window sizes haven't changed since the original average
  % that was used to calc the FSC -- that would break things here.
  imgNAME = sprintf('class_%d_Locations_REF_%s', className, 'EVE');
  weightNAME = sprintf('class_%d_Locations_REF_%s_Wgt', className, 'EVE');

  [ refIMG ] = BH_unStackMontage4d(1, ...
                                   masterTM.(cycleNumber).(imgNAME){1}, ...
                                   masterTM.(cycleNumber).(imgNAME){2},...
                                   sizeWindow);

  [ refWDG ] = BH_unStackMontage4d(1, ...
                                masterTM.(cycleNumber).(weightNAME){1},...
                                masterTM.(cycleNumber).(weightNAME){2},...
                                sizeCalc);
                     
% % % % % % %   [ peakMask] = gather(BH_mask3d('sphere', sizeMask, peakSearch, maskCenter));  
  [ peakMask ]  = gather(EMC_maskShape('sphere', sizeMask, peakSearch, 'gpu', {'shift', maskCenter}));

  peakBinary = single(find(peakMask > 0.01));
                                
end


% % % % % % % interpMask = gather(BH_mask3d('sphere',sizeWindow,(sizeWindow./2)-6, maskCenter));
[ interpMask ]  = gather(EMC_maskShape('sphere', sizeWindow, (sizeWindow./2)-6, 'gpu', {'shift', maskCenter}));

% % % % % % % interpMaskWdg = gather(BH_mask3d('sphere',sizeCalc,(sizeCalc./2), maskCenter));
[ interpMaskWdg ]  = gather(EMC_maskShape('sphere', sizeCalc, (sizeCalc./2), 'gpu', {'shift', maskCenter}));

% interpMask = (interpMask > 0.9);
interpMaskWdg = single(find(interpMaskWdg > 0.9));
avgResults=cell(nParProcesses);
avgTomoResults=cell(nParProcesses);
wgtResults=cell(nParProcesses);
wgtTomoResults=cell(nParProcesses);
extResults=cell(nParProcesses);
geoResults=cell(nParProcesses);
cntResults=cell(nParProcesses);
maxClasses = max(size(classVector{1},2),size(classVector{2},2));


minDIM = 64;
padDIM = max([minDIM,minDIM,minDIM], sizeMask);
padVal = minDIM - sizeMask;
padVal = padVal .* (padVal > 0);
fscPAD = padCalc;%[floor(padVal./2); ceil((padVal)./2)]

if ~(use_v2_SF3D)
  
  delete(gcp('nocreate'));
  parpool(nGPUs);

  % TODO need some way of only sending out the command for one tilt or
  % something to prevent collisions that result in no complete recon.

  tiltNameList = fieldnames(masterTM.mapBackGeometry);
  tiltNameList = tiltNameList(~ismember(tiltNameList,{'tomoName','viewGroups'}));
  wgtList = tomoList;

  for iGPU = 1:nGPUs
    nThisGPU = 0;
    for iParProc = iGPU:nGPUs:nParProcesses
      for iTomo = iterList{iParProc}
        iTilt = masterTM.mapBackGeometry.tomoName.(wgtList{iTomo}).tiltName;
        if (any(ismember(tiltNameList,iTilt)))
          tiltNameList{ismember(tiltNameList,iTilt)} = 'continue';
        else
          wgtList{iTomo} = 'continue';
        end
      end
    end
  end
  wgtList=wgtList(~ismember(wgtList,'continue'));
  maxPerGPU = ceil(length(wgtList)/nGPUs) + 1;


  parfor iGPU = 1:nGPUs
  % for iGPU = 1:nGPUs
    for iParProc = iGPU:maxPerGPU:length(wgtList)
      % Caclulating weights takes up a lot of memory, so do all that are necessary
      % prior to the main loop

        BH_multi_loadOrCalcWeight(masterTM,ctfGroupList,wgtList{iParProc},samplingRate ,...
                                  sizeCalc,geometry,cutPrecision,iGPU);



    end
  end
end

try
  parpool(nParProcesses)
catch
  delete(gcp('nocreate'))
  parpool(nParProcesses)
end


if (flgQualityWeight)
%get the average CCC for calculation of particle quality weighting.

    cccVect = [];
    for iParProc = 1:nParProcesses
      for iTomo = iterList{iParProc}              
        positionList = geometry.(tomoList{iTomo});
        cccVect = [cccVect ; positionList(positionList(:,26)~=-9999 & ...
                                          positionList(:,1) >= cccCutOff,1)];
      end
    end
    
    avgCCC = mean(cccVect);
    maxCCC = max(cccVect);
    median(cccVect)
    
    fprintf('\navgCCC %0.3f\nmaxCCC%0.3f\nminCCC%0.3f\n',avgCCC,maxCCC,min(cccVect));
    if abs(maxCCC) > 1
      error('\n|maxCCC(%3.3f) | > 1, either fix the problem or turn of QualityWeighting.\n', maxCCC);
    end
   
else
  maxCCC = [];
  avgCCC = [];
end

% % Clear all of the GPUs prior to entering the main processing loop
for iGPU = 1:nGPUs
  gpuDevice(iGPU);
end

parVect = 1:nParProcesses;
for iParProc = parVect
% for iParProc = parVect

    % Get the gpuIDX assigned to this process
    gpuIDXList = mod(parVect+nGPUs,nGPUs)+1;
    iGPUidx = gpuIDXList(iParProc);
    gpuDevice(iGPUidx);
    fprintf('parProc %d/%d assigned to GPU %d\n',iParProc,nParProcesses,iGPUidx);
    
    nExtracted_tmp = zeros(maxClasses,2);
    firstLoop = true;
    nIgnored = 0;
    nSubTomosTotal = 0;
    avgVolume_tmp = cell(maxClasses,2);
    avgWedge_tmp  = cell(maxClasses,2);  
    geometry_tmp = geometry;

    

    for iRow = 1:maxClasses
      for iCol = 1:2
        avgVolume_tmp{iRow, iCol} = zeros(sizeMask, 'single');
          avgWedge_tmp{iRow, iCol}= zeros(sizeCalc, 'single');

      end
    end
    

    
  if (flgQualityWeight)
  [cccWeight,~,~,~,~,~] = BH_multi_gridCoordinates(sizeCalc, ...
                                               'Cartesian','GPU',...
                                               {'none'},1,0,1);
   cccWeight = (cccWeight ./ pixelSize).^2;
   


  end
    
  if (eachTomo)
    tomoAvgStack = cell(length(iterList{iParProc}),1);
    tomoWgtStack = cell(length(iterList{iParProc}),1);
    for iTomo = 1:length(iterList{iParProc})
      tomoAvgStack{iTomo} = zeros(sizeMask, 'single');
      tomoWgtStack{iTomo} = zeros(sizeMask, 'single');
    end
  end  
  
  nTomos = 1;
  for iTomo = iterList{iParProc}
    
    peakMask_tmp = gpuArray(peakMask);
    peakBinary_tmp = gpuArray(peakBinary);
    
    interpPad = 0; interpTrim = 0 ;
    if interpOrder == 4
      interpPad = BH_multi_padVal(sizeWindow,2.*sizeWindow);
      interpTrim = BH_multi_padVal(2.*sizeWindow,sizeWindow);
      interpMask_tmpBinary = single(find(interpMask > 0.01));
      interpMask_tmp = interpMask;
      interpMaskWdg_tmp = interpMaskWdg;
    else
      interpMask_tmpBinary = gpuArray(single(find(interpMask > 0.01 )));
      interpMask_tmp = gpuArray(interpMask);
      interpMaskWdg_tmp = (gpuArray(interpMaskWdg));
    end

    if (eachTomo)
      tomoAvg = zeros(sizeMask, 'single', 'gpuArray');
      tomoWgt = zeros(sizeCalc , 'single', 'gpuArray');
      tomoCount = 0;
    end
  
    sprintf('gpu %d working on %d/%d volumes\n',iParProc,iTomo,nTomograms)

    tomoName = tomoList{iTomo};
                                


          
    tiltGeometry = masterTM.tiltGeometry.(tomoList{iTomo});
    tomoNumber = masterTM.mapBackGeometry.tomoName.(tomoList{iTomo}).tomoNumber;
    tiltName   = masterTM.mapBackGeometry.tomoName.(tomoList{iTomo}).tiltName;
    coords = masterTM.mapBackGeometry.(tiltName).coords(tomoNumber,1:4);
    [ binShift ] = [0,0,0];%BH_multi_calcBinShift( coords, samplingRate);

    % Load in the geometry for the tomogram, and get number of subTomos.
    positionList = geometry_tmp.(tomoList{iTomo});
    nSubTomos = sum(any(positionList(:,26:26:26*nPeaks) ~= -9999,2));
    
    nSubTomosTotal = nSubTomosTotal + nSubTomos;

   volumeData = [];
   %fprintf('loading tomo %d\n',iTomo);

   tomoNumber = masterTM.mapBackGeometry.tomoName.(tomoList{iTomo}).tomoNumber;
   tiltName = masterTM.mapBackGeometry.tomoName.(tomoList{iTomo}).tiltName;
   fprintf('Loading tomo %d from tilt %s \n',tomoNumber,tiltName);
   reconCoords = masterTM.mapBackGeometry.(tiltName).coords(tomoNumber,:);

   if (flgCutOutVolumes && ~doCut)
     volumeData = [];
   else

    [ volumeData, reconGeometry ] = BH_multi_loadOrBuild( tomoList{iTomo}, ...
                                     reconCoords, mapBackIter, ...
                                     samplingRate,iGPUidx,reconScaling,loadTomo); 
    
     if ( loadTomo )
       volHeader = struct();
       volHeader.('nX') = size(volumeData,1);
       volHeader.('nY') = size(volumeData,2);
       volHeader.('nZ') = size(volumeData,3);
     else
       volHeader = getHeader(volumeData);
     end
    end

                            

    iTiltName = masterTM.mapBackGeometry.tomoName.(tomoList{iTomo}).tiltName; 
    wgtName = sprintf('cache/%s_bin%d.wgt',iTiltName,samplingRate);
    nCtfGroups = ctfGroupList.(tomoList{1})(1);
     wedgeMask = BH_unStackMontage4d(1:nCtfGroups,wgtName,...
                                      ceil(sqrt(nCtfGroups)).*[1,1],'');     
                                    


    % Work on each class seperately pushing to main memory when finished.
     for iGold = 1:2-flgFinalAvg
       iGold
      size(classVector{iGold},2);
      for iClassPos = 1:size(classVector{iGold},2)
        
  
        iTempParticleODD = zeros(sizeMask,cutPrecision,'gpuArray');
        iTempWedgeODD = zeros(sizeCalc,cutPrecision,'gpuArray');
        
       
          iTempParticleEVE = zeros(sizeMask,cutPrecision,'gpuArray');      
          iTempWedgeEVE = zeros(sizeCalc,cutPrecision,'gpuArray');
         

        iClassIDX = classVector{iGold}(1,iClassPos);
    %     nClassPossible = nClassPossible + ...
    %                                   length(find(positionList(:,26) == iClassIDX));
        if classVector{iGold}(1,iClassPos) 
          % pick just the particular class to work with        
      
          if ( flgEstSNR )
            % When the class is for estimating SNR
            includeList = ( any(abs(positionList(:,1:26:26*nPeaks))  >= cccCutOff,2) & ...
                                positionList(:,10) ==  iClassIDX & ...
                                positionList(:,7)  == iGold );
          else
            % When the class is from statistical analysis
            includeList = ( any(abs(positionList(:,1:26:26*nPeaks)) >= cccCutOff,2)  & ...
                                positionList(:,26) ==  iClassIDX & ...
                                positionList(:,7)  == iGold );
          end

        else
          % if class is 0, pick all non-ignored particles

          includeList = ( any(abs(positionList(:,1:26:26*nPeaks))  >= cccCutOff,2)  & ...
                              any(positionList(:,26:26:26*nPeaks) ~= -9999,2)  & ...
                              positionList(:,7)  == iGold );
 
        end
        % Get the position index for each included particle.

        particleIndex = find(includeList);
        nClassPossible = length(particleIndex);



        for iSubTomo = particleIndex'
         iParticle = [];
         iCCCweight = [];
         iWedgeMask = [];
         % symmetry = classVector{iGold}(2, iClassPos);
          if (classSymmetry) && (classVector{iGold}(2, iClassPos) > 0)
            symmetry = classVector{iGold}(2, iClassPos);
          else
            symmetry = 1;
          end
            
          if ( nPeaks > 1 )
            % Calculate a relative weighting, normalize max score to one
            % and then raise to compressBy factor to downweight lower
            % scores.
            
             
            [ peakWgt, sortedList ] = BH_weightAngCheckPeaks( ...
                                               positionList(iSubTomo,:),...
                                               nPeaks, symmetry,...
                                               iSubTomo, tomoList{iTomo},...
                                               bh_global_ML_compressByFactor,...
                                               bh_global_ML_angleTolerance);        
          % Update any re-ordering or elimination                                 
            positionList(iSubTomo,:) = sortedList;                                 
          else
            peakWgt = 1;
          end


            
        for iPeak = 1:nPeaks
          if peakWgt(iPeak) == -9999
            % Skip this peak
            fprintf('skipping peak %d\n')
            continue
          end
           %Check that the given subTomo is not to be ignored


            % Get position and rotation info, angles stored as e1,e3,e2 as in AV3
            % and PEET. This also makes inplane shifts easier to see.

          center = positionList(iSubTomo,[11:13]+26*(iPeak-1))./samplingRate + binShift;
          angles = positionList(iSubTomo,[17:25]+26*(iPeak-1));
          wdgIDX = positionList(iSubTomo,9);
          
          if (flgFinalAvg)
            angles = reshape(angles,3,3)*oddRot;
          end
          
          TLT = masterTM.('tiltGeometry').(tomoList{iTomo});

          if iPeak == 1
            if (use_v2_SF3D)
              [ iSF3D ] = BH_weightMaskMex(sizeCalc, samplingRate, TLT, ...
                                                                center,reconGeometry);
            else
              iSF3D =  gpuArray(wedgeMask{wdgIDX});
            end
          end
                                                       

          
          if (flgQualityWeight)
            iCCC = positionList(iSubTomo,[1]+26*(iPeak-1));
          
            if iCCC < avgCCC
              % Downweight higher frequency in all subTomos with iCCC below the mean
              iBfactor = -1.*(flgQualityWeight.*(acosd(iCCC) - acosd(avgCCC)))^2;
              iCCCweight = exp(iBfactor.*cccWeight);   
             
            else
           
              iCCCweight=1;
             % iBfactor = -1.*flgQualityWeightAbove.*(acosd(iCCC) - acosd(avgCCC));
             % iCCCweight = exp(iBfactor.*cccWeight);
            end
            
            if ( any(flgFilterDefocus))
              iDef = abs(mean(tiltGeometry(:,15))*10^6);
              iDef = -1.*(flgFilterDefocus(1)*max(iDef-1,0.5))^flgFilterDefocus(2);
              fprintf('Using iDef %f\n',iDef);
              % Frequency is already squared so adjust to match iDef scale
              % factor.
              iCCCweight = iCCCweight.*exp(iDef.*cccWeight.^(flgFilterDefocus(2)/2));
            end
          end

          % Find range to extract, and check for domain error.
	  if (flgCutOutVolumes && ~doCut)
            [ indVAL, padVAL, shiftVAL ] = ...
                          BH_isWindowValid(2*CUTPADDING+sizeWindow, ...
                                            sizeWindow, maskRadius, center);
          else
            [ indVAL, padVAL, shiftVAL ] = ...
                          BH_isWindowValid([volHeader.nX,volHeader.nY,volHeader.nZ], ...
                                            sizeWindow, maskRadius, center);
          end
          

          if ~ischar(indVAL)

            
          if (flgCutOutVolumes && ~doCut)
            try
	            particleOUT_name = sprintf('cache/subtomo_%0.7d_%d.mrc',positionList(iSubTomo,4),iPeak);
              iParticle = gpuArray(getVolume(MRCImage(particleOUT_name),...
                                                        [indVAL(1,1),indVAL(2,1)], ...
                                                        [indVAL(1,2),indVAL(2,2)], ...
                                                        [indVAL(1,3),indVAL(2,3)],'keep'));
            catch
              fprintf('\n\nDid not load cut out volume subTomo %d\n\n',iSubTomo);
              continue;
            end
	  else
    
            if ( loadTomo )
              iParticle = gpuArray(volumeData(indVAL(1,1):indVAL(2,1), ...
                                              indVAL(1,2):indVAL(2,2), ...
                                              indVAL(1,3):indVAL(2,3)));

            else
              iParticle = gpuArray(getVolume(volumeData,[indVAL(1,1),indVAL(2,1)], ...
                                                        [indVAL(1,2),indVAL(2,2)], ...
                                                        [indVAL(1,3),indVAL(2,3)],'keep'));
            end

          end
  

            
            [ iParticle ] = BH_padZeros3d(iParticle,  padVAL(1,1:3), ...
                                          padVAL(2,1:3), 'GPU', cutPrecision);
          
       
          
          if (flgCutOutVolumes && doCut)
            % Test with some generic padding , only to be used on bin 1 at
            % first!!! TODO add a flag to check this.
            
            particleOUT = BH_padZeros3d(gather(iParticle), CUTPADDING.*[1,1,1], ...
                                          CUTPADDING.*[1,1,1], 'cpu', cutPrecision);
            
            particleOUT_name = sprintf('cache/subtomo_%0.7d_%d.mrc',positionList(iSubTomo,4),iPeak);
            positionList(iSubTomo,[11:13]+26*(iPeak-1)) = shiftVAL+CUTPADDING+ceil((sizeWindow+1)./2);
            SAVE_IMG(particleOUT,particleOUT_name,pixelSize);
            particleOUT =[];
          end


          % Each particle in the average is rotate about it's origin and then
          % shifted, and the the oddRot and shift is calculated between
          % half-sets after the final iteration. Combining those transforms is
          % not straightforward, so instead, apply the rotation, but calculate
          % each individual shift.
          iShift = shiftVAL;
          
          if (flgFinalAvg)
            % Only the odd half needs to be searched
            
            iRefIMG = BH_resample3d(refIMG{1},angles', iShift,'Bah', 'GPU', 'forward');
            iRefWdg = BH_resample3d(refWDG{1},angles', [0,0,0],'Bah', 'GPU', 'forward');
            
            [ ref_FT ] = BH_bandLimitCenterNormalize(iRefIMG.*peakMask_tmp, ...
                                                    fftshift(iSF3D), ...
                                                     peakBinary_tmp,padCalc,...
                                                     'single');
                                               
            [ part_FT ] = BH_bandLimitCenterNormalize(iParticle.*peakMask_tmp, ...
                                                      fftshift(iRefWdg), ...
                                                      peakBinary_tmp,padCalc,...
                                                      'single');
                                               
            [ peakCoord ] =  BH_multi_xcf_Translational(part_FT, ...
                                                        conj(ref_FT),...
                                                        peakMask_tmp, peakCOM);
                                                      
                                                      
            fprintf('%d peakShift is %f %f %f\n', gather(iSubTomo), gather(peakCoord));
            iRefIMG = [];
            iRefWdg = [];
            part_FT = [];
            ref_FT  = [];
            iShift = iShift + peakCoord;
          end
          


        
        
            interpM = 'linear';
            interpU = 'GPU';        
     
            if interpOrder ==4 
              fprintf('experimental FFT interp\n');
              iParticle = fftn(BH_padZeros3d(iParticle, interpPad(1,:), ...
                                          interpPad(2,:), 'GPU', cutPrecision));
            
            

    
              % No Interp mask, consider adding when it can be tested with
              % the Fourier interp. Memory requirements are too much to use
              % the fast symmetry approach. Need to shift each individually.
              if (symmetry > 1)
                % I was creating new arrays for each symmetry, fixed that,
                % test now, if it runs delete this.
% % %                tmpPart = zeros(sizeWindow,'single','gpuArray');
% % % 
% % %                 for iSym = 0:360/symmetry:360-360/symmetry
% % %                   Rsym = reshape(angles,3,3) * BH_defineMatrix([iSym,0,0],'Bah','inv');
% % %                  tmpPart = tmpPart + BH_padZeros3d(real(ifftn( ...
% % %                                 BH_resample3d(iParticle, ...
% % %                                     Rsym, iShift, ...             
% % %                                     {'Bah',1,interpM,1}, ...
% % %                                     interpU,'inv'))), ...
% % %                                             interpTrim(1,:),interpTrim(2,:),...
% % %                                             'GPU',cutPrecision);
% % % 
% % %                 end
% % %                 iParticle = tmpPart ; tmpPart = [];
              [ iParticle ] = BH_padZeros3d(real(ifftn( ...
                              BH_resample3d(iParticle, ...
                                  angles, iShift, ...             
                                  {'Bah',symmetry,interpM,1}, ...
                                  interpU,'inv'))), ...
                                          interpTrim(1,:),interpTrim(2,:),...
                                          'GPU',cutPrecision);  
              else
              [ iParticle ] = BH_padZeros3d(real(ifftn( ...
                              BH_resample3d(iParticle, ...
                                  angles, iShift, ...             
                                  {'Bah',1,interpM,1}, ...
                                  interpU,'inv'))), ...
                                          interpTrim(1,:),interpTrim(2,:),...
                                          'GPU',cutPrecision);              


              end
            else
                        %iParticle is already on GPU if it should be.
              [ iParticle ] = gpuArray( ...
                              BH_resample3d(iParticle, ...
                                  angles, iShift, ...             
                                  {'Bah',symmetry,interpM,1,interpMask_tmpBinary}, ...
                                  interpU,'inv'));

              [ iParticle ] = BH_padZeros3d(iParticle, ...
                                        -1.*padWindow(1,:),-1.*padWindow(2,:),...
                                        'GPU',cutPrecision);  
                         
            end
            
            % For now just leave linear interp, but test with spline
            [ iWedgeMask ] = BH_resample3d(iSF3D, ...
                                           angles, [0,0,0], ...
                                           {'Bah',symmetry,'linear', ...
                                           1,interpMaskWdg_tmp}, ...
                                           'GPU','inv');          


            iParticle = iParticle -  mean(double(iParticle(interpMask_tmpBinary)));
            iParticle = iParticle ./  rms(double(iParticle(interpMask_tmpBinary)));            
            iParticle = iParticle .* interpMask_tmp;
            if (flgQualityWeight)
                         
              if (singlePrecision)
                iParticle = real(ifftn(fftn(BH_padZeros3d(iParticle,...
                                 padCalc(1,:),padCalc(2,:),'GPU','singleTaper')).*iCCCweight));
                               
              else
                iParticle = real(ifftn(fftn(BH_padZeros3d(iParticle,...
                                 padCalc(1,:),padCalc(2,:),'GPU','doubleTaper')).*iCCCweight));
              end  
              
              iParticle = iParticle(padCalc(1,1)+1 : end - padCalc(2,1), ...
                                            padCalc(1,2)+1 : end - padCalc(2,2), ...
                                            padCalc(1,3)+1 : end - padCalc(2,3) );
               iWedgeMask = iWedgeMask .* fftshift(iCCCweight);
            end
           
        
          
         
          trimAvg =  mean(double(iParticle(:)));
          if isnan(trimAvg)
           fprintf('SubTomo %d from tomogram %s as Nan mean\n',...
           iSubTomo, tomoList{iTomo});
           iParticle(:,:,:) = 0;
           % Flag the particle as ignored
           positionList(iSubTomo, 26) = -9999;
          else
% % %             iParticle = iParticle -  mean(double(iParticle(:)));
% % %             iParticle = iParticle .*  rms(double(iParticle(:)));

              if positionList(iSubTomo,7) == 1 %%%%|| (flgGold == 0)
                iTempParticleODD = iTempParticleODD + iParticle.*peakWgt(iPeak); 

                  iTempWedgeODD    = iTempWedgeODD    + iWedgeMask.*peakWgt(iPeak);
         
                if (flgEstSNR )
                  for iSnr = [iClassPos:5:25]
                    nExtracted_tmp(iSnr, 1) =  nExtracted_tmp(iSnr, 1) + 1;
                  end
                else
                    nExtracted_tmp(iClassPos, 1) = nExtracted_tmp(iClassPos, 1) + 1;
                end
              elseif positionList(iSubTomo,7) == 2 %%%&& (flgGold)
                iTempParticleEVE = iTempParticleEVE + iParticle.*peakWgt(iPeak); 
                  iTempWedgeEVE    = iTempWedgeEVE    + iWedgeMask.*peakWgt(iPeak);
             
               if (flgEstSNR )
                  for iSnr = [iClassPos:5:25]
                    nExtracted_tmp(iSnr, 2) =  nExtracted_tmp(iSnr, 2) + 1;
                  end
                else
                    nExtracted_tmp(iClassPos, 2) = nExtracted_tmp(iClassPos, 2) + 1;
                end
              else
                error('Position List 7 must be a 1 or a 2')
              end

              if (eachTomo)
                tomoAvg = tomoAvg + iParticle.*peakWgt(iPeak);
                  tomoWgt = tomoWgt + iWedgeMask.*peakWgt(iPeak);
                tomoCount = tomoCount + 1;
              end

          end
          
          else

           fprintf('SubTomo %d from tomogram %s only sampled at %f\n',...
           iSubTomo, tomoList{iTomo}, 1-padVAL);
          
           % Flag the particle as ignored
           positionList(iSubTomo, 26:26:26*nPeaks) = -9999;
           nIgnored = nIgnored + 1;
           peakWgt(1:nPeaks) = -9999;
        
          end
    %       end

        end % loop over peaks
        end % end of the loop over subTomos


        if ( flgEstSNR )
          for iSnr = [iClassPos:5:25]
            avgVolume_tmp{ iSnr, 1} =  avgVolume_tmp{ iSnr, 1} + ...
                                                       gather(iTempParticleODD);
            avgVolume_tmp{ iSnr, 2} =  avgVolume_tmp{ iSnr, 2} + ...
                                                       gather(iTempParticleEVE);
                                                     
            avgWedge_tmp{  iSnr, 1} = avgWedge_tmp{  iSnr, 1} + ...
                                                       gather(iTempWedgeODD);
            avgWedge_tmp{  iSnr, 2} = avgWedge_tmp{  iSnr, 2} + ...
                                                       gather(iTempWedgeEVE);                                                     
          end
        else
          avgVolume_tmp{ iClassPos, 1} =  avgVolume_tmp{ iClassPos, 1} + ...
                                                             gather(iTempParticleODD);
          avgVolume_tmp{ iClassPos, 2} =  avgVolume_tmp{ iClassPos, 2} + ...
                                                             gather(iTempParticleEVE);
            avgWedge_tmp{  iClassPos, 1} = avgWedge_tmp{  iClassPos, 1} + ...
                                                             gather(iTempWedgeODD);
            avgWedge_tmp{  iClassPos, 2} = avgWedge_tmp{  iClassPos, 2} + ...
                                                             gather(iTempWedgeEVE);
        end
      end % end of loop over classes
    end% Update geometry to include information on ignored particles.
    geometry_tmp.(tomoList{iTomo})= positionList;

    if (eachTomo)
      tomoAvgStack{nTomos} = gather(tomoAvg./tomoCount);
      tomoWgtStack{nTomos} = gather(tomoWgt);
      nTomos = nTomos +1;
    end
    

    
  end % end of the loop over Tomograms
   

  avgResults{iParProc} = avgVolume_tmp;
  extResults{iParProc} = nExtracted_tmp;
  geoResults{iParProc} = geometry_tmp;
  cntResults{iParProc} = [nSubTomosTotal, nIgnored];
  wgtResults{iParProc} = avgWedge_tmp;
  

  
  if (eachTomo)
    avgTomoResults{iParProc} = tomoAvgStack;
    wgtTomoResults{iParProc} = tomoWgtStack;
  end
  
end % end parfor loop on gpus


      
nExtracted = zeros(maxClasses,2);
avgVolume = cell(maxClasses,2);
avgWedge  = cell(maxClasses,4); 
nSubTomosTotal = 0;
nIgnored = 0;
 

for iParProc = 1:nParProcesses  

  nExtracted = nExtracted + gather(extResults{iParProc});
  nSubTomosTotal = nSubTomosTotal + cntResults{iParProc}(1,1);
  nIgnored  = nIgnored  + cntResults{iParProc}(1,2);
  
  for iVol = 1:maxClasses
    iVol
    for iHalf = 1:2-flgFinalAvg
      iHalf
      if iParProc == 1
       avgVolume{iVol, iHalf} = avgResults{iParProc}{iVol, iHalf};
       avgWedge{iVol, iHalf} = wgtResults{iParProc}{iVol, iHalf};
 
      else       
        avgVolume{iVol, iHalf} = avgVolume{iVol, iHalf} + ...
                                                  avgResults{iParProc}{iVol, iHalf};
        avgWedge{iVol, iHalf} = avgWedge{iVol, iHalf} + ...
                                                  wgtResults{iParProc}{iVol, iHalf};
                                               
      end
    end
  end
  
  for iTomo = iterList{iParProc}
    geometry.(tomoList{iTomo}) = geoResults{iParProc}.(tomoList{iTomo});
  end
end



fprintf('\n%d / %d subTomos extracted\n',sum(nExtracted(:)), nSubTomosTotal)
                                       



masterTM.(cycleNumber).('nSubTomoAveraged') = gather(sum(nExtracted(:)));
masterTM.(cycleNumber).(sprintf('newIgnored_Avg%s',fieldPrefix)) = ...
                                                               gather(nIgnored);

                                                           
% get the total class average by combining eve/odd
classStorage = cell(maxClasses,2);
if (doNotTrim) && (flgClassify)
  filteredClass= cell(maxClasses,2);
  % low-pass to see class averages more clearly.
% %   [ bandpassFilt ] = BH_bandpass3d( sizeMask, 0.2, 300, 30, 'GPU',pixelSize);
end


       
% % % if (doNotTrim) && (flgClassify)
% % %   % reduce z dimension to area focused on in classification
% % %   % taken from PCA mask values.
% % %   zCenter = (sizeMask(3)+1)./2 + pcaMaskCenter(3);
% % %   zLow = floor(zCenter - pcaMaskRadius(3));
% % %   zTop = zLow + 2.* pcaMaskRadius(3);
% % %   if zLow < 1
% % %     fprintf('setting Z-low from %d to 1\n', zLow);
% % %     zLow = 1;
% % %   end
% % %   if zTop > sizeMask(3) %%%size(filteredClass{iClassPos,iGold},3)
% % %     fprintf('setting Z-top from %d to sizeFiltAvg,3\n', zTop);
% % %     zTop = sizeMask(3); %%%size(filteredClass{iClassPos,iGold},3);
% % %   end
% % %   zLow = zLow;
% % %   zTop = zTop ; % note not 2nd indx, just adding a shift
% % %   sizeFilteredClass = [sizeMask(1),sizeMask(2),( zTop - zLow +1)]
% % % end

for iClass = 1:maxClasses
  classStorage{iClass,1} = zeros(sizeMask, 'single');
  classStorage{iClass,2} = zeros(sizeMask, 'single');
% % %   if (doNotTrim) && (flgClassify)
% % %     filteredClass{iClass,1} = zeros(sizeFilteredClass, 'single');
% % %     filteredClass{iClass,2} = zeros(sizeFilteredClass, 'single');
% % %   end
end

if (eachTomo)
  [o,c] = system('mkdir -p initialTomoAvgs');
  % sizeWeight mask is sizeMask or 128^3 whichever is larger
%   bandpassFiltTomo = BH_bandpass3d( sizeCalc, lpTomo(1), lpTomo(2), lpTomo(3), 'GPU',pixelSize); 

  for iParProc = 1:nParProcesses
    nTomos = 1;
    
    for iTomo = iterList{iParProc}
      tomoName = sprintf('initialTomoAvgs/%s.mrc',tomoList{iTomo});
      classTmp = avgTomoResults{iParProc}{nTomos};
      avgTomoResults{iParProc}{nTomos} = [];

      SAVE_IMG(classTmp,tomoName,pixelSize);
      clear classTmp
      nTomos = nTomos + 1; 
    end
  end
end

classSum = cell(2,1);
classWgtSum = cell(2,1);
if (saveClassSum > -1)
  classSum{1} = zeros(size(avgVolume{1,1}),'single');
  classSum{2} = zeros(size(avgVolume{1,2}),'single');
  classWgtSum{1} = zeros(size(avgWedge{1,1}),'single');
  classWgtSum{2} = zeros(size(avgWedge{1,2}),'single');
end

for iClassPos = 1:maxClasses

  
  
    if (doNotTrim) && (flgClassify)
% % % % % % %       m = BH_mask3d('sphere',sizeMask,floor(sizeMask./2-6),pcaMaskCenter); 
      [ m ]  = EMC_maskShape('sphere', sizeMask,floor(sizeMask./2-6), 'gpu', {'shift', pcaMaskCenter});
    

    else
      m = 1;
    end
    
    % Re-weight both halves whether flgGold or not.
    
    for iGold = 1:2-flgFinalAvg

       avgVolume{iClassPos,iGold} = avgVolume{iClassPos,iGold} ./ ...
                                                sum(nExtracted(iClassPos,iGold));

       
    end
    
    fprintf('flgGold = %d\n',flgGold);
    for iGold = 1:2-flgFinalAvg

% % %       if( flgGold )
        if iGold == 1
          halfSet = 'ODD';
        else
          halfSet = 'EVE';
        end
% % %       else
% % %         halfSet = 'STD';
% % %       end
     

        classStorage{iClassPos,iGold} = gather(avgVolume{iClassPos,iGold} );


    
      if isnan( mean(classStorage{iClassPos,iGold}(:)) )
        clear classAVG 
        fprintf('zeroing out classavg because of NaN values detected.\n')
      else
% % %         if (doNotTrim) && (flgClassify)
% % %           % lowpass according to Kms bandpass
% % %           bandpassFilt = BH_bandpass3d( sizeMask, 0.1,300,30, 'GPU',pixelSize);
% % %           
% % %           tmpFilt = BH_bandLimitCenterNormalize( ...
% % %                                              classStorage{iClassPos,iGold}.* ...
% % %                                              m, bandpassFilt, (m>0.95), ...
% % %                                              [0,0,0;0,0,0],'single');
% % %           tmpFilt = m.*real(ifftn(tmpFilt));
% % %           tmpFilt = tmpFilt(:,:,zLow:zTop  );  
% % %           tmpFilt = tmpFilt - mean(tmpFilt(:));
% % %           tmpFilt = tmpFilt ./ rms(tmpFilt(:));
% % %                           
% % %           filteredClass{iClassPos,iGold} = gather(tmpFilt);
% % %         end

        % Normalize the regular averages
        classStorage{iClassPos,iGold} = classStorage{iClassPos,iGold} - ...
                                         mean(classStorage{iClassPos,iGold}(:));
                                
        classStorage{iClassPos,iGold} = classStorage{iClassPos,iGold} ./ ...
          rms(classStorage{iClassPos,iGold}(:));
        classStorage{iClassPos,iGold} = gather(classStorage{iClassPos,iGold});
        
        if (saveClassSum > -1)
          fprintf('adding class %d to %d\n',iClassPos,iGold);
          classSum{iGold} = classSum{iGold} + classStorage{iClassPos,iGold};
          classWgtSum{iGold} = classWgtSum{iGold} + avgWedge{iClassPos,iGold};
        end
      end
      
      
    
     
    end
    
end




% Using the filtered class average, calc real space CCC to reorder the even
% class to [most likely] match the corresponding odd class.
classListOut = 0;
% % % if (flgGold) && strcmpi(STAGEofALIGNMENT, 'Cluster')
if  strcmpi(STAGEofALIGNMENT, 'Cluster') && (flgClassify ~= -1)
% % %   [classListOut, geometry] = reorder_classes(filteredClass(:,1),filteredClass(:,2),maxClasses, geometry);

% % %   filteredClass(:,2) = filteredClass(classListOut(:,2), 2);
  [classListOut, geometry] = reorder_classes(avgVolume(:,1),avgVolume(:,2),maxClasses, geometry);

  avgVolume(:,2) = avgVolume(classListOut(:,2), 2);
  avgWedge(:,2) = avgWedge(classListOut(:,2), 2);
  
  classMatches = fopen(sprintf('%s_class%d_%s_matches.txt', ...
                                outputPrefix,className,fieldPrefix),'w');
  fprintf(classMatches,'Score\t\tOdd\tEVE\n');                            
  fprintf(classMatches,'%d\t%d\t%2.6f\n', classListOut');                            
  fclose(classMatches);
  
  masterTM.(cycleNumber).(sprintf('class_%d_%s_EveOddIdx',className,fieldPrefix)) = classListOut;
  

end

% Second option allows re-use of class designations to generate
% multi-reference alignment.
if strcmpi(STAGEofALIGNMENT, 'Cluster')
  masterTM.(cycleNumber).(ClusterGeomNAME) = geometry;
elseif strcmpi(STAGEofALIGNMENT, 'RawAlignment') && flgMultiRefAlignment && ~flgClassify
  masterTM.(cycleNumber).('ClusterRefGeom') = geometry;
else
  masterTM.(cycleNumber).('Avg_geometry') = geometry;
end
% Should this save differently depending on the stage of alignment?? I
% think so but leave alone for nw.

  
for iGold = 1:2-flgFinalAvg
%%%      if( flgGold )
        if iGold == 1
          halfSet = 'ODD';
        else
          halfSet = 'EVE';
        end
%%%      else
%%%        halfSet = 'STD';
%%%      end
 imgCounts = gather([classVector{iGold}(1,:) ; nExtracted(:,iGold)']);

      
      [montOUT, imgLocations] = BH_montage4d(avgVolume(:,iGold), '');
      
      
      imout = sprintf('%s_class%d_%s_%s_NoWgt.mrc',outputPrefix, ...
                                         className, fieldPrefix, halfSet);
      classOut = sprintf('class_%d_Locations_%s_%s_NoWgt', className,fieldPrefix, halfSet);                        

      masterTM.(cycleNumber).(classOut) = {imout,imgLocations,imgCounts};
      
      if (flgFinalAvg)
        system(sprintf('mv %s preHalfSetAli_%s',imout,imout));
      end
      SAVE_IMG(montOUT, imout, pixelSize);
      %%%%%%%%
      [montOUT, imgLocations] = BH_montage4d(avgWedge(:,iGold), '');
     
      imout = sprintf('%s_class%d_%s_%s_Wgt.mrc',outputPrefix, ...
                                         className, fieldPrefix, halfSet);
      classOut = sprintf('class_%d_Locations_%s_%s_Wgt', className,fieldPrefix, halfSet);                        

      % For the weight, instead of imgCounts save the padValues
      masterTM.(cycleNumber).(classOut) = {imout,imgLocations,fscPAD};
      
      if (flgFinalAvg)
        system(sprintf('mv %s preHalfSetAli_%s',imout,imout));
      end
      SAVE_IMG(montOUT, imout);
      %%%%%%%%
% % %       if (doNotTrim) && (flgClassify)
% % %         [montOUT, ~] = BH_montage4d(filteredClass(:,iGold), '');
% % %         imout = sprintf('%s_filtered%d_%s_%s.mrc',outputPrefix, ...
% % %                                          className, fieldPrefix, halfSet);
% % %         SAVE_IMG(MRCImage(gather(montOUT)), imout,pixelSize);
% % %         
% % %       end
      
      if (saveClassSum > -1)
        imgCounts = gather([classVector{iGold}(1,:) ; nExtracted(:,iGold)']);


        [montOUT, imgLocations] = BH_montage4d(classSum(iGold), '');

        imout = sprintf('%s_class%d_%s_%s_NoWgt.mrc',outputPrefix, ...
                                           saveClassSum, 'Raw', halfSet);
        classOut = sprintf('class_%d_Locations_%s_%s_NoWgt', saveClassSum,'Raw', halfSet);                        
        SAVE_IMG(montOUT, imout,pixelSize);
        masterTM.(cycleNumber).(classOut) = {imout,imgLocations,imgCounts};

        [montOUT, imgLocations] = BH_montage4d(classWgtSum(iGold), '');

        imout = sprintf('%s_class%d_%s_%s_Wgt.mrc',outputPrefix, ...
                                           saveClassSum, 'Raw', halfSet);
        classOut = sprintf('class_%d_Locations_%s_%s_Wgt', saveClassSum,'Raw', halfSet);
        SAVE_IMG(montOUT, imout,pixelSize);
        masterTM.(cycleNumber).(classOut) = {imout,imgLocations,imgCounts};

      end
end



    masterTM = gather(masterTM);
    classVector = gather(classVector);
    nExtracted = gather(nExtracted);      
      


subTomoMeta = gather(masterTM);
subTomoMeta.(cycleNumber).('SymmetryApplied').(STAGEofALIGNMENT).(sprintf('%s_className',fieldPrefix)) = ...
                                                            gather(classVector);
cycleNumber = gather(cycleNumber);                                                            
subTomoMeta.('currentCycle') = gather(CYCLE);


save(pBH.('subTomoMeta'), 'subTomoMeta'); 



fprintf('Total execution time : %f seconds\n', etime(clock, startTime));


% clean everything up, since this function is called from other functions.

delete(gcp('nocreate'));
for iGPU = 1:nGPUs
  gpuDevice(gpuList(iGPU));
end

% The final condition is setting up RefAlignment which with classification
% would have been averaging with 'Cluster' == Stageof
if strcmpi(STAGEofALIGNMENT, 'RawAlignment') || ...
   strcmpi(STAGEofALIGNMENT, 'NoAlignment')  || ...
   strcmpi(STAGEofALIGNMENT, 'RefAlignment') 
  BH_fscGold_class(PARAMETER_FILE, num2str(CYCLE), STAGEofALIGNMENT); 
elseif strcmpi(STAGEofALIGNMENT, 'Alignment') 
  BH_fscGold_class(PARAMETER_FILE, num2str(CYCLE), 'RawAlignment');
end


try
  parpool(nGPUs)
catch
  delete(gcp('nocreate'));
  parpool(nGPUs)
end

for iGPU = 1:nGPUs
  gpuDevice(gpuList(iGPU));
end




if ~( flgEstSNR )
    load(sprintf('%s.mat', pBH.('subTomoMeta')), 'subTomoMeta');
    masterTM = subTomoMeta;
    %%%%%%%%%%%%%55 Reweight now that the FSC is calculated

    for iGold = 1:2
  %%%    if( flgGold )
        if iGold == 1
          halfSet = 'ODD';
        else
          halfSet = 'EVE';
        end
   %%%   else
  %%%      halfSet = 'STD';
  %%%    end
        imgIN = sprintf('class_%d_Locations_%s_%s_NoWgt', ...
                                           className, fieldPrefix, halfSet);
        wgtIN = sprintf('class_%d_Locations_%s_%s_Wgt', ...
                                           className, fieldPrefix, halfSet); 

      [ refIMG{iGold} ] = BH_unStackMontage4d(1:maxClasses, ...
                                     masterTM.(cycleNumber).(imgIN){1}, ...
                                     masterTM.(cycleNumber).(imgIN){2},...
                                     sizeWindow);

      [ refWGT{iGold} ] = BH_unStackMontage4d(1:maxClasses, ...
                                  masterTM.(cycleNumber).(wgtIN){1},...
                                  masterTM.(cycleNumber).(wgtIN){2},...
                                  sizeCalc);
    end



  % This is slow ass when using cones and class averages and wouldn't be too
  % hard to put into parallel. Do that once the next manuscript is finished.
% % % %     if (flgMultiRefAlignment && ~flgClassify && ~strcmpi(STAGEofALIGNMENT, 'RefAlignment'))
    if (~flgMultiRefAlignment && ~flgClassify && ~strcmpi(STAGEofALIGNMENT, 'RefAlignment'))

	      nClassesReWgt = 1;
    else
        nClassesReWgt = maxClasses;
    end

    for iRef = 1:nClassesReWgt
      fprintf('Stage of alignment %s\niRef %d\n',STAGEofALIGNMENT,iRef);
      if strcmpi(STAGEofALIGNMENT, 'RawAlignment') || strcmpi(STAGEofALIGNMENT, 'NoAlignment')   
        savePrefix = 'Raw'
      else
        savePrefix = fieldPrefix
      end

      % When extracting class averages only a single FSC is initially available for 
      % calculation of the SPW filter. 
      if strcmpi(STAGEofALIGNMENT, 'Cluster')
        iRefPrev = 1;
      else
        iRefPrev = iRef;
      end

      if (flgGold) || (flgClassify < 0)
        flgCombine = 0;
        flgRefCutOff = 1;
      else
        flgCombine = 1;
        flgRefCutOff = 0;
      end





      % When switching from no Classification to classification, there will be no
      % previous savePrefix ref, it will be Raw. Another reason to get rid of the
      % (artificial) distinction. For now, use a try catch .

      try
        fscParams = masterTM.(cycleNumber).('fitFSC').(sprintf('%s%d',savePrefix,iRefPrev));
        aliParams = masterTM.(cycleNumber).('fitFSC').(sprintf('Resample%s%d',savePrefix,iRefPrev));
        mskParams = masterTM.(cycleNumber).('fitFSC').(sprintf('Mask%s%d',savePrefix,iRefPrev));
      catch
        fprintf('\nReverting from %s to Raw in loading fitFSC\n',savePrefix);
        fscParams = masterTM.(cycleNumber).('fitFSC').(sprintf('%s%d','Raw',iRefPrev));
        aliParams = masterTM.(cycleNumber).('fitFSC').(sprintf('Resample%s%d','Raw',iRefPrev));
        mskParams = masterTM.(cycleNumber).('fitFSC').(sprintf('Mask%s%d','Raw',iRefPrev));
      end

  %     if (flgGold) && strcmpi(STAGEofALIGNMENT, 'Cluster')   
  %       iOdd = cccFinal(iRef,2);
  %       iEve = cccFinal(iRef,3);
  %     else
        iOdd = iRef;
        iEve = iRef;
  %     end

      if (flgFinalAvg)
        % negative to combine but NOT apply the xform to the odd set
        flgCombine = -1;
        flgRefCutoff = 0;
   
      end

      % Only send the lowest Bfactor if not flgFinalAvg
      if (flgFinalAvg)
        bFactorSend = pBH.('Fsc_bfactor');
      else
        bFactorSend = pBH.('Fsc_bfactor')(1);
      end
      
      refTMP = gather(BH_multi_cRef_Vnorm(fscParams, aliParams, mskParams,...
                          {refIMG{1}{iOdd},refIMG{2}{iEve}}, ...
                          {refWGT{1}{iOdd},refWGT{2}{iEve}}, ...
                          flgCombine,flgRefCutOff, pixelSize, bFactorSend));

      if ~(flgFinalAvg)
        refIMG{1}{iOdd} = refTMP{1,1};
        if (flgGold) || (flgClassify < 0)
          refIMG{2}{iEve} = refTMP{1,2}; 
        end
        clear refTMP
      end

    end


    for iGold = 1:2-flgFinalAvg
      if( flgGold ) || (flgClassify < 0)
        if iGold == 1
          halfSet = 'ODD';
        else
          halfSet = 'EVE';
        end
      else
        halfSet = 'STD';
      end
      iRef
      imgIN = sprintf('class_%d_Locations_%s_%s_NoWgt', ...
                                           className, fieldPrefix, halfSet); 
      imgCounts = masterTM.(cycleNumber).(imgIN){3};
      % Save the unweighted, weighted imgs, weightes, optionally filtered.

      if (flgFinalAvg)
        for iBfactor = 1:length(pBH.('Fsc_bfactor'))
          imout = sprintf('%s_class%d_%s_bFact-%d.mrc',outputPrefix, ...
                               className, 'final',pBH.('Fsc_bfactor')(iBfactor));

          SAVE_IMG(refTMP{iBfactor}, imout, pixelSize);
        end
      else

        [montOUT, imgLocations] = BH_montage4d(refIMG{iGold}(:), '');
        imout = sprintf('%s_class%d_%s_%s.mrc',outputPrefix, ...
                                           className, fieldPrefix, halfSet);
        classOut = sprintf('class_%d_Locations_%s_%s', className,fieldPrefix, halfSet);                        

        masterTM.(cycleNumber).(classOut) = {imout,imgLocations,imgCounts};
        SAVE_IMG(montOUT, imout, pixelSize);
      end
      %%%%%%%
    end


  subTomoMeta = masterTM;
  subTomoMeta.('CUTPADDING') = CUTPADDING;
  if (flgCutOutVolumes && doCut)
    subTomoMeta.('volumesAreCutOut') = 1;
  end
  save(pBH.('subTomoMeta'), 'subTomoMeta');
  
end


end % end of average3d function

function [classListOut, geometry] = reorder_classes(odd,eve,maxClasses, geometry)


classListOut = zeros(maxClasses,3);
classRemain = 1:maxClasses;

oddNorm = zeros(1,maxClasses);
eveNorm = zeros(1,maxClasses);

for iClass = 1:maxClasses
  odd{iClass} = odd{iClass} - mean(odd{iClass}(:));
  oddNorm(iClass) = (numel(odd{iClass}).*rms(odd{iClass}(:)));
  eveNorm(iClass) = rms(eve{iClass}(:));
end

for iOdd = 1:maxClasses  
    ccc = zeros(1,length(classRemain));
    nCCC = 1;
    for iEve = classRemain
      ccc(nCCC) = sum(sum(sum(odd{iOdd}.*eve{iEve})))/(oddNorm(iOdd).*eveNorm(iEve));
      nCCC = nCCC + 1;
    end
    
    [v,c] = max(ccc);
    classListOut(iOdd,:) = [iOdd,classRemain(c),v];
    classRemain(c) = 0;
    classRemain = classRemain(classRemain ~= 0);   
end




tomoList = fieldnames(geometry);
for iTomo = 1:length(tomoList)
  % select included even half
  includeList = (geometry.(tomoList{iTomo})(:,26) ~= -9999) & ...
                 geometry.(tomoList{iTomo})(:,7) == 2; 
              
  % Get positions for all classes
  classIDX = zeros(length(includeList),maxClasses);
  for iClass = 1:maxClasses
    classIDX(:,iClass) = includeList & geometry.(tomoList{iTomo})(:,26) == iClass;
  end
  
  for iClass = 1:maxClasses
    classListOut(iClass,2)
    newClass = classIDX(:,classListOut(iClass,2));
    geometry.(tomoList{iTomo})(logical(newClass),26) = iClass;

  end

end

end
