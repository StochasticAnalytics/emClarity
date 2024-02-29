function [  ] = BH_pcaPub(PARAMETER_FILE, CYCLE, PREVIOUS_PCA)
%Extract and interpolate a subTomogram from a 3d volume.
%
%   Input variables:
%
%  CYCLE = 0,1,2 etc.
%
%
%
%   samplingRate = Binning factor, assumed to be integer value. Image is first
%              smoothed by an appropriate low-pass filter to reduce aliasing.
%
%   emc.Pca_randSubset = -1, count all non-ignored particles (class -9999).
%
%                float, randomly select this many paparticleBandpassrticles for the
%                decomposition, denote by updating the flag in column 8 to be 1
%                or 0.
%
%                string - indicates a mat file with prior decomposition.
%
%   emc.Pca_maxEigs = Maximum number of principle components to save, general 50
%                   has been plenty. This is a big memory saver.
%
%   bandpass =  [HIGH_THRESH, HIGH_CUT, LOW_CUT, PIXEL_SIZE]
%
%   AVERAGE_MOTIF = string pointing to an average to be used for wedge masked
%                   differences.
%
%   GEOMETRY = A structure with tomogram names as the field names, and geometry
%              information in a 28 column array.
%              Additionally, a field called 'source_path' has a value with the
%              absolute path to the location of the tomograms.
%
%              The input is a string 'Geometry_templatematching.mat' for
%              example, and it is expected that the structure is saved as the
%              variable named geometry.
%
%
%
%
%   Output variables:
%
%   None = files are written to disk in the current directory. This will be a
%   mat file that has the
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Goals & Limitations:
%
%   Extract, filter, and mask included subtomograms and perform a singular value
%   decomposition. Either all included (non -9999 class) or a subset may be
%   specified. If a subset is specified, then the decomposition will be done on
%   these while the full projection follows.
%
%   Assumed to run on GPU (extraction) and then CPU (Decomposition) the latter
%   can be very memory intensive.
%
%   To avoid errors and complications, the full data set will be extracted,
%   which is used to generate an average
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
%     - calculate max size possible to keep temp data matrix on the gpu,
%     sigfificantly reduces overhead.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETER_FILE = 'testParam.m';
%CYCLE = 2;
%PREVIOUS_PCA = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (nargin ~= 3)
  %  error('PARAMETER_FILE, CYCLE, PREVIOUS_PCA')
end


% FIXME:
% To test seeding the classification with existing classes, rather than always reverting to the global average,
% use the mechanism in place to handle  multiple references at different length scales derived from the global average,
% to instead be used for multiple distinct classes. If the results are promising, then expand so each ref may also be
% looked at over its own scale space.
test_multi_ref_diffmap = true;

startTime =  datetime("now");

CYCLE = EMC_str2double(CYCLE);
PREVIOUS_PCA = EMC_str2double(PREVIOUS_PCA);

global bh_global_binary_pcaMask_threshold;

use_new_interpolator = true;
% Previous_pca has two functions, when true and > 0 use the decomposition
% calculated from a random subset of the data to project the full data set
% onto each of the selected principle components. when true and < 0 run
% with the same parameters as before, but this time load in the adjusted
% variance maps to use as a mask for each scale space.
% -2 use variance map, -1 use stdDev instead

switch PREVIOUS_PCA
  case -3
    flgVarianceMap = 0;
    flgStdDev = 0;
    flgLoadMask = 1;
    PREVIOUS_PCA = 0;
  case -2
    flgVarianceMap = 1;
    flgStdDev = 1;
    PREVIOUS_PCA = 0;
    flgLoadMask = 0;
  case -1
    flgVarianceMap = 1;
    flgStdDev = 0.5;
    PREVIOUS_PCA = 0;
    flgLoadMask = 0;
  case 0
    flgVarianceMap = 0;
    flgStdDev = 0;
    flgLoadMask = 0;
  case 1
    % case 0 and 1 same except value of PREVIOUS_PCA
    flgVarianceMap = 0;
    flgStdDev = 0;
    flgLoadMask = 0;
  otherwise
    error('PREVIOUS_PCA should be 1,0,-1,-2,-3')
end

flgWMDs = 3;

cycleNumber = sprintf('cycle%0.3u', CYCLE);

emc = BH_parseParameterFile(PARAMETER_FILE);
reconScaling = 1;
%%% Put this in the param file later - the input values should be in angstrom
%%% and are the relevant scale spaces for classification.


samplingRate   = emc.('Cls_samplingRate');
refSamplingRate= emc.('Ali_samplingRate');

pixelSize = emc.pixel_size_angstroms .* samplingRate;
refPixelSize = emc.pixel_size_angstroms .* refSamplingRate;

% FIMXE: Probably remove this incomplete idea
if (refSamplingRate ~= samplingRate)
  error('refSamplingRate ~= samplingRate')
end

nCores  = BH_multi_parallelWorkers(emc.('nCpuCores'));
pInfo = parcluster();

nParticlesOnGpuBeforeTransfer = emc.('PcaGpuPull');


outputPrefix   = sprintf('%s_%s', cycleNumber, emc.('subTomoMeta'));
%%%flgGold      = emc.('flgGoldStandard');


flgNorm = 1;% emc.('flgNormalizeWMDs');
try
  flgPcaShapeMask = emc.('flgPcaShapeMask');
catch
  flgPcaShapeMask = 1;
end

% The defaults used in fscGold are modified here to make a more permissive
% mask since we are concerned with densities that are likely damped during
% averaging due to low occupancy.
emc.shape_mask_lowpass = emc.shape_mask_lowpass + 10;
emc.shape_mask_threshold = emc.shape_mask_threshold - 0.4;


try
  tmpVal = emc.('whitenPS');
  if (numel(tmpVal) == 3)
    wiener_constant = tmpVal(3);
  else
    error('flgWhitenPS should be a 3 element vector');
  end
catch
  wiener_constant = 0.0;
end



% Removed flgGold everywhere else, but keep ability to classify full data set at
% the end (after all alignment is finished.)


if (emc.classification)
  flgGold = 0;
else
  flgGold = 1;
end


load(sprintf('%s.mat', emc.('subTomoMeta')), 'subTomoMeta');
mapBackIter = subTomoMeta.currentTomoCPR;


if (test_multi_ref_diffmap && ~emc.multi_reference_alignment)
  test_multi_ref_diffmap = false;
  fprintf("WARNING: test_multi_ref_diffmap is incompatible with ~emc.multi_reference_alignment, disabling\n");
end

geom_name=''
if (emc.multi_reference_alignment )
  geom_name='ClusterClsGeom';
  
else
  geom_name='Avg_geometry';
end
geometry = subTomoMeta.(cycleNumber).(geom_name);


try
  CUTPADDING = subTomoMeta.('CUTPADDING')
catch
  CUTPADDING=20
end



[ useGPU ] = BH_multi_checkGPU( -1 );
gDev = gpuDevice(useGPU);


% Get the number of tomograms to process.
tomoList = fieldnames(geometry);
nTomograms = length(tomoList);


[ maskType, maskSize, maskRadius, maskCenter ] = ...
  BH_multi_maskCheck(emc, 'Cls', pixelSize);

[ preMaskType, preMaskSize, preMaskRadius, preMaskCenter ] = ...
  BH_multi_maskCheck(emc, 'Ali', refPixelSize);

% This is prob not a good way to make sure the mask size matches::w
maskSize=preMaskSize;

% Make sure everthing matches the extracted average and wedge
cpuVols = struct;

[ preSizeWindow, preSizeCalc, preSizeMask, prePadWindow, prePadCalc ] = ...
  BH_multi_validArea(preMaskSize,preMaskRadius, emc.scale_calc_size )


[ sizeWindow, sizeCalc, sizeMask, padWindow, padCalc ] = ...
  BH_multi_validArea( maskSize, maskRadius, emc.scale_calc_size )


if (test_multi_ref_diffmap)
  refName = emc.('Raw_className');
else
  refName = 0;
end

% If emc.classification is negative combine the data for clustering, but don't set
% any of the alignment changes to be persistant so that extracted class
% averages are still independent half-sets.
if (flgGold)
  oddRot = eye(3);
else
  iRefPrev = 1;
  
  aliParams = subTomoMeta.(cycleNumber).('fitFSC').(sprintf('Resample%s%d','Ref',iRefPrev))
  oddRot = reshape(aliParams(1,:),3,3)';
  % refine the translation per particle.

  clear iRefPrev
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Added for test_multi_ref_diffmap %%
refVector = cell(2,1);
refGroup = cell(2,1);
classVector{1}  = emc.('Raw_classes_odd')(1,:);
classVector{2}  = emc.('Raw_classes_eve')(1,:);
refVectorFull{1}= [emc.('Raw_classes_odd');classVector{1} ]
refVectorFull{2}= [emc.('Raw_classes_eve');classVector{2} ]
for iGold = 1:2
  % Sort low to high, because order is rearranged as such unstack
  refVectorFull{iGold} = sortrows(refVectorFull{iGold}', 1)';
  % class id corresponding to membership in ???_refName
  refVector{iGold} = refVectorFull{iGold}(1,:)
  % reference id, so multiple classes can be merged into one
  refGroup{iGold}  = refVectorFull{iGold}(3,:)
end


% make sure the number of references match the unique groups in the classVector
% and also that the class/group pairs match the class/ref pairs.
nReferences(1:2) = [length(unique(refGroup{1})),length(unique(refGroup{1}))];
nReferences = nReferences .* [~isempty(refGroup{1}),~isempty(refGroup{2})];

averageMotif = cell(2,1);

if (test_multi_ref_diffmap)
  fprintf('emc.n_scale_spaces = %d\n',emc.n_scale_spaces);
  fprintf('nReferences = %d\n',nReferences);
  emc.n_scale_spaces = nReferences(1);
  pause(3);
else
  nReferences = [1,1];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iGold = 1:2
  
  if iGold == 1;
    halfSet = 'ODD';
  else
    halfSet = 'EVE';
  end
  
  
  
  
  imgNAME = sprintf('class_%d_Locations_Ref_%s', refName, halfSet);
  
  
  [ averageMotif{iGold} ] = BH_unStackMontage4d(1:nReferences(iGold), ...
    subTomoMeta.(cycleNumber).(imgNAME){1}, ...
    subTomoMeta.(cycleNumber).(imgNAME){2},...
    preSizeMask);
  
  
  
  if ~(test_multi_ref_diffmap)
    averageMotif{iGold} = averageMotif{iGold}{1};
  end
  if (flgLoadMask) && (iGold == 1)
    fprintf('\n\nLoading external mask\n');
    externalMask = OPEN_IMG('single',sprintf('%s-pcaMask',subTomoMeta.(cycleNumber).(imgNAME){1}));
  end
end

% IF combining for analysis, resample prior to any possible binning.
if ~(flgGold)
  if (nReferences(1) ~= nReferences(2))
    error('When combining half sets, the number of references must match')
  end
  size(averageMotif)
  if (nReferences(1) > 1)
    for iRef = 1:nReferences(1)
      averageMotif{1}{iRef} = averageMotif{2}{iRef} + ...
        BH_resample3d(gather(averageMotif{1}{iRef}), ...
        oddRot, ...
        aliParams(2,1:3), ...
        {'Bah',1,'spline'}, 'cpu', ...
        'forward');
      averageMotif{2}{iRef} = [];
    end
  else
    averageMotif{1} = averageMotif{2} + ...
      BH_resample3d(gather(averageMotif{1}), ...
      oddRot, ...
      aliParams(2,1:3), ...
      {'Bah',1,'spline'}, 'cpu', ...
      'forward');
    averageMotif{2} = [];
  end

end

%%% incomplete, the idea is to generate an antialiased scaled volume for PCA
if ( refSamplingRate ~= samplingRate )
  fprintf('Resampling from %d refSampling to %d pcaSampling\n',refSamplingRate,samplingRate);
  for iGold = 1:1+flgGold
    averageMotif{iGold} = BH_reScale3d(averageMotif{iGold},'',sprintf('%f',1/samplingRate),'GPU');
  end
  
  if (flgLoadMask)
    externalMask = BH_reScale3d(externalMask,'',sprintf('%f',1/samplingRate),'GPU');
  end
end



prevVarianceMaps = struct();
if (flgVarianceMap)
  
  for iGold = 1:1+flgGold
    
    if (flgGold)
      if iGold == 1;
        halfSet = 'ODD';
      else
        halfSet = 'EVE';
      end
    else
      halfSet = 'STD';
    end
    
    % For randomsubset (PREVIOUS_PCA = 0) the suffix is *_pcaPart.mat) but
    % presumably we could have also just done full, so try that first
    
    try
      load(sprintf('%s_%s_pcaFull.mat',outputPrefix,halfSet))
    catch
      fprintf('\nDid not find, %s_%s_pcaFUll.mat, trying *_pcaPart.mat\n',outputPrefix,halfSet);
      load(sprintf('%s_%s_pcaPart.mat',outputPrefix,halfSet));
    end
    
    % In most cases, this is the number of "features" specified in the
    % parameter file, but in some data not even this may non-zero singluar
    % values are found, so the number could be different (lower)
    for iScale = 1:emc.n_scale_spaces
      eigsFound = size(coeffs{iScale},1);
      fname = sprintf('%s_varianceMap%d-%s-%d.mrc', ...
        outputPrefix, eigsFound, halfSet, iScale);
      
      prevVarianceMaps.(sprintf('h%d',iGold)).(sprintf('s%d',iScale)) = ...
        OPEN_IMG('single', fname).^flgStdDev;
    end
    clear v coeffs eigsFound idxList
  end
end


if (PREVIOUS_PCA)
  volumeMask = gpuArray(OPEN_IMG('single', sprintf('%s_pcaVolMask.mrc',outputPrefix)));
else
  
  if (emc.Pca_constrain_symmetry)
    gridSearch = eulerSearch(emc.symmetry,180,5,360,5,0.0,1,true);
    [ volumeMask ]    = BH_mask3d(maskType, sizeMask, maskRadius, maskCenter, ...
      '3d', gridSearch.number_of_asymmetric_units);
  else
    [ volumeMask ]    = BH_mask3d(maskType, sizeMask, maskRadius, maskCenter);
  end
  
  if ( flgPcaShapeMask )
    % For testing we won't handle this block
    if (test_multi_ref_diffmap)
      error('test_multi_ref_diffmap is incompatible with flgPcaShapeMask')
    end
    % when combining the addition is harmless, but is a convenient way to
    % include when sets are left 100% separate.
    %       volumeMask = volumeMask .* BH_mask3d(averageMotif{1}+averageMotif{1+flgGold}, pixelSize, '','');
    volumeMask = volumeMask .* EMC_maskReference(averageMotif{1}+averageMotif{1+flgGold}, pixelSize, ...
      {'pca', true; 'lowpass', emc.shape_mask_lowpass; 'threshold', emc.shape_mask_threshold});
    
  end
  
  if (flgLoadMask)
    % For testing we won't handle this block
    if (test_multi_ref_diffmap)
      error('test_multi_ref_diffmap is incompatible with flgLoadMask')
    end
    volumeMask = volumeMask .* externalMask;
  end
  
  SAVE_IMG(MRCImage(gather(volumeMask)),sprintf('%s_pcaVolMask.mrc',outputPrefix),pixelSize);
end

volMask = struct();
nPixels = zeros(2,emc.n_scale_spaces);
for iScale = 1:emc.n_scale_spaces
  for iGold = 1:1+flgGold
    stHALF = sprintf('h%d',iGold);
    stSCALE = sprintf('s%d',iScale);
    if (flgVarianceMap)
      if (test_multi_ref_diffmap)
        error('test_multi_ref_diffmap is incompatible with flgVarianceMap')
      end
      volTMP = gather(volumeMask.*prevVarianceMaps.(stHALF).(stSCALE));
    else
      volTMP = gather(volumeMask);
    end
    
    masks.('volMask').(stHALF).(stSCALE) = (volTMP);
    masks.('binary').(stHALF).(stSCALE)  = (volTMP >= bh_global_binary_pcaMask_threshold);
    masks.('binary').(stHALF).(stSCALE)  = ...
      masks.('binary').(stHALF).(stSCALE)(:);
    masks.('binaryApply').(stHALF).(stSCALE)  = (volTMP >= 0.01);
    
    nPixels(iGold,iScale) = gather(sum(masks.('binary').(stHALF).(stSCALE)));
    clear volTMP stHALF stSCALE
  end
end
clear volumeMask



% radius, convert Ang to pix , denom = equiv stdv from normal to include, e.g.
% for 95% use 1/sig = 1/2
%stdDev = 1/2 .* (emc.pca_scale_spaces ./ pixelSize - 1)  .* 3.0./log(emc.pca_scale_spaces)
if ~(test_multi_ref_diffmap)
  threeSigma = 1/3 .* (emc.pca_scale_spaces ./ pixelSize);
  for iScale = 1:emc.n_scale_spaces
    
    kernelSize = ceil(threeSigma(iScale).*3) + 3;
    kernelSize = kernelSize + (1-mod(kernelSize,2));
    % masks.('scaleMask').(sprintf('s%d',iScale))  = EMC_gaussianKernel([1,kernelSize],  threeSigma(iScale), 'cpu', {});
    masks.('scaleMask').(sprintf('s%d',iScale))  = EMC_gaussianKernel([1,1,1].*kernelSize, 2*threeSigma(iScale), 'gpu', {});

  
    % SAVE_IMG( masks.('scaleMask').(sprintf('s%d',iScale)), ...
    %   sprintf('%s_scaleMask_s%d.mrc',outputPrefix,iScale),pixelSize);
    
    
  end
end

avgMotif_FT = cell(1+flgGold,emc.n_scale_spaces);
avgFiltered = cell(1+flgGold,emc.n_scale_spaces);
% Here always read in both, combine if flgGold = 0
for iGold = 1:1+flgGold
  for iScale = 1:emc.n_scale_spaces
    
    if (test_multi_ref_diffmap)
      tmp_avg = averageMotif{iGold}{iScale};
    else
      tmp_avg = averageMotif{iGold};
    end
    
    
    tmp_avg = tmp_avg - mean(tmp_avg(masks.('binaryApply').(sprintf('h%d',iGold)).(sprintf('s%d',iScale))));
    tmp_avg = tmp_avg ./ rms(tmp_avg(masks.('binaryApply').(sprintf('h%d',iGold)).(sprintf('s%d',iScale))));
    tmp_avg = tmp_avg .* masks.('volMask').(sprintf('h%d',iGold)).(sprintf('s%d',iScale));
    % FIXME: ideally we would do both, but for testing I am stealing scaleSpace for iRef
    if ~(test_multi_ref_diffmap)
      tmp_avg = EMC_convn(single(gpuArray(tmp_avg)) , single(gpuArray(masks.('scaleMask').(sprintf('s%d',iScale))) ));
    end
    avgMotif_FT{iGold, iScale} = ...
      BH_bandLimitCenterNormalize(tmp_avg,...
      BH_bandpass3d(sizeMask,1e-6,400,2.2*pixelSize,'GPU',pixelSize), ...
      masks.('binaryApply').(sprintf('h%d',iGold)).(sprintf('s%d',iScale)),...
      [0,0,0;0,0,0],'single');

    avgFiltered{iGold, iScale} = real(ifftn(avgMotif_FT{iGold, iScale}));
    
    avgFiltered{iGold, iScale} = avgFiltered{iGold, iScale} - mean(avgFiltered{iGold, iScale}(masks.('binary').(sprintf('h%d',iGold)).(sprintf('s%d',iScale))));
    avgFiltered{iGold, iScale} = gather(avgFiltered{iGold, iScale} ./rms(avgFiltered{iGold, iScale}(masks.('binaryApply').(sprintf('h%d',iGold)).(sprintf('s%d',iScale)))) .* ...
      masks.('binaryApply').(sprintf('h%d',iGold)).(sprintf('s%d',iScale)));
  end
end



montOUT = BH_montage4d(avgFiltered(1,:),'');
SAVE_IMG(MRCImage(montOUT), sprintf('%s_filt.mrc', outputPrefix),pixelSize);
clear montOUT


% If emc.Pca_randSubset is string with a previous matfile use this, without any
% decomposition.

for iGold = 1:1+flgGold
  flgRefIsPadded = 0;
  if (flgGold)
    if iGold == 1;
      halfSet = 'ODD';
      stHALF = sprintf('h%d',iGold);
      randSet =1;
    else
      stHALF = sprintf('h%d',iGold);
      halfSet = 'EVE';
      randSet = 2;
    end
  else
    stHALF = sprintf('h%d',iGold);
    halfSet = 'STD';
    randSet = [1,2];
  end
  
  if (PREVIOUS_PCA)
    previousPCA = sprintf('%s_%s_pcaPart.mat',outputPrefix,halfSet);
    emc.Pca_randSubset = -1;
    [ geometry, nTOTAL, nSUBSET ] = BH_randomSubset( geometry,'pca', -1 , randSet);
  else
    if (emc.Pca_randSubset)
      previousPCA = false;
      [ geometry, nTOTAL, nSUBSET ] = BH_randomSubset( geometry,'pca', emc.Pca_randSubset, randSet );
    else
      previousPCA = false;
      [ geometry, nTOTAL, nSUBSET ] = BH_randomSubset( geometry,'pca', -1 , randSet);
    end
  end
  
  % Extend the random subset to each peak if needed
  % FIXME: this should just be in random subset
  if (emc.nPeaks > 1)
    for iTomo = 1:nTomograms
      selectedList = geometry.(tomoList{iTomo})(:,8) > 0;
      geometry.(tomoList{iTomo})(selectedList,8+26:26:emc.nPeaks*26) = 1;
      clear selectedList
    end
    nTOTAL = nTOTAL*emc.nPeaks;
    nSUBSET = nSUBSET*emc.nPeaks;
  end
  
  % Initialize array in main memory for pca
  clear dataMatrix tempDataMatrix
  dataMatrix = cell(3,1);
  tempDataMatrix = cell(3,1);
  for iScale = 1:emc.n_scale_spaces
    dataMatrix{iScale} = zeros(nPixels(iGold,iScale), nSUBSET, 'single');
    tempDataMatrix{iScale} = zeros(nPixels(iGold,iScale), nParticlesOnGpuBeforeTransfer, 'single', 'gpuArray');
  end
  
  % Pull masks onto GPU (which are cleared along with everything else when
  % the device is reset at the end of each loop.)
  gpuMasks = struct();
  
  for iScale = 1:emc.n_scale_spaces
    stSCALE = sprintf('s%d',iScale);
    
    gpuMasks.('volMask').(stSCALE) = ...
      gpuArray(masks.('volMask').(stHALF).(stSCALE));
    gpuMasks.('binary').(stSCALE) = ...
      gpuArray(masks.('binary').(stHALF).(stSCALE));
    gpuMasks.('binaryApply').(stSCALE)  = ...
      gpuArray(masks.('binaryApply').(stHALF).(stSCALE));
    if ~(test_multi_ref_diffmap)
      gpuMasks.('scaleMask').(stSCALE) = gpuArray(masks.('scaleMask').(stSCALE));
    end
    
    
    gpuMasks.('highPass').(stSCALE) = BH_bandpass3d(sizeMask,1e-6,400,2.2*pixelSize,'GPU',pixelSize);
  end
  
  % % %   for iGold_inner = 1:1+flgGold
  % % %     for iScale = 1:emc.n_scale_spaces
  % % %       avgMotif_FT{iGold_inner, iScale} = ...
  % % %                gpuArray(cpuVols.('avgMotif_FT').(sprintf('g%d_%d',iGold_inner,iScale)));
  % % %     end
  % % %   end
  
  
  nExtracted = 1;
  nTemp = 1;
  nTempPrev = 0;
  idxList = zeros(1,nSUBSET);
  peakList = zeros(1,nSUBSET);
  
  firstLoop = true;  sI = 1;
  nIgnored = 0;
  for iTomo = 1:nTomograms
    
    
    tomoName = tomoList{iTomo};

    iGPU = 1;
    tomoIdx = subTomoMeta.mapBackGeometry.tomoName.(tomoList{iTomo}).tomoIdx;
    tiltName = subTomoMeta.mapBackGeometry.tomoName.(tomoList{iTomo}).tiltName;
    reconCoords = subTomoMeta.mapBackGeometry.tomoCoords.(tomoList{iTomo});

    TLT = subTomoMeta.('tiltGeometry').(tomoList{iTomo});
    
    
    if (emc.flgCutOutVolumes)
      volumeData = [];
    else
      do_load = false;
      [ volumeData ] = BH_multi_loadOrBuild(tomoList{iTomo}, ...
                                            mapBackIter, ...
                                            samplingRate, ...
                                            iGPU, ...
                                            do_load);
      volHeader = getHeader(volumeData);
    end
    
    
    iTiltName = subTomoMeta.mapBackGeometry.tomoName.(tomoName).tiltName;
    
    tiltGeometry = subTomoMeta.tiltGeometry.(tomoList{iTomo});
    
    fprintf('Working on %d/%d volumes %s\n',iTomo,nTomograms,tomoName);
    % Load in the geometry for the tomogram, and get number of subTomos.
    positionList = geometry.(tomoList{iTomo});
    
    % Loop over peaks inside each tomo to limit wedge mask xfer
    positionList = positionList(positionList(:,26) ~= -9999,:);
    nSubTomos = size(positionList,1);
    
    

    
    % reset for each tomogram
    wdgIDX = 0;
    radialMask = '';
    if (flgNorm)
      
      
      % bins = 1./[1000,800,600,400,300,200,150,100,80,60,50,40,35,30,28,26,24,22,20,18,16,14,12,10,8,6,4,2];
      % bins = [0, bins];
      
      [radialGrid,~,~,~,~,~] = BH_multi_gridCoordinates(size(avgMotif_FT{iGold, iScale}),'Cartesian',...
        'GPU',{'none'},1,0,1);
      
      bins = radialGrid(1:floor(size(avgMotif_FT{iGold, iScale},1)/2),1,1);
      bins = bins(bins < 0.5);
      
      radialMask = cell(length(bins)-1,1);
      
      for iBin = 1:length(bins)-1
        radialMask{iBin} = find(radialGrid >= bins(iBin) & radialGrid < bins(iBin+1));
      end
      
      radialGrid = '';
      
    end
    
    wdgBP = ifftshift(gpuMasks.('highPass').(sprintf('s%d',iScale)));
    for iSubTomo = 1:nSubTomos
      %%%%% %%%%%
      
      
      % Check that the given subTomo is not to be ignored - for now, treat
      % all peaks as included. The assumption is that using this will be
      % for initializing the project to get a good starting model. "True"
      % classification will be done at a later stage after reducing to some
      % subset of peaks. FIXME
      includeParticle = positionList(iSubTomo, 8);
      
      iPeak=0; % make sure this exists if we are no including the particle
      if (includeParticle)
        make_sf3d = true;
        for iPeak = 0:emc.nPeaks-1
          
          % Get position and rotation info, angles stored as e1,e3,e2 as in AV3
          % and PEET. This also makes inplane shifts easier to see.
          
          center = positionList(iSubTomo,[11:13]+26*iPeak)./samplingRate;
          angles = positionList(iSubTomo,[17:25]+26*iPeak);
          
          
          if ( make_sf3d )
            make_sf3d = false;
            radialGrid = '';
            padWdg = [0,0,0;0,0,0];
            [ wedgeMask ] = BH_weightMaskMex(sizeWindow, samplingRate, TLT, center, reconCoords, wiener_constant);
          end
          
          % If flgGold there is no change, otherwise temporarily resample the
          % eve halfset to minimize differences due to orientaiton
          if positionList(iSubTomo,7) == 1 % This is true for all peaks
            % TODO FIXME should this be the transpose of oddRot?
            angles = reshape(angles,3,3) * oddRot;
          end
          wedgeMask = wedgeMask .* wdgBP;
          
          % Find range to extract, and check for domain error.
          if (emc.flgCutOutVolumes)
            [ indVAL, padVAL, shiftVAL ] = ...
              BH_isWindowValid(2*CUTPADDING+sizeWindow, ...
              sizeWindow, maskRadius, center);
          else
            [ indVAL, padVAL, shiftVAL ] = ...
              BH_isWindowValid([volHeader.nX,volHeader.nY,volHeader.nZ], ...
              sizeWindow, maskRadius, center);
          end
          
          if ~(flgGold)
            shiftVAL = shiftVAL + aliParams(2,1:3)./samplingRate;
          end
          
          using_this_subtomo = true;
          particleIDX = positionList(iSubTomo, 4); % Same for all peaks
          if ~ischar(indVAL)
            % Read in and interpolate at single precision as the local values
            % in the interpolant suffer from any significant round off errors.
            
            
            
            if (emc.flgCutOutVolumes)
              
              particleOUT_name = sprintf('cache/subtomo_%0.7d_%d.mrc',positionList(iSubTomo,4),iPeak+1);
              iParticle = gpuArray(OPEN_IMG('single',particleOUT_name,...
                [indVAL(1,1),indVAL(2,1)], ...
                [indVAL(1,2),indVAL(2,2)], ...
                [indVAL(1,3),indVAL(2,3)],'keep'));
              
            else
              
              iParticle = gpuArray(OPEN_IMG('single', volumeData,[indVAL(1,1),indVAL(2,1)], ...
                [indVAL(1,2),indVAL(2,2)], ...
                [indVAL(1,3),indVAL(2,3)],'keep'));
            end
            
            
            if any(padVAL(:))
              [ iParticle ] = BH_padZeros3d(iParticle,  padVAL(1,1:3), ...
                padVAL(2,1:3), 'GPU', 'single');
            end
            
            if ( use_new_interpolator )
              % Pulling in the newer interpolater from alignRaw3d_v2. There, I instantiate a new interpolator every subtomo, but it is generally
              % being used many times, over the angle loop. It may be more efficient to do this outside the for subtomo loop here, but
              % to start, just do it the same way.
              use_only_once = true;
              [ ~, iParticle ] = interpolator(gpuArray(iParticle),angles, shiftVAL, 'Bah', 'inv', emc.symmetry, use_only_once);
              
              [ ~, iWedge ] = interpolator(gpuArray(wedgeMask),angles,[0,0,0], 'Bah', 'inv', emc.symmetry, use_only_once);
              
            else
              % Transform the particle, and then trim to motif size
              
              [ iParticle ] = BH_resample3d(iParticle, angles, shiftVAL, ...
                'Bah', 'GPU', 'inv');
              
              [ iWedge ] = BH_resample3d(wedgeMask, angles, [0,0,0], ...
                'Bah', 'GPU', 'inv');
            end
            
            
            
            iTrimParticle = iParticle(padWindow(1,1)+1 : end - padWindow(2,1), ...
              padWindow(1,2)+1 : end - padWindow(2,2), ...
              padWindow(1,3)+1 : end - padWindow(2,3));
            
            
            
            
            for iScale = 1:emc.n_scale_spaces
              if (test_multi_ref_diffmap)
                iPrt = iTrimParticle;
              else
                iPrt = EMC_convn(iTrimParticle , gpuMasks.('scaleMask').(sprintf('s%d',iScale)));
              end
              
              iPrt = BH_bandLimitCenterNormalize( ...
                iPrt .* ...
                gpuMasks.('volMask').(sprintf('s%d',iScale)), ...
                gpuMasks.('highPass').(sprintf('s%d',iScale)),...
                gpuMasks.('binary').(sprintf('s%d',iScale)),...
                [0,0,0;0,0,0],'single');
              
              [iWmd,~] = BH_diffMap(avgMotif_FT{iGold, iScale},iPrt,ifftshift(iWedge),...
                flgNorm,pixelSize,radialMask, padWdg);
              
              
              % using_this_subtomo is set true for each particle, but if we are not finite for ANY scale space, we want to skip this particle for all
              % scale spaces.
              if all(isfinite(iWmd(gpuMasks.('binary').(sprintf('s%d',iScale)))))
                tempDataMatrix{iScale}(:,nTemp) = single(iWmd(gpuMasks.('binary').(sprintf('s%d',iScale))));
              else
                fprintf('inf or nan in subtomo %d scalePace %d',particleIDX, iScale);
                using_this_subtomo = false;
              end
              
            end % loop on scale spaces
            clear iAvg iWmd  iTrimParticle
            
            
            
            if (using_this_subtomo)
              idxList(1, nExtracted) = particleIDX;
              peakList(1,nExtracted) = iPeak+1;
              nExtracted = nExtracted + 1;
              nTemp = nTemp + 1;
              
              % pull data of the gpu every 1000 particls (adjust this to max mem)
              if nTemp - 1 == nParticlesOnGpuBeforeTransfer
                for iScale = 1:emc.n_scale_spaces
                  dataMatrix{iScale}(:,1+nTempPrev:nTemp+nTempPrev-1) = ...
                    gather(tempDataMatrix{iScale}(:,1:nTemp-1));
                end
                
                nTempPrev = nTempPrev + nTemp - 1;
                nTemp = 1;
              end
            else
              nIgnored = nIgnored + 1;
              fprintf('Ignoring subtomo %d from %s\n',particleIDX, tomoList{iTomo});
              eraseIDX = subTomoMeta.(cycleNumber).(geom_name).(tomoList{iTomo})(:,4) == particleIDX;
              subTomoMeta.(cycleNumber).(geom_name).(tomoList{iTomo})(eraseIDX, 26+iPeak*26) = -9999;
            end
            
            
          else
            nIgnored = nIgnored + 1;
            fprintf('Ignoring inside subtomo %d from %s\n',particleIDX, tomoList{iTomo});
            eraseIDX = subTomoMeta.(cycleNumber).(geom_name).(tomoList{iTomo})(:,4) == particleIDX;
            subTomoMeta.(cycleNumber).(geom_name).(tomoList{iTomo})(eraseIDX, 26+iPeak*26) = -9999;
            
          end % end of update data matrix or if ignoring update metaData
          
        end % end of loop over peaks
      else
        nIgnored = nIgnored + 1;
        fprintf('Ignoring outside subtomo %d from %s\n',particleIDX, tomoList{iTomo});
        eraseIDX = subTomoMeta.(cycleNumber).(geom_name).(tomoList{iTomo})(:,4) == particleIDX;
        subTomoMeta.(cycleNumber).(geom_name).(tomoList{iTomo})(eraseIDX, 26+iPeak*26) = -9999;
      end % end of ignore if statment from extracted window out of bounds (ischar(indVAL))
      if ~rem(iSubTomo,100)
        fprintf('\nworking on %d/%d subTomo peak %d/%d from %d/%d Tomo\n', ...
          iSubTomo, nSubTomos,iPeak+1,emc.nPeaks, iTomo,nTomograms);
        
        fprintf('Total nExtracted = %d\n', nExtracted-1);
        fprintf('Total nIgnored = %d\n', nIgnored);
        
      end
    end % end of the loop over subTomos
    
    clear volumeData
  end % end of the loop over Tomograms,
  
  % % %   volBinaryMask = reshape(gather(volBinaryMask),sizeMask);
  for iScale = 1:emc.n_scale_spaces
    masks.('binary').(stHALF).(sprintf('s%d',iScale)) = ...
      reshape(masks.('binary').(stHALF).(sprintf('s%d',iScale)),sizeMask);
  end
  
  
  subTomoMeta.(cycleNumber).('newIgnored_PCA').(halfSet) = gather(nIgnored);
  
  subTomoMeta = subTomoMeta;
  save(sprintf('%s.mat', emc.('subTomoMeta')), 'subTomoMeta');
  
  for iScale = 1:emc.n_scale_spaces
    dataMatrix{iScale}(:,1+nTempPrev:nTemp-1+nTempPrev) = ...
      gather(tempDataMatrix{iScale}(:,1:nTemp-1));
  end
  
  clear tempDataMatrix
  % Get rid of any zero vals from newly ignored particles which are there due to
  % pre-allocation. Assuming no zeros have found their way in anywhere else which
  % would be a major problem.
  cleanIDX = (idxList~=0)';
  idxList = idxList(cleanIDX);
  peakList = peakList(cleanIDX);
  
  for iScale = 1:emc.n_scale_spaces
    dataMatrix{iScale} = dataMatrix{iScale}(:,1:size(idxList,2));
    % Center the rows
    for row = 1:size(dataMatrix{iScale},1)
      dataMatrix{iScale}(row,:) = dataMatrix{iScale}(row,:) -  mean(double(dataMatrix{iScale}(row,:)));
    end
  end
  
  
  
  
  %save('preparpoolSave.mat');
  try
    EMC_parpool(nCores);
  catch
    delete(gcp('nocreate'));
    EMC_parpool(nCores);
  end
  
  if (previousPCA)
    % Read the matrix of eigenvectors from the prior PCA.
    oldPca = load(previousPCA);
    U = oldPca.U;
    clear oldPca;
    sDiag = cell(emc.n_scale_spaces,1);
    coeffs = cell(emc.n_scale_spaces,1);
    
    for iScale = 1:emc.n_scale_spaces
      % Sanity checks on the dimensionality
      numEigs = size(U{iScale}, 2);
      if nPixels(iGold,iScale) ~= size(U{iScale}, 1)
        error('Image size %d does not match that of previous PCA %d!', nPixels(iGold,iScale), size(U{iScale},1));
      end
      
      coeffs{iScale} = U{iScale}' * dataMatrix{iScale};
    end
  else
    
    U = cell(emc.n_scale_spaces,1);
    V = cell(emc.n_scale_spaces,1);
    S = cell(emc.n_scale_spaces,1);
    sDiag = cell(emc.n_scale_spaces,1);
    coeffs = cell(emc.n_scale_spaces,1);
    varianceMap = cell(emc.n_scale_spaces,1);
    for iScale = 1:emc.n_scale_spaces
      
      krylovScalar = 5;
      % Calculate the decomposition
      [ U{iScale},S{iScale},V{iScale}, convergenceFlag ] = svds(double(dataMatrix{iScale}), ...
        emc.Pca_maxEigs, 'largest', ...
        'MaxIterations',1000, ... % default 300
        'SubspaceDimension',max(krylovScalar*emc.Pca_maxEigs,30),... % default max(3*emc.Pca_maxEigs,15)
        'Display',true); % Diagnostics default false (will this work in compiled?)
      U{iScale} = single(U{iScale});
      S{iScale} = single(S{iScale});
      V{iScale} = single(V{iScale});
      
      %             [U{iScale},S{iScale},V{iScale}] = svd(dataMatrix{iScale}, 0);
      
      sDiag{iScale} = diag(S{iScale});
      numNonZero = find(( sDiag{iScale} ~= 0 ), 1, 'last');
      % For Method 1, save eigenvectors 1-4 (or user-specified max) as images
      eigsFound = min(emc.Pca_maxEigs, numNonZero);
      
      fprintf('Found %d / %d non-zero eigenvalues sum = %4.4f, in set %s.\n All singular values converged is t/f ( %d ) ', ...
        numNonZero, size(S{iScale}, 1), sum(sDiag{iScale}), halfSet, convergenceFlag);
      
      coeffs{iScale} = S{iScale} * V{iScale}'
      
      % Can be GB-TB if calculated full
      %varianceMap{iScale} = (U{iScale}*S{iScale}.^2*V{iScale} ./ numel(U{iScale}-1));
      fprintf('Size S, %d %d  Size U %d %d \n', size(S{iScale},1),size(S{iScale},2), size(U{iScale},1),size(U{iScale},2));
      
      % We want the diagnol of US^2U'/ n-1
      % This will be emc.Pca_maxEigs * Nvoxels matrix (U is Nvoxels * emc.Pca_maxEigs)
      rightSide = S{iScale}(1:numNonZero,1:numNonZero).^2*U{iScale}';
      varianceMap = zeros(nPixels(iGold,iScale),1);
      for k = 1:nPixels(iGold,iScale)
        %         varianceMap(k) = U{iScale}(k,1)*rightSide(1,k) + ...
        %                          U{iScale}(k,2)*rightSide(2,k) + ...
        %                          U{iScale}(k,3)*rightSide(3,k);
        varianceMap(k) = U{iScale}(k,:)*rightSide(:,k);
      end
      varianceMap = varianceMap ./ (numel(varianceMap) - 1);
      
      tmpReshape = zeros(prod(sizeMask),1);
      tmpReshape(masks.('binary').(stHALF).(sprintf('s%d',iScale)) ) = varianceMap(:);
      
      fname = sprintf('%s_varianceMap%d-%s-%d.mrc',outputPrefix, eigsFound, halfSet, iScale);
      SAVE_IMG(MRCImage(single(gather(reshape(tmpReshape, sizeMask)))), fname,pixelSize);
      
      
      
      eigList = cell(eigsFound,1);
      eigList_SUM = cell(eigsFound,1);
      
      for iEig = 1:eigsFound
        tmpReshape = zeros(prod(sizeMask),1);
        tmpReshape(masks.('binary').(stHALF).(sprintf('s%d',iScale)) ) = U{iScale}(:, iEig);
        eigenImage = reshape(tmpReshape, sizeMask);
        eigenImage = eigenImage - mean(eigenImage(masks.('binaryApply').(stHALF).(sprintf('s%d',iScale))));
        eigenImage = eigenImage ./rms(eigenImage(masks.('binaryApply').(stHALF).(sprintf('s%d',iScale)))).* masks.('binary').(stHALF).(sprintf('s%d',iScale)) ;
        eigList{iEig,1} = gather(eigenImage);
        eigList_SUM{iEig,1} = gather((eigenImage + avgFiltered{iGold, iScale} )./2);
      end
      
      
      [ eigMont ] = BH_montage4d(eigList, 'eigMont');
      [ eigMont_SUM ] = BH_montage4d(eigList_SUM, 'eigMont_SUM');
      fname = sprintf('%s_eigenImage%d-%s-mont_%d.mrc',outputPrefix, eigsFound, halfSet, iScale);
      fname_SUM = sprintf('%s_eigenImage%d-SUM-%s-mont_%d.mrc',outputPrefix, eigsFound, halfSet, iScale);
      SAVE_IMG(MRCImage(single(gather(eigMont))), fname,pixelSize);
      SAVE_IMG(MRCImage(single(gather(eigMont_SUM))), fname_SUM,pixelSize);
      
      
      % If requested, limit the number of principal components and coeffs saved
      if emc.Pca_maxEigs < size(S{iScale}, 1)
        fprintf('Saving only the first %d principal components.\n',          ...
          emc.Pca_maxEigs);
        if ~isempty(U{iScale}) % U will not exist for pcaMethods 2 or 3
          U{iScale} = U{iScale}(:, 1:emc.Pca_maxEigs);
        end
        if ~(previousPCA)
          S{iScale} = S{iScale}(1:emc.Pca_maxEigs, 1:emc.Pca_maxEigs);
          V{iScale} = V{iScale}(:, 1:emc.Pca_maxEigs);
        end
        coeffs{iScale} = coeffs{iScale}(1:emc.Pca_maxEigs, :); %
      end
    end
  end
  
  
  % Only U is needed for further analysis, so save only this, unless
  % troubleshooting.
  if (previousPCA)
    for iScale = 1:emc.n_scale_spaces
      
      % If requested, limit the number of principal components and coeffs saved.
      if emc.Pca_maxEigs < size(U{iScale}, 2)
        fprintf('Saving only the first %d principal components.\n',          ...
          p.pcaMaxNumComponents);
        U{iScale} = U{iScale}(:, 1:emc.Pca_maxEigs);
        coeffs{iScale} = coeffs{iScale}(1:emc.Pca_maxEigs, :);
      end
      
    end
    save(sprintf('%s_%s_pcaFull.mat',outputPrefix,halfSet), 'nTOTAL', 'coeffs','idxList','peakList', 'sDiag');
  else
    if (emc.Pca_randSubset)
      save(sprintf('%s_%s_pcaPart.mat',outputPrefix,halfSet),'U', 'idxList','peakList');
    else
      save(sprintf('%s_%s_pcaFull.mat',outputPrefix,halfSet), 'nTOTAL', 'coeffs','idxList','peakList', 'sDiag');
    end
  end
  
  fprintf('Total execution time on %s set: %f seconds\n', halfSet, seconds(datetime("now") - startTime));
  
  close all force;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  clear dataMatrix U S V coeffs eigMont eigMontSum
  
  delete(gcp('nocreate'));
  
  % after resetting the device, bring back masks etc.
  
  
end % end of loop over halfsets
gpuDevice(1);
delete(gcp('nocreate'));
end % end of pca function


