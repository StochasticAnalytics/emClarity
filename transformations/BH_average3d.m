
function [  ] = BH_average3d(PARAMETER_FILE, CYCLE, STAGEofALIGNMENT)



% FIXME: hacking in a test
test_fuzz=false;
test_multi_ref_diffmap= true;


if (nargin ~= 3)
  error('args = PARAMETER_FILE, CYCLE, STAGEofALIGNMENT')
end


startTime =  datetime("now");
CYCLE = EMC_str2double(CYCLE);

if strcmpi(STAGEofALIGNMENT, 'RawAlignment')
  % Ensure we don't have any duplicates: TODO: add an override flag
  % This modifies the RawAlign geometry, so should be cycle -1
  if (CYCLE > 0)
    BH_removeDuplicates(PARAMETER_FILE,sprintf('%d', CYCLE-1));
  end
end

cycleNumber = sprintf('cycle%0.3u', CYCLE);

emc = BH_parseParameterFile(PARAMETER_FILE);
load(sprintf('%s.mat', emc.('subTomoMeta')), 'subTomoMeta');


mapBackIter = subTomoMeta.currentTomoCPR;

if (CYCLE == 0)
  fprintf('No quality weighting in the initial cycle after template matching\n');
  emc.flgQualityWeight = 0;
end

if (emc.projectVolumes && ~emc.flgCutOutVolumes)
  emc.flgCutOutVolumes = true;
end

volumesNeedToBeExtracted = 0;
if (emc.flgCutOutVolumes)
  if isfield(subTomoMeta,'volumesAreCutOut')
    if ~(subTomoMeta.volumesAreCutOut)
      volumesNeedToBeExtracted = 1;
    end
  else
    volumesNeedToBeExtracted = 1;
  end
end

rotConvention = 'Bah';
if ( emc.doHelical )
  rotConvention = 'Helical';
end

% The weights are only re-estimated for an out of plane search. Until this
% happens, they are not valid.
if (emc.track_stats)
  if isfield(subTomoMeta,'updatedWeights')
    if ~(subTomoMeta.updatedWeights)
      emc.track_stats = false;
    end
  else
    emc.track_stats = false;
  end
end

%%% For general release, I've disabled class average alignment and
%%% multi-reference alignment, so set the default to OFF. If either of
%%% these features are re-introduced, this will need to be reverted.
% FIXME: get rid of the -1

flgGold=1;

% Optionally specify gpu idxs
if (numel(emc.nGPUs) == 1)
  gpuList = 1:emc.nGPUs;
else
  gpuList = emc.nGPUs;
  emc.nGPUs = length(gpuList);
end

if (emc.nPeaks > 1)
  fprintf('For ML approach:\nUsing a compression factor %3.3f\nUsing an angulare tolerance of %3.3f degrees\n', ...
    emc.ML_compressByFactor, emc.ML_angleTolerance);
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

% FIXME: It doesn' tmake any sense to average the classes then also average those together.
% This is a waste of memory as we could just create a single average. I think this must have been
% a hack override for multi-ref alignment with classification or something.
saveClassSum = -1;

% The prefix of the variable read in from the parameter file.
parameterPrefix = '';
% The prefix used in file names and subTomoMeta fields
savePrefix = '';
switch STAGEofALIGNMENT
  case 'RawAlignment'

    fieldPrefix = 'Ref';
    savePrefix = 'Ref';
    
    % FIXME: put these into parser, change to deprecated and use Ref
    classVector{1}  = emc.(sprintf('%s_classes_odd','Raw'));
    classVector{2}  = emc.(sprintf('%s_classes_eve','Raw'));
    
    className    = emc.(sprintf('%s_className','Raw'));
    samplingRate = emc.('Ali_samplingRate');
    
    if (emc.multi_reference_alignment && (test_multi_ref_diffmap || ~emc.classification))
      fprintf('\n\nMutliRef enabled.\n');
      pause(2);
      className    = emc.(sprintf('Raw_className'));
    else
      if (emc.multi_reference_alignment && emc.classification)
        fprintf('\n\nMutliRef and Classify enabled.\n');
        className = 0;
        saveClassSum = 0;
        classVector{1} = [0;1];
        classVector{2} = [0;1];
      end
    end
    
    
  case 'FinalAlignment'
    % Special case for the final cycle.
    % Assuming RawAlignment already run for this cycle and FSC is calculated
    % Goal is to re-extract odd-half, applying the xform found in fscGold
    fieldPrefix = 'Ref';
    savePrefix = 'Ref';
    classVector{1}  = emc.(sprintf('%s_classes_odd',fieldPrefix));
    classVector{2}  = emc.(sprintf('%s_classes_eve',fieldPrefix));
    
    className    = emc.(sprintf('%s_className',fieldPrefix));
    samplingRate = emc.('Ali_samplingRate');
    
    flgFinalAvg = 1;
    % Update at some point to handle multiple classes, for now just test on the
    % "global" or whatever requested class
    iRefPrev = 1;
    
    aliParams = subTomoMeta.(cycleNumber).('fitFSC').(sprintf('Resample%s%d','Ref',iRefPrev));
    oddRot = reshape(aliParams(1,:),3,3);
    % refine the translation per particle.

    clear iRefPrev
    
    
    
  case 'Cluster_cls'
    STAGEofALIGNMENT = 'Cluster';
    ClusterGeomNAME = 'ClusterClsGeom';
    fieldPrefix = 'Cls';
    savePrefix = 'Cls';
    
    classVector{1}  = emc.(sprintf('%s_classes_odd',fieldPrefix));
    classVector{2}  = emc.(sprintf('%s_classes_eve',fieldPrefix));
    
    classCoeffs{1} =  emc.('Pca_coeffs');
    classCoeffs{2} =  emc.('Pca_coeffs');
    
    samplingRate = emc.(sprintf('Cls_samplingRate'));
    className    = emc.(sprintf('%s_className',fieldPrefix));
    if (emc.classification)
      flgGold = 0;
    end
    
  case 'SnrEstimate'
    
    fieldPrefix = 'Ref';
    savePrefix = 'Ref';

    classVector{1}  = [1:25;ones(1,25)];
    classVector{2}  = [1:25;ones(1,25)];
    
    className    = 25;
    samplingRate = emc.(sprintf('%s_samplingRate','Ali'));
    
  otherwise
    error('STAGEofALIGNMENT incorrect');
end


outputPrefix = sprintf('%s_%s',cycleNumber, emc.('subTomoMeta'));

emc.pixel_size_angstroms = emc.pixel_size_angstroms .* samplingRate;
peakSearch   = floor(0.85.*emc.('particleRadius')./emc.pixel_size_angstroms);
peakCOM      = [1,1,1].*3;




if ~(isfield(subTomoMeta,'currentCycle'))
  subTomoMeta.('currentCycle') = 0;
end

% Somewhere I am saving currentCycle as a string. Haven't taken the time to
% track it down, but probably in this function.
if isa(subTomoMeta.('currentCycle'), 'char')
  subTomoMeta.('currentCycle') = EMC_str2double(subTomoMeta.('currentCycle'))
end

if subTomoMeta.('currentCycle') == CYCLE - 1
  
  cycleRead = sprintf('cycle%0.3u', CYCLE - 1);
  % Save a backup of the cycles total geometry
  save(sprintf('%s_%s_backup.mat',cycleRead,emc.('subTomoMeta')), ...
    'subTomoMeta');
end

if strcmpi(STAGEofALIGNMENT, 'RawAlignment')
  if ( CYCLE )
    cycleRead = sprintf('cycle%0.3u', CYCLE - 1);
  else
    emc.eucentric_fit = false; % No possible updates on cycle 0
    cycleRead = sprintf('cycle%0.3u', CYCLE);
  end
else
  emc.eucentric_fit = false; % No possible updates for other stages of alignments
  cycleRead = sprintf('cycle%0.3u', CYCLE );
end

% leave averages at size appropriate for interpolation when extracting to use
% in class avg alignment.
doNotTrim = false;
eachTomo  = false;
flgEstSNR = 0;
switch STAGEofALIGNMENT
  case 'RawAlignment'
    if ( CYCLE )
      geometry = subTomoMeta.(cycleRead).RawAlign;
    else
      geometry = subTomoMeta.(cycleRead).geometry;
      eachTomo = false;%true;
    end
    if ~(emc.classification)
      doNotTrim = true;
    end
    
  case 'FinalAlignment'
    geometry = subTomoMeta.(cycleRead).Avg_geometry;
    if ~(emc.classification)
      doNotTrim = true;
    end
    
  case 'Cluster'
    cN = cell(2,1);
    
    if (flgGold)
      cN{1} = sprintf('%s_%d_%d_nClass_%d_ODD',outputPrefix,classCoeffs{1}(1,1), classCoeffs{1}(1,end), className);
      
      geometry{1} = subTomoMeta.(cycleRead).ClusterResults.(cN{1});
      cN{2} = sprintf('%s_%d_%d_nClass_%d_EVE',outputPrefix,classCoeffs{2}(1,1), classCoeffs{2}(1,end), className);
      geometry{2} = subTomoMeta.(cycleRead).ClusterResults.(cN{2});
      
      geometry = BH_mergeClassGeometry(geometry{1}, geometry{2});
    else
      if (test_fuzz)
        cN{1} = sprintf('%s_%d_%d_nClass_%d_STD','cycle002_full_2',classCoeffs{1}(1,1), classCoeffs{1}(1,end), className);
      else
        cN{1} = sprintf('%s_%d_%d_nClass_%d_STD',outputPrefix,classCoeffs{1}(1,1), classCoeffs{1}(1,end), className);
      end
      
      geometry = subTomoMeta.(cycleRead).ClusterResults.(cN{1});
    end
    
    
    subTomoMeta.(cycleRead).('KmsSampling') = samplingRate;
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





class_idx='';
class_weights='';
if (test_fuzz)
  % FIXME: testing fuzzy classification, this is hardcoded
  class_idx = subTomoMeta.(cycleRead).ClusterResults.cycle002_full_2_0_64_nClass_36_STD_idxList;
  class_weights = subTomoMeta.(cycleRead).ClusterResults.cycle002_full_2_0_64_nClass_36_STD_p;
end

if isfield(subTomoMeta,('tomoCPR_run_in_cycle'))
  
  if (emc.eucentric_fit && ~isfield(subTomoMeta.(sprintf('%s',cycleRead)), 'eucentric_shifts'))
    cycle_to_update = subTomoMeta.('tomoCPR_run_in_cycle')(find(subTomoMeta.('tomoCPR_run_in_cycle')(:,1) == subTomoMeta.currentTomoCPR),2);
    if (cycle_to_update == cycleRead)
      error('You specified eucentric_fit=1, and you are averaging cycle %d and no shifts are found from cycle %d\n',cycleNumber,cycleRead);
    else
      emc.eucentric_fit = false;
    end
  end
else
  emc.eucentric_fit = false;
end
% Get the number of tomograms to process.
tomoList = fieldnames(geometry);
nTomograms = length(tomoList);


if (emc.classification)
  
  [ ~, maskSize, maskRadius, maskCenter ] = ...
    BH_multi_maskCheck(emc, 'Ali', emc.pixel_size_angstroms);
  % These are used when 'Cluster' is called, to take the masking parameters
  % from focused PCA/Classification, to produce a montage with reduced
  % Z-dimension & low pass filtering to be used in decision making but not
  % alignment. The X,Y dimensions must not change compared to the full
  % version or else graphical deletion of classes will fail.
  
  
  [~, ~, ~, pcaMaskCenter ] = ...
    BH_multi_maskCheck(emc, 'Cls', emc.pixel_size_angstroms);
else
  
  [ ~, maskSize, maskRadius, maskCenter ] = ...
    BH_multi_maskCheck(emc, 'Ali', emc.pixel_size_angstroms);
  
end

[ sizeWindow, sizeCalc, sizeMask, padWindow, padCalc] = ...
  BH_multi_validArea( maskSize, maskRadius, emc.scale_calc_size );

[ nParProcesses, iterList] = BH_multi_parallelJobs(nTomograms, emc.nGPUs, sizeCalc(1), emc.nCpuCores);

if any( (sizeCalc - sizeWindow) < 0 )
  sizeCalc = BH_multi_iterator( sizeWindow, 'fourier' )
end
sizeMask = sizeWindow;

padCalc = BH_multi_padVal(sizeWindow, sizeCalc);

if (flgFinalAvg)
  % Load in the even (low-passed) ref and use this only for a (non-iterative)
  % translational alignment.
  
  % Assuming the mask/window sizes haven't changed since the original average
  % that was used to calc the FSC -- that would break things here.
  imgNAME = sprintf('class_%d_Locations_REF_%s', className, 'EVE');
  weightNAME = sprintf('class_%d_Locations_REF_%s_Wgt', className, 'EVE');
  
  [ refIMG ] = BH_unStackMontage4d(1, ...
    subTomoMeta.(cycleNumber).(imgNAME){1}, ...
    subTomoMeta.(cycleNumber).(imgNAME){2},...
    sizeWindow);
  
  [ refWDG ] = BH_unStackMontage4d(1, ...
    subTomoMeta.(cycleNumber).(weightNAME){1},...
    subTomoMeta.(cycleNumber).(weightNAME){2},...
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



try
  EMC_parpool(nParProcesses+1)
catch
  delete(gcp('nocreate'))
  EMC_parpool(nParProcesses+1)
end

maxCCC = 0;

spike_info = struct();
spike_info.('std_dev') = nan;
if (emc.flgQualityWeight)
  %get the average CCC for calculation of particle quality weighting.
  
  cccVect = [];
  wgtVect = [];
  angVect = [];
  chiVect = [];
  
  if (emc.spike_prior)
    
    tiltList_tmp = fieldnames(subTomoMeta.mapBackGeometry);
    tiltList_tmp = tiltList_tmp(~ismember(tiltList_tmp,{'viewGroups','tomoName'}));
    nST = 1; tiltList = {};
    % First make sure this tilt actualy has tomos. Why is this here/
    for iStack = 1:length(tiltList_tmp)
      if subTomoMeta.mapBackGeometry.(tiltList_tmp{iStack}).nTomos
        tiltList{nST} = tiltList_tmp{iStack};
        nST = nST +1;
      end
    end
    
    
    % Now loop over all of the tomograms. In the first loop get the
    % distribution characterizing the sphericity of the data (if that's a
    % word?) i.e. make sure no principle axes are way to big, due to
    % points from adjacent virions that were not removed in
    % cleanTemplateSearch.
    f = fieldnames(subTomoMeta.mapBackGeometry.tomoName);
    
    for iTomo = 1:length(f)
      
      tmpTomo = [];
      spike_info.(f{iTomo}).('angular_diff') = zeros(size(geometry.(f{iTomo}) , 1),emc.nPeaks,'single');
      spike_info.(f{iTomo}).('normal_distance') = zeros(size(geometry.(f{iTomo}) , 1),emc.nPeaks,'single');
      
      spike_info.(f{iTomo}).('angular_prob') = zeros(size(geometry.(f{iTomo}) , 1),emc.nPeaks,'single');
      spike_info.(f{iTomo}).('angular_weight') = zeros(size(geometry.(f{iTomo}) , 1),emc.nPeaks,'single');
      
      nSubTomos = size(geometry.(f{iTomo}) , 1);
      particle_coords = zeros(nSubTomos .* emc.nPeaks,8,'single');
      nVol = 1;
      for iSubTomo = 1:nSubTomos
        for iPeak = 1:emc.nPeaks
          particle_coords(nVol,1:5) = geometry.(f{iTomo})(iSubTomo,[26,4,11:13]+(iPeak-1)*26);
          nVol = nVol + 1;
        end
      end
      % Logical size nsubtomos x emc.nPeaks
      positions_to_analyze = particle_coords(:,1) ~= -9999;
      display_fit = false;
      radial_shrink_factor = 2;
      
      
      [ normal_vect, chi2 ] = BH_fit_ellipsoidal_prior(emc.pixel_size_angstroms .* particle_coords(positions_to_analyze,3:5), ...
        emc.('particleRadius')(3), ...
        radial_shrink_factor, ...
        display_fit);
      
      chiVect = [chiVect chi2];
      particle_coords(positions_to_analyze,[6:8]) = [ normal_vect];
      nVol = 1;
      for iSubTomo = 1:nSubTomos
        for iPeak = 1:emc.nPeaks
          if (particle_coords(nVol,1) ~= -9999)
            particleAxis = reshape(geometry.(f{iTomo})(iSubTomo,[17:25]+(iPeak-1)*26),3,3)*[0;0;1];
            angularDiff = dot(particle_coords(nVol,6:8), particleAxis);
            if (abs(angularDiff) > 1)
              angularDiff = fix(angularDiff);
            end
            angularDiff = acosd(angularDiff);
            spike_info.(f{iTomo}).('angular_diff')(iSubTomo,iPeak) = angularDiff;
            tmpTomo = [tmpTomo, angularDiff];
            
          end
          nVol = nVol + 1;
        end % iPeak
      end % iSubTomo
      
      angVect = [angVect tmpTomo];
      
    end % iTomo
    
    
    spike_info.('std_dev') = std(angVect);
    h = histogram(angVect,'Normalization','probability','BinMethod','fd');
    hv = h.Values;
    %       [~,mc] = max(hv);
    %       hv(1:mc-1) = hv(mc);
    hv = (hv ./ max(hv(:))) .^ 0.5; %(mean(chiVect)./std(chiVect).^2);
    x = h.BinWidth:h.BinWidth:h.BinLimits(2);
    v = 0:h.BinLimits(2)./1000:h.BinLimits(2);
    spike_info.('angular_pdf') = griddedInterpolant(x,hv,'makima','linear');
    figure('Visible','off'), bar(x,h.Values,'w'); hold on
    plot(v,spike_info.('angular_pdf')(v),'b','linewidth',2)
    sprintf(' %2.3f degrees',spike_info.('std_dev'))
    
    title(sprintf('Spike Angle Prior, std-dev %2.3f degrees',spike_info.('std_dev')));
    file_out = sprintf('%s-spike-angle-prior.pdf', cycleNumber);
    saveas(gcf, file_out,'pdf')
    hold off;
    
    
    save('spike_hist.mat','angVect','chiVect');
    %       figure,
  end % if spike_prior
  
  nVolumes = 0;
  addedWeight = 0;
  for iParProc = 1:nParProcesses
    for iTomo = iterList{iParProc}
      if (emc.track_stats)
        geometry.(tomoList{iTomo})(:,1:26:26*emc.nPeaks) = geometry.(tomoList{iTomo})(:,1:26:26*emc.nPeaks)./geometry.(tomoList{iTomo})(:,2:26:26*emc.nPeaks);
      end
      
      min_weight = 1e-6;
      if (emc.spike_prior)
        for iSubTomo = 1:size(geometry.(tomoList{iTomo}) , 1)
          peakList = false(emc.nPeaks,1);
          for iPeak = 1:emc.nPeaks
            if (geometry.(tomoList{iTomo})(iSubTomo,26*iPeak)~=-9999)
              iWeight =    ...
                spike_info.('angular_pdf')(spike_info.(tomoList{iTomo}).('angular_diff')(iSubTomo,iPeak));
              peakList(iPeak) = true;
              if (iWeight < min_weight || ~isfinite(iWeight))
                iWeight = min_weight;
              end
              
              spike_info.(tomoList{iTomo}).('angular_prob')(iSubTomo,iPeak) = iWeight;
              
              nVolumes = nVolumes + 1;
              addedWeight = addedWeight + iWeight;
            end
          end
          spike_info.(tomoList{iTomo}).('angular_prob')(iSubTomo,peakList) = ...
            spike_info.(tomoList{iTomo}).('angular_prob')(iSubTomo,peakList) ./ ...
            sum(spike_info.(tomoList{iTomo}).('angular_prob')(iSubTomo,peakList));
          for iScoreMod = 1:emc.nPeaks
            if (peakList(iScoreMod))
              geometry.(tomoList{iTomo})(iSubTomo,2 + 26*(iScoreMod-1)) = ...
                spike_info.(tomoList{iTomo}).('angular_prob')(iSubTomo,iScoreMod);
            end
          end
          
          
          
        end
        
      end
      
      keepVect = geometry.(tomoList{iTomo})(:,26:26:26*emc.nPeaks)~=-9999 ;
      
      tmpVect = geometry.(tomoList{iTomo})(:,1:26:26*emc.nPeaks);
      
      
      
      cccVect = [cccVect ; reshape(tmpVect(keepVect),[],1)];
      tmpVect = geometry.(tomoList{iTomo})(:,2:26:26*emc.nPeaks);
      
      wgtVect = [wgtVect ; reshape(tmpVect(keepVect),[],1)];
    end
  end
  
  
  if (emc.ccc_cutoff > 1.0)
    sorted_ccc = sort(cccVect);
    reqVol = int32(round(emc.ccc_cutoff))
    length(sorted_ccc) - reqVol
    emc.ccc_cutoff = sorted_ccc(length(sorted_ccc) - reqVol);
    fprintf('Removing all volumes with score < %2.2f to return the requested %d volumes\n\n',emc.ccc_cutoff,reqVol);
  elseif (emc.ccc_cutoff > 0.0)
    sorted_ccc = sort(cccVect);
    reqVol = emc.ccc_cutoff;
    emc.ccc_cutoff = sorted_ccc(floor(length(cccVect).*(1 - reqVol)));
    fprintf('Removing all volumes with score < %2.2f to return the requested percent %2.2f of possible volumes\n\n',emc.ccc_cutoff,reqVol);
  end
  
  subTomoMeta.(cycleNumber).('score_sigma') = std(cccVect);
  if (emc.spike_prior)
    %       spike_info.('normalization_factor') = 1;%nVolumes ./ (emc.nPeaks * addedWeight);
    %       fprintf('From %d possible volumes the total weight is %3.3e\n',nVolumes,addedWeight);
  end
  avgCCC = mean(cccVect);
  
  
  
  mean(wgtVect)
  std(wgtVect)
  maxCCC = max(cccVect);
  mean(wgtVect)
  median(wgtVect)
  
  %     figure, hist(cccVect,29)
  %     figure, hist(wgtVect,29)
  %     figure, hist((wgtVect./median(wgtVect)).^weightScale,29)
  %     error('asdf')
  
  if (emc.track_stats)
    fprintf('Avgerage score is %3.3f, using a quality weight of %2.2f\n\n',avgCCC,emc.flgQualityWeight);
  else
    fprintf('Avgerage CCC is %3.3f, using a quality weight of %2.2f\n\n',avgCCC,emc.flgQualityWeight);
  end
  
  
else
  maxCCC = [];
  avgCCC = [];
  subTomoMeta.(cycleNumber).('score_sigma') = 1;
end

% % Clear all of the GPUs prior to entering the main processing loop
for iGPU = 1:emc.nGPUs
  gpuDevice(iGPU);
end

parVect = 1:nParProcesses;
parfor iParProc = parVect
  % for iParProc = parVect
  
  % Get the gpuIDX assigned to this process
  gpuIDXList = mod(parVect+emc.nGPUs,emc.nGPUs)+1;
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
  
  
  
  if (emc.flgQualityWeight)
    [cccWeight,~,~,~,~,~] = BH_multi_gridCoordinates(sizeCalc, ...
      'Cartesian','GPU',...
      {'none'},1,0,1);
    cccWeight = (cccWeight ./ emc.pixel_size_angstroms).^2;
    
    
    
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
    
    if (emc.eucentric_fit)
      try
        geometry_tmp.(tomoList{iTomo})(:,13) = geometry_tmp.(tomoList{iTomo})(:,13) + ...
          subTomoMeta.(sprintf('%s',cycleRead)).('eucentric_shifts').(tomoList{iTomo}) ;
      catch
        fprintf('WARNING, did not find the eucentric shift for tomo %s\n', tomoList{iTomo});
      end
    end
    
    peakMask_tmp = gpuArray(peakMask);
    peakBinary_tmp = gpuArray(peakBinary);
    

    interpMask_tmpBinary = gpuArray(single(find(interpMask > 0.01 )));
    interpMask_tmp = gpuArray(interpMask);
    interpMaskWdg_tmp = (gpuArray(interpMaskWdg));
    
    if (eachTomo)
      tomoAvg = zeros(sizeMask, 'single', 'gpuArray');
      tomoWgt = zeros(sizeCalc , 'single', 'gpuArray');
      tomoCount = 0;
    end
    
    fprintf('gpu %d working on %d/%d volumes\n',iParProc,iTomo,nTomograms);
    
    tomoName = tomoList{iTomo};
    
    
    
    
    tiltGeometry = subTomoMeta.tiltGeometry.(tomoList{iTomo});
    tomoNumber = subTomoMeta.mapBackGeometry.tomoName.(tomoList{iTomo}).tomoNumber;
    tiltName   = subTomoMeta.mapBackGeometry.tomoName.(tomoList{iTomo}).tiltName;
    coords = subTomoMeta.mapBackGeometry.(tiltName).coords(tomoNumber,1:4);
    [ binShift ] = [0,0,0];%BH_multi_calcBinShift( coords, samplingRate);
    
    % Load in the geometry for the tomogram, and get number of subTomos.
    positionList = geometry_tmp.(tomoList{iTomo});
    nSubTomos = sum(any(positionList(:,26:26:26*emc.nPeaks) ~= -9999,2));
    
    nSubTomosTotal = nSubTomosTotal + nSubTomos;
    
    volumeData = [];
    %fprintf('loading tomo %d\n',iTomo);
    
    tomoNumber = subTomoMeta.mapBackGeometry.tomoName.(tomoList{iTomo}).tomoNumber;
    tiltName = subTomoMeta.mapBackGeometry.tomoName.(tomoList{iTomo}).tiltName;
    reconCoords = subTomoMeta.mapBackGeometry.(tiltName).coords(tomoNumber,:);
    reconGeometry = (subTomoMeta.reconGeometry.(tomoList{iTomo}) ./ samplingRate);

    if (emc.flgCutOutVolumes && ~volumesNeedToBeExtracted)
      volumeData = [];
    else

      reconScaling = 1;
      [ volumeData, ~ ] = BH_multi_loadOrBuild( ...
        tomoList{iTomo}, ...
        reconCoords, mapBackIter, ...
        samplingRate,iGPUidx,reconScaling,0);
      
      volHeader = getHeader(volumeData);
    end
    
    
    
    iTiltName = subTomoMeta.mapBackGeometry.tomoName.(tomoList{iTomo}).tiltName;
    wgtName = sprintf('cache/%s_bin%d.wgt',iTiltName,samplingRate);
    
    
    % Work on each class seperately pushing to main memory when finished.
    for iGold = 1:2-flgFinalAvg
      for iClassPos = 1:size(classVector{iGold},2)
        
        
        iTempParticleODD = zeros(sizeMask,'single','gpuArray');
        iTempWedgeODD = zeros(sizeCalc,'single','gpuArray');
        
        
        iTempParticleEVE = zeros(sizeMask,'single','gpuArray');
        iTempWedgeEVE = zeros(sizeCalc,'single','gpuArray');
        
        
        iClassIDX = classVector{iGold}(1,iClassPos);
        %     nClassPossible = nClassPossible + ...
        %                                   length(find(positionList(:,26) == iClassIDX));
        if classVector{iGold}(1,iClassPos)
          % pick just the particular class to work with
          
          if ( flgEstSNR )
            % When the class is for estimating SNR
            includeList = ( any(abs(positionList(:,1:26:26*emc.nPeaks))  >= emc.ccc_cutoff,2) & ...
              positionList(:,10) ==  iClassIDX & ...
              positionList(:,7)  == iGold );
          else
            % When the class is from statistical analysis
            includeList = ( any(abs(positionList(:,1:26:26*emc.nPeaks)) >= emc.ccc_cutoff,2)  & ...
              positionList(:,26) ==  iClassIDX & ...
              positionList(:,7)  == iGold );
          end
          
        else
          % if class is 0, pick all non-ignored particles
          
          includeList = ( any(abs(positionList(:,1:26:26*emc.nPeaks))  >= emc.ccc_cutoff,2)  & ...
            any(positionList(:,26:26:26*emc.nPeaks) ~= -9999,2)  & ...
            positionList(:,7)  == iGold );
          
        end
        % Get the position index for each included particle.
        
        particleIndex = find(includeList);
        nClassPossible = length(particleIndex);
        
        
        
        for iSubTomo = particleIndex'
          iParticle = [];
          iCCCweight = [];
          iWedgeMask = [];
          
          if ( emc.nPeaks > 1 )
            % Calculate a relative weighting, normalize max score to one
            % and then raise to compressBy factor to downweight lower
            % scores.
            
            [ peakWgt, sortedList ] = BH_weightAngCheckPeaks( ...
              positionList(iSubTomo,:),...
              emc.nPeaks,  ...
              subTomoMeta.(cycleNumber).('score_sigma') ,...
              iSubTomo, tomoList{iTomo},...
              emc.track_stats);
            % Update any re-ordering or elimination
            positionList(iSubTomo,:) = sortedList;
          else
            peakWgt = 1;
          end
          
          % %           if (emc.spike_prior)
          % %             peakWgt = peakWgt.*spike_info.(tomoList{iTomo}).('angular_prob')(iSubTomo,1).*spike_info.('normalization_factor');
          % %           end
          
          make_sf3d = true;
          
          for iPeak = 1:emc.nPeaks
            
            if peakWgt(iPeak) == -9999
              positionList(iSubTomo, 26*iPeak) = -9999;
              % Skip this peak
              continue
            end
            %Check that the given subTomo is not to be ignored
            
            
            % Get position and rotation info, angles stored as e1,e3,e2 as in AV3
            % and PEET. This also makes inplane shifts easier to see.
            
            center = positionList(iSubTomo,[11:13]+26*(iPeak-1))./samplingRate + binShift;
            angles = positionList(iSubTomo,[17:25]+26*(iPeak-1));
            wdgIDX = positionList(iSubTomo,9);
            
            % tmpang = BH_defineMatrix([0,0,-14],'Bah','inv');
            % angles = reshape(angles,3,3)*tmpang;
            % if (flgFinalAvg)
            %   angles = reshape(angles,3,3)*oddRot;
            % end
            
            TLT = subTomoMeta.('tiltGeometry').(tomoList{iTomo});
            
            if (make_sf3d)
              [ iSF3D ] = BH_weightMaskMex(sizeCalc, samplingRate, TLT, center,reconGeometry, emc.wiener_constant);
              make_sf3d = false;
            end
            
            
            
            if (emc.flgQualityWeight)
              iCCC = positionList(iSubTomo,[1]+26*(iPeak-1));
              
              if (emc.track_stats)
                % Downweight higher frequency in all subTomos with iCCC below the mean
                iBfactor = (emc.flgQualityWeight.*(iCCC - maxCCC)./4)
                iCCCweight = exp(iBfactor.*cccWeight);
              else
                if iCCC < avgCCC
                  % Downweight higher frequency in all subTomos with iCCC below the mean
                  iBfactor = -1.*(emc.flgQualityWeight.*(acosd(iCCC) - acosd(avgCCC)))^2;
                  iCCCweight = exp(iBfactor.*cccWeight);
                  
                else
                  iCCCweight=1;
                end
                
              end
              
              if ( any(emc.filterDefocus))
                iDef = abs(mean(tiltGeometry(:,15))*10^6);
                iDef = (emc.filterDefocus(1)*max(iDef-1,0.5))^emc.filterDefocus(2);
                fprintf('Using iDef %f\n',iDef);
                % Frequency is already squared so adjust to match iDef scale
                % factor.
                iCCCweight = iCCCweight.*exp(iDef.*cccWeight.^(emc.filterDefocus(2)/2));
              end
            end
            
            % Find range to extract, and check for domain error.
            if (emc.flgCutOutVolumes && ~volumesNeedToBeExtracted)
              [ indVAL, padVAL, shiftVAL ] = ...
                BH_isWindowValid(2*CUTPADDING+sizeWindow, ...
                sizeWindow, maskRadius, center);
            else
              [ indVAL, padVAL, shiftVAL ] = ...
                BH_isWindowValid([volHeader.nX,volHeader.nY,volHeader.nZ], ...
                sizeWindow, maskRadius, center);
            end
            
           
            if ~ischar(indVAL)
 
              if (emc.flgCutOutVolumes && ~volumesNeedToBeExtracted)
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

                iParticle = gpuArray(getVolume(volumeData,[indVAL(1,1),indVAL(2,1)], ...
                  [indVAL(1,2),indVAL(2,2)], ...
                  [indVAL(1,3),indVAL(2,3)],'keep'));
              end
              
              
              
              [ iParticle ] = BH_padZeros3d(iParticle,  padVAL(1,1:3), ...
                padVAL(2,1:3), 'GPU', 'single');
              
              
              
              if (emc.flgCutOutVolumes && volumesNeedToBeExtracted)
                % Test with some generic padding , only to be used on bin 1 at
                % first!!! TODO add a flag to check this.
                
                
                
                particleOUT_name = sprintf('cache/subtomo_%0.7d_%d.mrc',positionList(iSubTomo,4),iPeak);
                positionList(iSubTomo,[11:13]+26*(iPeak-1)) = shiftVAL+CUTPADDING+ceil((sizeWindow+1)./2);
                if (emc.projectVolumes)
                  SAVE_IMG(sum(iParticle,3),particleOUT_name,emc.pixel_size_angstroms);
                else
                  particleOUT = BH_padZeros3d(gather(iParticle), CUTPADDING.*[1,1,1], ...
                    CUTPADDING.*[1,1,1], 'cpu', 'single');
                  SAVE_IMG(particleOUT,particleOUT_name,emc.pixel_size_angstroms);
                end
                
                
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
                [~, iRefWdg] = interpolator(gpuArray(refWDG{1}),angles',[0,0,0], rotConvention , 'forward', 'C1', true);
                [~, iRefIMG] = interpolator(gpuArray(refIMG{1}),angles',iShift, rotConvention , 'forward', 'C1', true);
                %             iRefIMG = BH_resample3d(refIMG{1},angles', iShift,rotConvention , 'GPU', 'forward');
                %             iRefWdg = BH_resample3d(refWDG{1},angles', [0,0,0],rotConvention , 'GPU', 'forward');
                
                [ ref_FT ] = BH_bandLimitCenterNormalize(iRefIMG.*peakMask_tmp, ...
                  ifftshift(iSF3D), ...
                  peakBinary_tmp,padCalc,...
                  'single');
                
                [ part_FT ] = BH_bandLimitCenterNormalize(iParticle.*peakMask_tmp, ...
                  ifftshift(iRefWdg), ...
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
              
              
              %iParticle is already on GPU if it should be.
              [~, iParticle] = interpolator(gpuArray(iParticle), angles, iShift, rotConvention , 'inv', emc.symmetry, true);
              
              
              [~, iWedgeMask] = interpolator(gpuArray(iSF3D), angles, [0,0,0], rotConvention , 'inv', emc.symmetry, true);
                            
              
              iParticle = iParticle -  mean(iParticle(interpMask_tmpBinary));
              iParticle = iParticle ./  rms(iParticle(interpMask_tmpBinary));
              iParticle = iParticle .* interpMask_tmp;
              if (emc.flgQualityWeight && numel(iCCCweight) > 1)
                
                iParticle = real(ifftn(fftn(BH_padZeros3d(iParticle,...
                  'fwd',padCalc,'GPU','singleTaper')).*iCCCweight));
                
                iParticle = BH_padZeros3d(iParticle,...
                  'inv',padCalc,'GPU','singleTaper');
                iWedgeMask = iWedgeMask .* fftshift(iCCCweight);
              end
              
              
              
              
              trimAvg =  mean(iParticle(:));
              
              if ~isfinite(trimAvg)
                fprintf('SubTomo %d from tomogram %s as Nan mean\n',...
                  iSubTomo, tomoList{iTomo});
                iParticle(:,:,:) = 0;
                % Flag the particle as ignored
                positionList(iSubTomo, 26:26:emc.nPeaks*26) = -9999;
              else
                
                if (test_fuzz)
                  % Find the right row in the weight array
                  iProb_row = find(class_idx == positionList(iSubTomo,4));
                  iProb = class_weights(iProb_row,:);
                  % iProb_weight = 4;
                  % iProb = iProb.^iProb_weight ./ sum(iProb.^iProb_weight);
                  iParticle = gather(iParticle .* peakWgt(iPeak));
                  iWedgeMask = gather(iWedgeMask .* peakWgt(iPeak));
                  
                  % Results are weird so override and see if it replicates the "normal" behavior by only taking the class with the minimum
                  % distance
                  [tmin,tidx] = min(iProb);
                  iProb = iProb.*0;
                  iProb(tidx) = 1;
                  
                  for iWeight = 1:length(iProb)
                    avgVolume_tmp{iWeight, positionList(iSubTomo,7)} =  ...
                      avgVolume_tmp{iWeight, positionList(iSubTomo,7)}  + (iParticle .* iProb(iWeight));
                    
                    avgWedge_tmp{ iWeight, positionList(iSubTomo,7)} = ...
                      avgWedge_tmp{ iWeight, positionList(iSubTomo,7)}  + (iWedgeMask .* iProb(iWeight));
                    nExtracted_tmp(iWeight, positionList(iSubTomo,7)) = ...
                      nExtracted_tmp(iWeight, positionList(iSubTomo,7)) + iProb(iWeight);
                  end
                  
                else
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
                  
                end % end test fuzz
                
              end
              
            else
              
              fprintf('SubTomo %d from tomogram %s only sampled at %f\n',...
                iSubTomo, tomoList{iTomo}, 1-padVAL);
              
              % Flag the particle as ignored
              positionList(iSubTomo, 26:26:26*emc.nPeaks) = -9999;
              nIgnored = nIgnored + 1;
              peakWgt(1:emc.nPeaks) = -9999;
              
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
          % for the fuzz each particle goes into every class, so we need to pull earlier
          if ~(test_fuzz)
            avgVolume_tmp{ iClassPos, 1} =  avgVolume_tmp{ iClassPos, 1} + ...
              gather(iTempParticleODD);
            avgVolume_tmp{ iClassPos, 2} =  avgVolume_tmp{ iClassPos, 2} + ...
              gather(iTempParticleEVE);
            avgWedge_tmp{  iClassPos, 1} = avgWedge_tmp{  iClassPos, 1} + ...
              gather(iTempWedgeODD);
            avgWedge_tmp{  iClassPos, 2} = avgWedge_tmp{  iClassPos, 2} + ...
              gather(iTempWedgeEVE);
          end
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
    for iHalf = 1:2-flgFinalAvg
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




subTomoMeta.(cycleNumber).('nSubTomoAveraged') = gather(sum(nExtracted(:)));
subTomoMeta.(cycleNumber).(sprintf('newIgnored_Avg%s',fieldPrefix)) = ...
  gather(nIgnored);


% get the total class average by combining eve/odd
classStorage = cell(maxClasses,2);
if (doNotTrim) && (emc.classification)
  filteredClass= cell(maxClasses,2);
  % low-pass to see class averages more clearly.
  % %   [ bandpassFilt ] = BH_bandpass3d( sizeMask, 0.2, 300, 30, 'GPU',emc.pixel_size_angstroms);
end

for iClass = 1:maxClasses
  classStorage{iClass,1} = zeros(sizeMask, 'single');
  classStorage{iClass,2} = zeros(sizeMask, 'single');
end

if (eachTomo)
  system('mkdir -p initialTomoAvgs');
  % sizeWeight mask is sizeMask or 128^3 whichever is larger
  %   bandpassFiltTomo = BH_bandpass3d( sizeCalc, lpTomo(1), lpTomo(2), lpTomo(3), 'GPU',emc.pixel_size_angstroms);
  
  for iParProc = 1:nParProcesses
    nTomos = 1;
    
    for iTomo = iterList{iParProc}
      tomoName = sprintf('initialTomoAvgs/%s.mrc',tomoList{iTomo});
      classTmp = avgTomoResults{iParProc}{nTomos};
      avgTomoResults{iParProc}{nTomos} = [];
      
      SAVE_IMG(classTmp,tomoName,emc.pixel_size_angstroms);
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
  
  if (doNotTrim) && (emc.classification)
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
  
  for iGold = 1:2-flgFinalAvg
    if iGold == 1
      halfSet = 'ODD';
    else
      halfSet = 'EVE';
    end   
    
    classStorage{iClassPos,iGold} = gather(avgVolume{iClassPos,iGold} );
    
    if ~isfinite( mean(classStorage{iClassPos,iGold}(:)) )
      clear classAVG
      fprintf('zeroing out classavg because of NaN values detected.\n')
    else
      
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
if  strcmpi(STAGEofALIGNMENT, 'Cluster') && ~(emc.classification)
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
  
  subTomoMeta.(cycleNumber).(sprintf('class_%d_%s_EveOddIdx',className,fieldPrefix)) = classListOut;
  
  
end

% Second option allows re-use of class designations to generate
% multi-reference alignment.
if strcmpi(STAGEofALIGNMENT, 'Cluster')
  subTomoMeta.(cycleNumber).(ClusterGeomNAME) = geometry;
elseif strcmpi(STAGEofALIGNMENT, 'RawAlignment') && emc.multi_reference_alignment
  if (emc.classification)
    subTomoMeta.(cycleNumber).('ClusterClsGeom') = geometry;
  else
    subTomoMeta.(cycleNumber).('ClusterRefGeom') = geometry;
  end
else
  subTomoMeta.(cycleNumber).('Avg_geometry') = geometry;
end
% Should this save differently depending on the stage of alignment?? I
% think so but leave alone for nw.


for iGold = 1:2-flgFinalAvg
  if iGold == 1
    halfSet = 'ODD';
  else
    halfSet = 'EVE';
  end
  imgCounts = gather([classVector{iGold}(1,:) ; nExtracted(:,iGold)']);
  
  
  [montOUT, imgLocations] = BH_montage4d(avgVolume(:,iGold), '');
  
  
  imout = sprintf('%s_class%d_%s_%s_NoWgt.mrc',outputPrefix, ...
    className, fieldPrefix, halfSet);
  classOut = sprintf('class_%d_Locations_%s_%s_NoWgt', className,fieldPrefix, halfSet);
  
  subTomoMeta.(cycleNumber).(classOut) = {imout,imgLocations,imgCounts};
  
  if (flgFinalAvg)
    system(sprintf('mv %s preHalfSetAli_%s',imout,imout));
  end
  SAVE_IMG(montOUT, imout, emc.pixel_size_angstroms);
  %%%%%%%%
  [montOUT, imgLocations] = BH_montage4d(avgWedge(:,iGold), '');
  
  imout = sprintf('%s_class%d_%s_%s_Wgt.mrc',outputPrefix, ...
    className, fieldPrefix, halfSet);
  classOut = sprintf('class_%d_Locations_%s_%s_Wgt', className,fieldPrefix, halfSet);
  
  % For the weight, instead of imgCounts save the padValues
  subTomoMeta.(cycleNumber).(classOut) = {imout,imgLocations,fscPAD};
  
  if (flgFinalAvg)
    system(sprintf('mv %s preHalfSetAli_%s',imout,imout));
  end
  SAVE_IMG(montOUT, imout);

  if (saveClassSum > -1)
    imgCounts = gather([classVector{iGold}(1,:) ; nExtracted(:,iGold)']);
    
    
    [montOUT, imgLocations] = BH_montage4d(classSum(iGold), '');
    
    imout = sprintf('%s_class%d_%s_%s_NoWgt.mrc',outputPrefix, ...
      saveClassSum, savePrefix, halfSet);
    classOut = sprintf('class_%d_Locations_%s_%s_NoWgt', saveClassSum,savePrefix, halfSet);
    SAVE_IMG(montOUT, imout,emc.pixel_size_angstroms);
    subTomoMeta.(cycleNumber).(classOut) = {imout,imgLocations,imgCounts};
    
    [montOUT, imgLocations] = BH_montage4d(classWgtSum(iGold), '');
    
    imout = sprintf('%s_class%d_%s_%s_Wgt.mrc',outputPrefix, ...
      saveClassSum, savePrefix, halfSet);
    classOut = sprintf('class_%d_Locations_%s_%s_Wgt', saveClassSum, savePrefix, halfSet);
    SAVE_IMG(montOUT, imout,emc.pixel_size_angstroms);
    subTomoMeta.(cycleNumber).(classOut) = {imout,imgLocations,imgCounts};
    
  end
end



subTomoMeta = gather(subTomoMeta);
classVector = gather(classVector);
nExtracted = gather(nExtracted);



subTomoMeta.(cycleNumber).('SymmetryApplied').(STAGEofALIGNMENT) = emc.symmetry;
subTomoMeta.(cycleNumber).('ClassVector').(STAGEofALIGNMENT) = classVector;
cycleNumber = gather(cycleNumber);
subTomoMeta.('currentCycle') = gather(CYCLE);


save(emc.subTomoMeta, 'subTomoMeta');



fprintf('Total execution time : %f seconds\n', seconds(datetime("now")-startTime));

% clean everything up, since this function is called from other functions.

delete(gcp('nocreate'));
for iGPU = 1:emc.nGPUs
  gpuDevice(gpuList(iGPU));
end


if strcmpi(STAGEofALIGNMENT, 'RawAlignment')
  BH_fscGold_class(PARAMETER_FILE, num2str(CYCLE), STAGEofALIGNMENT);
end


try
  EMC_parpool(emc.nGPUs)
catch
  delete(gcp('nocreate'));
  EMC_parpool(emc.nGPUs)
end

for iGPU = 1:emc.nGPUs
  gpuDevice(gpuList(iGPU));
end




if ~( flgEstSNR )
  load(sprintf('%s.mat', emc.('subTomoMeta')), 'subTomoMeta');
  masterTM = subTomoMeta;
  %%%%%%%%%%%%%55 Reweight now that the FSC is calculated
  
  for iGold = 1:2
    if iGold == 1
      halfSet = 'ODD';
    else
      halfSet = 'EVE';
    end
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
  if (emc.multi_reference_alignment || emc.classification )
    nClassesReWgt = maxClasses;
  else
    nClassesReWgt = 1;
  end
  
  for iRef = 1:nClassesReWgt
    fprintf('Stage of alignment %s\niRef %d\n',STAGEofALIGNMENT,iRef);
    if strcmpi(STAGEofALIGNMENT, 'RawAlignment')
      savePrefix = 'Ref'
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
    
    % FIXME: logic?
    if (flgGold || emc.classification)
      flgCombine = 0;
      flgRefCutOff = 1;
    else
      flgCombine = 1;
      flgRefCutOff = 0;
    end
    
    try
      fscParams = masterTM.(cycleNumber).('fitFSC').(sprintf('%s%d',savePrefix,iRefPrev));
      aliParams = masterTM.(cycleNumber).('fitFSC').(sprintf('Resample%s%d',savePrefix,iRefPrev));
      mskParams = masterTM.(cycleNumber).('fitFSC').(sprintf('Mask%s%d',savePrefix,iRefPrev));
    catch
      fprintf('\nReverting from %s to REf in loading fitFSC\n',savePrefix);
      fscParams = masterTM.(cycleNumber).('fitFSC').(sprintf('%s%d','Ref',iRefPrev));
      aliParams = masterTM.(cycleNumber).('fitFSC').(sprintf('Resample%s%d','Ref',iRefPrev));
      mskParams = masterTM.(cycleNumber).('fitFSC').(sprintf('Mask%s%d','Ref',iRefPrev));
    end
    iOdd = iRef;
    iEve = iRef;
    
    if (flgFinalAvg)
      % negative to combine but NOT apply the xform to the odd set
      flgCombine = -1;
    end
    
    % Only send the lowest Bfactor if not flgFinalAvg
    if (flgFinalAvg)
      bFactorSend = emc.Fsc_bfactor;
    else
      bFactorSend = emc.Fsc_bfactor(1);
    end
    
    refTMP = gather(BH_multi_cRef_Vnorm(fscParams, aliParams, mskParams,...
      {refIMG{1}{iOdd},refIMG{2}{iEve}}, ...
      {refWGT{1}{iOdd},refWGT{2}{iEve}}, ...
      flgCombine,flgRefCutOff, emc.pixel_size_angstroms, bFactorSend));
    
    if ~(flgFinalAvg)
      refIMG{1}{iOdd} = refTMP{1,1};
      if (flgGold || emc.classification)
        refIMG{2}{iEve} = refTMP{1,2};
      end
      clear refTMP
    end
    
  end
  
  
  for iGold = 1:2-flgFinalAvg
    if( flgGold || emc.classification )
      if iGold == 1
        halfSet = 'ODD';
      else
        halfSet = 'EVE';
      end
    else
      halfSet = 'STD';
    end
    imgIN = sprintf('class_%d_Locations_%s_%s_NoWgt', ...
      className, fieldPrefix, halfSet);
      imgCounts = masterTM.(cycleNumber).(imgIN){3};
    % Save the unweighted, weighted imgs, weightes, optionally filtered.
    
    if (flgFinalAvg)
      for iBfactor = 1:length(emc.Fsc_bfactor)
        imout = sprintf('%s_class%d_%s_bFact-%d.mrc',outputPrefix, ...
          className, 'final',emc.Fsc_bfactor(iBfactor));
        
        SAVE_IMG(refTMP{iBfactor}, imout, emc.pixel_size_angstroms);
      end
    else
      
      [montOUT, imgLocations] = BH_montage4d(refIMG{iGold}(:), '');
      imout = sprintf('%s_class%d_%s_%s.mrc',outputPrefix, ...
        className, fieldPrefix, halfSet);
      classOut = sprintf('class_%d_Locations_%s_%s', className,fieldPrefix, halfSet);
      
      masterTM.(cycleNumber).(classOut) = {imout,imgLocations,imgCounts};
      SAVE_IMG(montOUT, imout, emc.pixel_size_angstroms);
    end
    %%%%%%%
  end
  
  
  subTomoMeta = masterTM;
  subTomoMeta.('CUTPADDING') = emc.CUTPADDING;
  if (emc.flgCutOutVolumes && volumesNeedToBeExtracted)
    subTomoMeta.('volumesAreCutOut') = 1;
  end
  save(emc.('subTomoMeta'), 'subTomoMeta');
  
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
