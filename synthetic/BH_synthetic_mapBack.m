function [ ] = BH_synthetic_mapBack(PARAMETER_FILE, CYCLE, tiltStart)

% Map back and align using the subtomograms as fiducial markers.

% If multiple classes are being mapped back, pass in the className (refName
% really b/c these will be the only re-weighted volumes) and then each class
% will be given a unique density value in a seperate volume used to visualize
% color in Chimera.
%
% Otherwise, pass a string that points at a single volume to use.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some flags that are worth keeping as options, but not accessible
% directlyCT
% by the users (private methods-ish)

% FIXME: is this even relevant any more?
buildTomo=1;% % % % % % %

save_diagnostic_ccf=0;
% If false, this is faster, simplifies the code and permits defocus estimation
% This will likely be removed in favor of deleting all the blockes under its control
use_background_estimate = false;
% Default true, we don't need this after projection
delete_background_estimate = true;



emc = BH_parseParameterFile(PARAMETER_FILE);

CYCLE = EMC_str2double(CYCLE);
cycle_numerator = '';
cycle_denominator ='';

% When tomoCPR is run on one node (start to finish) or is being finalized on a many node run using the [cycle, nodeIDX, totalNodes] syntax with we [ cycle, 0, 0 ]
% this is set to true. For the many node alignment, we want to defer running the final alignment until all the nodes have finished their work, which we do by setting this to false.
flgRunAlignments = true;
multi_node_run = false;
skip_to_the_end_and_run = false;

if numel(CYCLE) == 3
  multi_node_run = true;
  % After splitting, run the alignments while skipping everything else
  if CYCLE(2) == 0 && CYCLE(3) == 0
    skip_to_the_end_and_run = true;
  else
    flgRunAlignments = false;
  end
  cycle_numerator = CYCLE(2);
  cycle_denominator = CYCLE(3);
  CYCLE = CYCLE(1);
end

EMC_assert_numeric(CYCLE, 1, [0, inf]);

% skip_to_the_end_and_run is only relevant when running on multiple nodes
if (skip_to_the_end_and_run && ~multi_node_run)
  error('You are trying to skip to the end and run, but you are not running on multiple nodes');
end


cycleNumber = sprintf('cycle%0.3u', CYCLE);

samplingRate = emc.('Ali_samplingRate');
% used to determine the number of fiducials/patch for local area.
MOL_MASS = emc.('particleMass');
molMass = MOL_MASS.*(25/samplingRate);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Parameters I am currently experimenting with as of Jan 2018

try
  rmsScale = emc.('rmsScale');
catch
  % Larger RMS downweights the contribution of the tomogram. The use of
  % the molecular mass in MDa fits with experiments from the ribosome that
  % seem to be best ~ 4 and with HIV ~ 0.5 however it may be that the
  % number of subtomograms contributing to the average should also be
  % considered to then work back to an estimate of the SNR in the particle
  % (tomogram) volume.
  rmsScale = sqrt(MOL_MASS);
end


try
  nFiducialsPerPatch = emc.('n_fiducials_per_patch');
catch
  % TODO how smooth should the solutions really be - should multiple
  % results be run and compared?
  nFiducialsPerPatch = ceil(100./sqrt(molMass));
end

nFiducialsPerPatch = min(nFiducialsPerPatch, 32);
% Used to calc defocus values using tilt instead of manually. Convention
% diff.
flgInvertTiltAngles = 0;

% We expect our particles to be distributed throught the thickness.
% Setting to 2 can cause the alignment to silently die
n_surfaces=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
  use_PCF =  emc.('use_PCF');
catch
  use_PCF = 0;
end

if (use_PCF)
  error('The PCF scaling is not working correctly, please set use_PCF=0');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if (skip_to_the_end_and_run)
  % The fractional runs have already copied everything to cache/mapback%d,
  % so override the tmpCache.
  tmpCache = '';
else
  tmpCache= emc.('fastScratchDisk');
end

[tmpCache, flgCleanCache, CWD] = EMC_setup_tmp_cache(tmpCache, '', 'tomoCPR', true);

nGPUs = emc.('nGPUs');
pInfo = parcluster();
gpuScale=3*samplingRate;
nWorkers = min(nGPUs*gpuScale,emc.('nCpuCores')); % 18
fprintf('Using %d workers as max of %d %d*nGPUs and %d nWorkers visible\n', ...
  nWorkers,gpuScale,nGPUs*gpuScale,pInfo.NumWorkers);


load(sprintf('%s.mat', emc.('subTomoMeta')), 'subTomoMeta');
mapBackIter = subTomoMeta.currentTomoCPR;

% FIXME: rename this to something more descriptive
ctfRange = emc.('tomo_cpr_defocus_range')*10^10;
ctfInc = emc.('tomo_cpr_defocus_step')*10^10;

calcCTF = emc.('tomo_cpr_defocus_refine');


[tiltNameList, nTiltSeries] = BH_returnIncludedTilts( subTomoMeta.mapBackGeometry );

% if (multi_node_run)
%   [ nParProcesses, iterList] = BH_multi_parallelJobs(nTiltSeries, nGPUs, sizeCalc(1), emc.nCpuCores, [cycle_numerator,cycle_denominator]);
% else
%   [ nParProcesses, iterList] = BH_multi_parallelJobs(nTiltSeries, nGPUs, sizeCalc(1), emc.nCpuCores);
% end

if (multi_node_run && ~skip_to_the_end_and_run)
  nParts = ceil(nTiltSeries ./ cycle_denominator);
  tiltStart = 1+(cycle_numerator - 1)*nParts;
  nTotal = nTiltSeries;
  nTiltSeries = min(cycle_numerator*nParts,nTiltSeries);
  fprintf('Running a subset of your tiltSeries %d - %d (of %d total)\n',tiltStart,nTiltSeries,nTotal);
end

% Cycle 0 is named differently - I'll be deleting this in an overhaul of the way
% the subTomoMeta is written.
if (CYCLE)
  try
    fprintf('Using Alignment geometry %s\n',cycleNumber);
    geometry = subTomoMeta.(cycleNumber).RawAlign;
  catch
    fprintf('Using Average geometry %s\n',cycleNumber);
    geometry = subTomoMeta.(cycleNumber).Avg_geometry;
  end
else
  try
    fprintf('Using Alignment geometry %s\n',cycleNumber);
    geometry = subTomoMeta.(cycleNumber).RawAlign;
  catch
    fprintf('Using Average geometry %s\n',cycleNumber);
    geometry = subTomoMeta.(cycleNumber).geometry;
  end
end

%%%%%%%%%%% From align Raw -- should this be its own function? %%%%%%%%
% Load in the reference images.
refVol = cell(2,1);


refName = emc.('Raw_className');


classVector{1}  = emc.('Raw_classes_odd')(1,:);
classVector{2}  = emc.('Raw_classes_eve')(1,:);

nRefs = length(classVector{1});

particleMask = cell(nRefs,1);
for iGold = 1:2
  if iGold == 1
    halfSet = 'ODD';
  else
    halfSet = 'EVE';
  end
  
  imgNAME = sprintf('class_%d_Locations_Ref_%s', refName, halfSet);
  
  iHeader = getHeader(MRCImage(subTomoMeta.(cycleNumber).(imgNAME){1},0));
  sizeWindow = iHeader.nZ.*[1,1,1];
  [ refVol{iGold} ] = BH_unStackMontage4d(1:nRefs, ...
    subTomoMeta.(cycleNumber).(imgNAME){1}, ...
    subTomoMeta.(cycleNumber).(imgNAME){2},...
    sizeWindow);
  
end


% If there are multiple classes, create a volume that can be used to color
% the mapped back tomo
% Color map to reproject for something like the ribosome.
if (refName)
  flgClassAvg = 1;
else
  flgClassAvg = 0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


try
  conserveDiskSpace = emc.('conserveDiskSpace');
catch
  conserveDiskSpace = 0;
end
tiltGeometry = subTomoMeta.tiltGeometry;

% TODO: this is a bit of an aritfact, can be removed.
outCTF = '_ctf';

is_first_run = true;

mbOUT = {[tmpCache],[mapBackIter+1],'dummy'};

for iTiltSeries = tiltStart:nTiltSeries
  if (skip_to_the_end_and_run)
    continue;
  end


  
  mapBackRePrjSize = min(64,subTomoMeta.mapBackGeometry.(tiltNameList{iTiltSeries}).('tomoCprRePrjSize'))
  
  % % %   iViewGroup = subTomoMeta.mapBackGeometry.viewGroups.(tiltNameList{iTiltSeries});
  nTomograms = subTomoMeta.mapBackGeometry.(tiltNameList{iTiltSeries}).nTomos
  if nTomograms == 0
    % No points were saved after template matching so skip this tilt series
    % altogether.
    continue
  end

  skip_this_tilt_series_because_it_is_empty = false(nTomograms,1);
  
  % tomoList = fieldnames(subTomoMeta.mapBackGeometry.tomoName);
  tomoList = subTomoMeta.mapBackGeometry.(tiltNameList{iTiltSeries}).tomoList;
  tilt_filename = sprintf('%saliStacks/%s_ali%d.fixed', CWD, tiltNameList{iTiltSeries}, mapBackIter + 1);
  [~, tltName, tltExt] = fileparts(tilt_filename);
  tilt_binned_filename = sprintf('%scache/%s_bin%d%s', CWD, tltName, samplingRate, tltExt);

  
  [~,tiltBaseName,~] = fileparts(tilt_filename);
  mbOUT{3} = tiltBaseName;
  
  
  if (mapBackIter)
    localFile = sprintf('%smapBack%d/%s_ali%d_ctf.local', CWD,mapBackIter,tiltNameList{iTiltSeries},mapBackIter);
  else
    localFile = sprintf('%sfixedStacks/%s.local',CWD,tiltNameList{iTiltSeries});
  end
  
  if exist(localFile,'file')
    fprintf('Found local file %s\n', localFile);
  else
    fprintf('No local transforms requested.\n');
    localFile = 0;
  end
  

  % The model is scaled to full sampling prior to passing to tiltalign,
  % make sure the header in the synthetic stack is set appropriately.
  unsampled_pixel_size = emc.pixel_size_angstroms;
  pixel_size = unsampled_pixel_size .* samplingRate;
  
  try
    eraseMaskType = emc.('Peak_mType');
    eraseMaskRadius = emc.('Peak_mRadius') ./ pixel_size;
    % fprintf('Further restricting peak search to radius of [%f %f %f] pixels\n', eraseMaskRadius);
    eraseMask = 1;
  catch
    eraseMask = 0;
    fprintf('\n');
  end
  
  
  particle_radius = floor(max(emc.('particleRadius')./pixel_size));
  
  % TODO, is this too restricted?
  % current default peak_mask_fraction = 0.4
  peak_search_radius = floor(emc.peak_mask_fraction .* particle_radius .* [1,1]);

  % FIXME: this should be in parseParameterFile

  lowPassCutoff = emc.('tomoCprLowPass');

  % FIXME: this can also be in parseParameterFile 
  if lowPassCutoff < 2* pixel_size
    fprintf('Psych, the cutoff is being set to Nyquist\n');
    lowPassCutoff = 2*pixel_size;
  end
  
  % FIXME: this should be in parseParameterFile
  min_res_for_ctf_fitting = 10.0;
  if (calcCTF)
    try
      min_res_for_ctf_fitting = emc.('min_res_for_ctf_fitting');
    catch
    end
    
    if sqrt(2)*pixel_size > min_res_for_ctf_fitting
      fprintf('Warning the current resolution is too low to refine the defocus. Turning off this feature\n');
      calcCTF = false;
    end
  end
  
% TODO: these defaults should be re-examined
  if any(emc.tomoCPR_target_n_patches_x_y)
    % This will be re-calculated once the tilt-series size is known.
    targetPatchSize = emc.tomoCPR_target_n_patches_x_y;
  else
    targetPatchSize = ceil(max(500, ceil(2.*(particle_radius).*sqrt(nFiducialsPerPatch))));
  end
    

  % The binned stacks should already exist, if not, this will re-create it in the cache dir.
  % Note that this will also be checked when reconstructing the full 3d background tomo.
  if (samplingRate > 1)
    for iTomo = 1:nTomograms
      BH_multi_loadOrBin(tilt_filename, samplingRate, 2, false);
    end
  end
  
  % Check that an existing mapBack dir doen't exist, if so move to backup
  if exist(sprintf('mapBack%d',mapBackIter+1), 'dir')
    [y,m,d] = ymd(datetime);
    [h,mi,s] = hms(datetime);
    system(sprintf('mv mapBack%d mapBack%d_%d%0.2d%0.2d_%d_%d_%d',mapBackIter+1,mapBackIter+1,y,m,d,h,mi,floor(s)));
    clear y m d h mi s
  end

  system(sprintf('mkdir -p %smapBack%d',tmpCache,mapBackIter+1));
  system(sprintf('mkdir -p cache/mapBack%d',mapBackIter+1));
  
  % re-initialize the parpool for each tilt series to free up mem.
  if ~isempty(gcp('nocreate'))
    delete(gcp('nocreate'))
    EMC_parpool(nWorkers);
  else
    EMC_parpool(nWorkers);
  end
  fprintf('init with %d workers\n',nWorkers);
  
   
  % Get the thickest for recon
  maxZ = 0;

  % The 
  tiltHeader = getHeader(MRCImage(tilt_binned_filename, 0));
  tiltName = subTomoMeta.mapBackGeometry.tomoName.(tomoList{1}).tiltName;
  [ maxZ ] = emc_get_max_specimen_NZ( ...
                                     subTomoMeta.mapBackGeometry.tomoCoords,  ...
                                     tomoList, ...
                                     nTomograms, ...
                                     samplingRate);

  fprintf('combining thickness and shift, found a maxZ of %d\n',maxZ);
  
  % xyzproj assumes centered in Z, so add extra height for z offsets to create
  % the true "in microsope" dimension
  
  reconstruction_size = [tiltHeader.nX, tiltHeader.nY, maxZ];
  binned_specimen_origin_in_specimen_frame = emc_get_origin_index(reconstruction_size);
  avgTomo = cell(3,1);
  
  % These two are mutually exclusive for now, but not enforced.
  if (flgClassAvg)
    avgColor = zeros(reconstruction_size, 'int16');
  end
  
  if (emc.save_mapback_classes)
    avgColor = zeros(reconstruction_size, 'int16');
  end
  

  TLT = tiltGeometry.(tomoList{1});
  
  
  iRawTltName = sprintf('%smapBack%d/%s_align.rawtlt',mbOUT{1:3})
  iTiltFile = fopen(iRawTltName, 'w');
  rawTLT = sortrows(TLT(:,[1,4]),1);
  fprintf(iTiltFile,'%f\n',rawTLT(:,2)');
  fclose(iTiltFile);
  % Test this out with the full reconstruction, should enforce zeroing
  % past the first CTF zero. For now just flip blindly. No Offsets should
  % be needed.
  
  % There is a gpu clear inside that is a prob. Not sure how to handle.
  % I could run outside the loop, but that would be disk space heave
  for iRef = 1:nRefs
    refVol{1}{iRef} = gather(refVol{1}{iRef});
    refVol{2}{iRef} = gather(refVol{2}{iRef});
    particleMask{iRef} = gather(particleMask{iRef});
  end
  
  sprintf('[%d,%d]',maxZ,samplingRate);
  tiltNameList{iTiltSeries};
  
  if (use_background_estimate)
    backgroundName = sprintf('%scache/%s_%d_bin%d_backgroundEst.rec',CWD,tiltNameList{iTiltSeries},1, samplingRate);
    fprintf('In tomocpr, using background estimate %s\n\n',backgroundName);
    send_backgroundLowPassResolution = 28;

    % TODO: investigate deviations from the default, which is to shut off the phakePhasePlate and to use a backgroundLowPassResolution of 28
    % Default false, we don't apply this filter
    % if enabled, it currently only saves the filtered background estimate for visualization in addition to the normal version
    % if (emc.save_mapback_classes)
    %   BH_ctf_Correct3d(PARAMETER_FILE,sprintf('[%d,%d]',maxZ,samplingRate),tiltNameList{iTiltSeries}, 1, 3, tmpCache);
    % end
    % FIXME: calling like this does not use the surface fit for the background
    send_phakePhasePlateOption = [0,0];
    BH_ctf_Correct3d(PARAMETER_FILE,sprintf('[%d,%d]',maxZ,samplingRate),tiltNameList{iTiltSeries}, send_phakePhasePlateOption, send_backgroundLowPassResolution, tmpCache);
    
    % re-initialize the parpool for each tilt series to free up mem.
    delete(gcp('nocreate'))
    EMC_parpool(nWorkers);
  
  
    avgTomo{1} = OPEN_IMG('single',backgroundName);
    avgTomo{1} = avgTomo{1} ./ (rmsScale*rms(avgTomo{1}(:)));
    if (delete_background_estimate)
      system(sprintf('rm -f %s',backgroundName));
    end
  else
    avgTomo{1} = zeros(reconstruction_size, 'single');
  end % if (use_background_estimate)
  
  for iRef = 1:nRefs
    refVol{1}{iRef} = gpuArray(refVol{1}{iRef});
    refVol{2}{iRef} = gpuArray(refVol{2}{iRef});
    particleMask{iRef} = gpuArray(particleMask{iRef});
  end
  

  
  if (emc.save_mapback_classes)
    avgColor = zeros(reconstruction_size, 'int16');
  end
  
  if (buildTomo)
    coordOUT = fopen(sprintf('%smapBack%d/%s.coord',mbOUT{1:3}),'w');
    coordSTART = fopen(sprintf('%smapBack%d/%s.coord_start',mbOUT{1:3}),'w');
    if (emc.save_mapback_classes)
      coordCLASS = fopen(sprintf('%smapBack%d/%s.coord_class',mbOUT{1:3}),'w');
      coordCLASS_perTomo = cell(nTomograms,1);
      for iCoordClassPerTomo = 1:nTomograms
        coordClass_perTomo{iCoordClassPerTomo} = fopen(sprintf('%smapBack%d/%s.coord_class',mbOUT{1:2},tomoList{iCoordClassPerTomo}),'w');
      end
    end    
    defOUT   = fopen(sprintf('%smapBack%d/%s.defAng',mbOUT{1:3}),'w');
  end
  
  % Track the number of fiducials in order to scale the K-factor to more or less
  % aggressivley downweight outliers in the alignment
  nFidsTotal = 0;
  fidIDX = 0;
  for iTomo = 1:nTomograms
    
    TLT = tiltGeometry.(tomoList{iTomo});
    
    doseList = TLT(:,[1,11]);
    postExposure = doseList(:,2)';
    [sorted_doseList, doseIDX] = sortrows(doseList,2);
    preExposure = diff(sorted_doseList(:,2));
    preExposure = [preExposure; preExposure(end)];
    preExposure = postExposure - preExposure(doseIDX)';
    
    % Extract a "defocus file" for tilt to calculate the defocus for each
    % fiducial also considering the local alignment. If this works, I can
    % get rid of defAng
    iDefocusFileName = sprintf('%smapBack%d/%s_align.defocus',mbOUT{1:3});
    iDefocusFile = fopen(iDefocusFileName,'w');
    defTLT = sortrows(TLT(:,[1,15]),1);
    % Imod expects nanometers and underfocus positive (origin on specimen)
    % whereas I let the origin be the focal plane such that underfocus is
    % negative.
    fprintf(iDefocusFile,'%f\n',abs(defTLT(:,2)') .* 10^9);
    fclose(iDefocusFile);
    
    % We also need the transform from the microscope frame in order to
    % get an accurate defocus value. Not sure if I should be binning?
    % Additionally, we do NOT want the model for alignment in the
    % microscope frame,
    iXFName = sprintf('%smapBack%d/%s_align.XF',mbOUT{1:3});
    iXF = fopen(iXFName,'w');
    
    %     % 20190509 - I think this is royally screwing things up FIXME
    %     % Commenting this out invalidates the defocus vals
    %     xfTLT = sortrows(TLT(:,[1,7:10,2,3],1));
    %     fprintf(iXF,'%f %f %f %f %f %f\n',xfTLT(:,2:7)');
    %     fclose(iXF);
    
    xfTLT = zeros(size(TLT,1),6);
    xfTLT(:,[1,4]) = 1.0;
    fprintf(iXF,'%f %f %f %f %f %f\n',xfTLT');
    fclose(iXF);
    
    positionList = geometry.(tomoList{iTomo});
    
    positionList = positionList(positionList(:,26) ~= -9999,:);
    nFidsTotal = nFidsTotal + size(positionList,1);
    
    % Need to store tilt name/path explicity in meta deta
    tiltName    = tilt_binned_filename;
    
    
    tiltHeader = getHeader(MRCImage(tiltName,0));
    
    % FIXME, this isn't necessarily going to be the correct size
    fullTiltSizeXandY = [tiltHeader.nX,tiltHeader.nY].*samplingRate;
    
    if any(emc.tomoCPR_target_n_patches_x_y)
      % This will be re-calculated once the tilt-series size is known.
      targetPatchSize = floor([tiltHeader.nX,tiltHeader.nY].*samplingRate ./ emc.tomoCPR_target_n_patches_x_y);
      fprintf('\nUsing targetPatchSize of [%d,%d] for %s\n',targetPatchSize, tomoList{iTomo});
    end
    
    sTX = floor(tiltHeader.nX );
    sTY = floor(tiltHeader.nY );
    iTLT = floor(tiltHeader.nZ);

    tomoIdx = subTomoMeta.mapBackGeometry.tomoName.(tomoList{iTomo}).tomoIdx;
    tiltName = subTomoMeta.mapBackGeometry.tomoName.(tomoList{iTomo}).tiltName;
    reconGeometry = subTomoMeta.mapBackGeometry.tomoCoords.(tomoList{iTomo});
    tomo_origin_wrt_tilt_origin = [reconGeometry.dX_specimen_to_tomo, ...
                                    reconGeometry.dY_specimen_to_tomo, ...
                                    reconGeometry.dZ_specimen_to_tomo];              
    tomo_origin_in_tomo_frame = emc_get_origin_index([reconGeometry.NX, ...
                                                      reconGeometry.NY, ...
                                                      reconGeometry.NZ]); 
    nPrjs = size(TLT,1);
    nSubTomos = size(positionList,1);

    if (nSubTomos == 0)
      % No points were saved after template matching so skip this tilt series
      % altogether.
      skip_this_tilt_series_because_it_is_empty(iTomo) = true;
      continue;
    end

    
    
    % % TODO need to update this.
    % if (emc.save_mapback_classes)
    %   colorMap = OPEN_IMG('single',COLOR_MAP);
    %   % should be the same size as the average
      
    %   if any(size(refVol{1})-size(colorMap))
    %     error('Color map and average vol must be the same size.\n');
    %   end
    %   colorMap = colorMap(avgOrigin(1)-maxRad:avgOrigin(1)+maxRad,...
    %     avgOrigin(2)-maxRad:avgOrigin(2)+maxRad,...
    %     avgOrigin(3)-maxRad:avgOrigin(3)+maxRad);
    % end
    
    sizeAvgVol = size(refVol{1}{1});
    for iRef = 1:nRefs
      
      % FIXME: change to EMC_maskreference
      refVol{1}{iRef} = gpuArray(refVol{1}{iRef});
      refVol{2}{iRef} = gpuArray(refVol{2}{iRef});
      particleMask{iRef} = BH_mask3d('sphere',sizeAvgVol,particle_radius.*[1,1,1],[0,0,0]).* ...
        BH_mask3d(refVol{1}{iRef} + refVol{2}{iRef} ,pixel_size,'','');
    end
    
    if (buildTomo)
      
      % We need to rotate the model 90 degrees around X to match the "natural" reconstruction reference frame of imod
      % that is [x,z,-y]
      modelRot = BH_defineMatrix([0,90,0],'Bah','fwdVector');
      
      for iSubTomo = 1:nSubTomos

        subtomo_rot_matrix = reshape(positionList(iSubTomo,17:25),3,3);
        subtomo_origin_in_tomo_frame = positionList(iSubTomo,11:13);
        subtomo_origin_wrt_tilt_origin = subtomo_origin_in_tomo_frame - tomo_origin_in_tomo_frame + tomo_origin_wrt_tilt_origin;
        
        iRefIDX = 1;
        iClassIDX = positionList(iSubTomo,26);
        if (nRefs > 1)
          % Assuming generally there are fewer classes seleceted as references than there are total classes
          % For those that aren't one of the select ones, we could try to track the best matched reference from the most recent
          % alignment
          % FIXME: having a class occupancy factor would be better than just picking a random one.
          
          if ~(ismember(iClassIDX,classVector{1}) || ismember(iClassIDX,classVector{2}))
            use_class = datasample(classVector{1},1);
          else
            use_class = iClassIDX;
          end
          iRefIDX = find(classVector{1} == use_class);
        end
        
        
        % FIXME: This extra shift came from experiments with real data but is both annoying and not understood.
        subtomo_origin_wrt_tilt_origin = subtomo_origin_wrt_tilt_origin - emc.flgPreShift;

        % subTomo origin relative to reconLowerLeft
        subtomo_origin_in_sample = binned_specimen_origin_in_specimen_frame + subtomo_origin_wrt_tilt_origin./samplingRate; 
        
        % Resample a copy of the average to match the position in the tomogram
        % The third entry is a dummy, normally used to make sure at least the
        % particle was being extracted even if the surrounding density (where
        % some delocalized values may be located) are not.
        [ indVAL, padVAL, shiftVAL ] =  BH_isWindowValid(reconstruction_size, sizeAvgVol, sizeAvgVol./5, subtomo_origin_in_sample);
        
                
        if ischar(indVAL)
          fprintf('ignoring subTomo %d for out of bounds conditions.\n', iSubTomo);
        else
          if positionList(iSubTomo,7) == 1
            iAvgResamp = BH_resample3d(refVol{1}{iRefIDX},subtomo_rot_matrix',shiftVAL,'Bah','GPU','forward');
          elseif positionList(iSubTomo,7) ==2
            iAvgResamp = BH_resample3d(refVol{2}{iRefIDX},subtomo_rot_matrix',shiftVAL,'Bah','GPU','forward');
          else
            error('positionList iSubtomo %d col 7 is %d',iSubTomo,positionList(iSubTomo,7));
          end
          iMaskResamp = BH_resample3d(particleMask{iRefIDX},subtomo_rot_matrix',shiftVAL,'Bah','GPU','forward');
          
          iAvgResamp = gather(iMaskResamp.*iAvgResamp);

          
          if (emc.save_mapback_classes)
            if ~(emc.save_mapback_classes)
              iColorMap = gather(int16(iMaskResamp.* BH_resample3d(colorMap, ...
                subtomo_rot_matrix',shiftVAL,'Bah','GPU','forward')));
            else
              % Set value to class average number
              iColorMap = iMaskResamp;
              iColorMap(iColorMap < 0.05) = 0;
              iColorMap(iColorMap >= 0.05) = iRefIDX;
              iColorMap = gather(int16(iColorMap));
            end

            avgColor(indVAL(1,1):indVAL(2,1), ...
              indVAL(1,2):indVAL(2,2), ...
              indVAL(1,3):indVAL(2,3)) = avgColor(indVAL(1,1):indVAL(2,1), ...
              indVAL(1,2):indVAL(2,2), ...
              indVAL(1,3):indVAL(2,3)) + ...
              iColorMap(1+padVAL(1,1):end-padVAL(2,1),...
              1+padVAL(1,2):end-padVAL(2,2),...
              1+padVAL(1,3):end-padVAL(2,3));
            
          end
          
          try
            avgTomo{1}(indVAL(1,1):indVAL(2,1), ...
              indVAL(1,2):indVAL(2,2), ...
              indVAL(1,3):indVAL(2,3)) =  ...
              avgTomo{1}(indVAL(1,1):indVAL(2,1), ...
              indVAL(1,2):indVAL(2,2), ...
              indVAL(1,3):indVAL(2,3)) .* ...
              gather((1 -  iMaskResamp(1+padVAL(1,1):end-padVAL(2,1),... % zeros out region being replaced
              1+padVAL(1,2):end-padVAL(2,2),...
              1+padVAL(1,3):end-padVAL(2,3)))) + ...
              iAvgResamp(1+padVAL(1,1):end-padVAL(2,1),...
              1+padVAL(1,2):end-padVAL(2,2),...
              1+padVAL(1,3):end-padVAL(2,3));
          catch
            fprintf('Warning, subTomo %d appears to be out of bounds in mapBack?\n', iSubTomo);
            indVAL
            size(avgTomo{1})
            continue
          end
          
          % Reproject using tilt, so just save the 3d coords.
    
          fprintf(coordOUT,'%0.4f %0.4f %0.4f %d\n', (modelRot * subtomo_origin_wrt_tilt_origin')./samplingRate + ...
                                                    [binned_specimen_origin_in_specimen_frame(1),binned_specimen_origin_in_specimen_frame(3),binned_specimen_origin_in_specimen_frame(2)]' - ...
                                                    emc.prjVectorShift([1,2,3])', ...
                                                    fidIDX);
          
          % Save a non-rotated model with each class on its own object for visualization
          if (emc.save_mapback_classes)
            fprintf(coordCLASS,'%d 1 %0.4f %0.4f %0.4f\n', iClassIDX, subtomo_origin_wrt_tilt_origin'./samplingRate + binned_specimen_origin_in_specimen_frame'- emc.prjVectorShift');
            % Also save one for each tomo
            fprintf(coordClass_perTomo{iTomo},'%d 1 %0.4f %0.4f %0.4f\n', iClassIDX, subtomo_origin_in_tomo_frame'./samplingRate - emc.prjVectorShift');
          
            
          end

          for iPrj = 1:nPrjs
            
            iPrj_nat = find(TLT(:,1) == iPrj);
            % imod is indexing from zero
            zCoord = iPrj_nat;
            % For a positive angle, this will rotate the positive X axis farther from the focal plane (more underfocus)% For a positive angle, this will rotate the positive X axis farther from the focal plane (more underfocus)
            rTilt = BH_defineMatrix(TLT(iPrj_nat,4),'TILT','fwdVector');
            
            prjCoords = rTilt*subtomo_origin_wrt_tilt_origin';
            
            fprintf(defOUT,'%d %d %6.6e\n', fidIDX, zCoord, prjCoords(3).*unsampled_pixel_size.*10^-10+abs(TLT(iPrj_nat,15)));
            
            d1 = (abs(TLT(iPrj_nat,15)) - samplingRate.*subtomo_origin_wrt_tilt_origin(3).*unsampled_pixel_size.*10^-10) * 10^9; % Defocus value adjusted for Z coordinate in the tomogram. nm
            d2 = TLT(iPrj_nat,12)*10^9; % half astigmatism value
            
            fprintf(coordSTART,'%d %d %d %3.3f %3.3f %3.3f %3.3f %3.3f %3.3f %3.3f %3.3f %3.3f %3.3f %3.3f %3.3f %3.3f %3.3f %d\n', ...
                                fidIDX, tomoIdx,positionList(iSubTomo,4),d1,d2,180./pi.*TLT(iPrj_nat,13),reshape(subtomo_rot_matrix,1,9) , preExposure(iPrj_nat), postExposure(iPrj_nat),positionList(iSubTomo,7));
            
            % These shifts are a record of transformation from the raw data, but here
            % we are comparing with [CTF] corrected data, from which the
            % reconstructino was made directly
            
          end % loop over tilt projections
          
          fidIDX = fidIDX + 1;
        end % if condition for valid subTomo windowing
      end % loop over subtomos
    end %%%% temp condition to skip building full tomo
    
    
  end % end of loop over tomograms on this tilt-series
  
  % No subtomos remain
  if all( skip_this_tilt_series_because_it_is_empty )
    continue;
  end
    
  if (buildTomo)

    fclose(coordOUT);
    fclose(coordSTART);
    if (emc.save_mapback_classes)
      fclose(coordCLASS);
      p2m = sprintf(['point2model -sphere 6 -thick 6 -scat  ',...
        '%smapBack%d/%s.coord_class %smapBack%d/%s_classIdx.3dfid'], ...
        mbOUT{1:3},mbOUT{1:3});
      system(p2m);
      for iCoordClassPerTomo = 1:nTomograms
        fclose(coordClass_perTomo{iCoordClassPerTomo});
        p2m = sprintf(['point2model -sphere 6 -thick 6 -scat  ',...
          '%smapBack%d/%s.coord_class %smapBack%d/%s_classIdx.3dfid'], ...
          mbOUT{1:2},tomoList{iCoordClassPerTomo},mbOUT{1:2},tomoList{iCoordClassPerTomo});
        system(p2m);
      end
    end
    p2m = sprintf(['point2model -zero -circle 3 -color 0,0,255 -values -1 ',...
                  '%smapBack%d/%s.coord %smapBack%d/%s.3dfid'], ...
                  mbOUT{1:3},mbOUT{1:3});
    system(p2m);
    
    for iSave = 1
      SAVE_IMG(avgTomo{iSave},{sprintf('%smapBack%d/%s.tmpTomo%d', mbOUT{1:3},iSave), 'half'},pixel_size);
      avgTomo{iSave} = [];
    end
    clear avgTomo
    
    if (emc.save_mapback_classes || flgClassAvg)
      SAVE_IMG(single(avgColor),{sprintf('%smapBack%d/%s.tmpTomoColor', mbOUT{1:3}),'half'},pixel_size);
      clear avgColor
    end

    tmpTomoBin = 1;
    if (emc.save_mapback_classes || flgClassAvg && tmpTomoBin > 1)
      system(sprintf(['binvol -bin %d %smapBack%d/%s.tmpTomoColor ',...
        '%smapBack%d/%s.bin%dTomoColor.mrc'], ...
        tmpTomoBin,mbOUT{1:3},mbOUT{1:3},tmpTomoBin));
      system(sprintf('rm %smapBack%d/%s.tmpTomoColor ', mbOUT{1:3}));
      
    end
    
    rotSize = [tiltHeader.nX,maxZ,tiltHeader.nY];
    
    for iSave = 1
      rotCMD = sprintf(['rotatevol -angles 0,0,90 -size %d,%d,%d ',...
        '%smapBack%d/%s.tmpTomo%d %smapBack%d/%s.tmpRot%d'], ...
        rotSize, mbOUT{1:3},iSave,mbOUT{1:3},iSave);
      
      system(rotCMD);
      
      system(sprintf('rm %smapBack%d/%s.tmpTomo%d',  mbOUT{1:3},iSave));
    end
    
    clear avgTomo{1}  wgt
  end

  % % %   % It may be faster to work with a rotated vol since the reading in may cause
  % % %   % problems, but the projection is so slow, that this isn't worth dealing with
  % % %   % now.
  if (emc.n_tilt_workers > 1)
    chunkSize = ceil(sTY./emc.n_tilt_workers);
    chunkInc = zeros(emc.n_tilt_workers,3);
    for iWorker = 1:emc.n_tilt_workers-1
      chunkInc(iWorker,:) = [iWorker,(iWorker-1)*chunkSize+1,iWorker*chunkSize];
    end
    iWorker = iWorker +1 ;
    chunkInc(end, :) = [iWorker, (iWorker-1)*chunkSize+1, sTY];
  else
    chunkSize = sTY;
    chunkInc = [1,1,sTY];
  end
  

  % % % % % % %
  if (buildTomo)
    
    
    taStr = [sprintf('%f',rawTLT(1,2))];
    for iTa = 2:length(rawTLT(:,2))
      taStr = [taStr sprintf(',%f',rawTLT(iTa,2))];
    end
    
    
    
    if (localFile)
      lastLine1 = sprintf('LOCALFILE %s', localFile)
      % Used if GPU fails
    else
      lastLine1 = '';
    end
    
    if strcmpi('GPU', 'GPU')
      if (lastLine1)
        lastLine2 = 'UseGPU 0';
        lastLine3 = 'ActionIfGPUFails 2,2';
      else
        lastLine1 = 'UseGPU 0';
        lastLine2 = 'ActionIfGPUFails 2,2';
        lastLine3 = '';
      end
    else
      lastLine2 = '';
      lastLine3 = '';
    end
    
    
    % Break this up into chunks since things hang even with the
    % ActionIfGPUFails option. Try 3 times 512,256,128
    
    for iSave = 1
      keepItRunning = 1;
      outputStackName = sprintf('%smapBack%d/%s_%d_mapBack.st',mbOUT{1:3},iSave);
      
      while (keepItRunning)
        
        inc = 0:mapBackRePrjSize:sTY-1;
        if inc(end) < sTY-1
          inc = [inc,sTY-1];
        end
        nChunks = length(inc)-1;
        
        for iChunk = 1:nChunks
          if iChunk == nChunks
            extend_by = 1;
          else
            extend_by = 0;
          end
          if iChunk == 1
            if exist(outputStackName,'file')
              fprintf('removing %s\n',outputStackName);
              system(sprintf('rm %s',outputStackName));
            end
            % Special case, initialize the full sized volume and the
            % header but don't actually reproject anything.
            fprintf('Initializing volume %d/%d with size %d\n', iChunk, nChunks, mapBackRePrjSize);

            rePrjFileName = sprintf('%smapBack%d/%s_%d_rePrj.sh',mbOUT{1:3},iSave);
            reModFileName = sprintf('%smapBack%d/%s_%d_reMod.sh',mbOUT{1:3},iSave);
            reProjFile = fopen(rePrjFileName,'w');
            reModFile = fopen(reModFileName,'w');
            fprintf(reProjFile,['#!/bin/bash\n\n',...
              'tilt -StandardInput << EOF\n',...
              'input %s\n', ...
              'output %s\n', ...
              'COSINTERP 0\n', ...
              'THICKNESS %d\n', ...
              'TILTFILE %smapBack%d/%s_align.rawtlt\n', ...
              'REPROJECT %s\n', ...
              'RecFileToReproject %smapBack%d/%s.tmpRot%d\n',...
              'TOTALSLICES %d,%d\n',...
              'ZMinAndMaxReproj %d,%d\n',...
              '%s\n', ...
              '%s\n', ...
              '%s\n',...
              'EOF'],tilt_binned_filename ,outputStackName, maxZ, ...
              mbOUT{1:3},...
              taStr, mbOUT{1:3},iSave,...
              0,sTY-1,...
              -1,-1,...
              lastLine1,lastLine2,...
              lastLine3);
            
            fclose(reProjFile);
            system(sprintf('chmod a=wrx %s',rePrjFileName));
            [failedToRun,~] = system(sprintf('%s',rePrjFileName));
            
            if (failedToRun)
              error('failed to initialize reporojection %s\n',outputStackName);
            end
          end % if iChunk == 1, initialize the projection
          
     
          fprintf('Reprojecting volume %d/%d with size %d\n',...
            iChunk,nChunks,mapBackRePrjSize);
          rePrjFileName = sprintf('%smapBack%d/%s_%d_rePrj.sh',mbOUT{1:3},iSave);
          reModFileName = sprintf('%smapBack%d/%s_%d_reMod.sh',mbOUT{1:3},iSave);
          reProjFile = fopen(rePrjFileName,'w');
          reModFile = fopen(reModFileName,'w');
          fprintf(reProjFile,['#!/bin/bash\n\n',...
            'tilt -StandardInput << EOF\n',...
            'input %s\n', ...
            'output %s\n', ...
            'COSINTERP 0\n', ...
            'THICKNESS %d\n', ...
            'TILTFILE %smapBack%d/%s_align.rawtlt\n', ...
            'REPROJECT %s\n', ...
            'RecFileToReproject %smapBack%d/%s.tmpRot%d\n',... 
            'TOTALSLICES %d,%d\n',... 
            'ZMinAndMaxReproj %d,%d\n',...
            '%s\n', ...
            '%s\n', ...
            '%s\n',...
            'EOF'],tilt_binned_filename ,outputStackName, maxZ, ...
            mbOUT{1:3},...
            taStr, mbOUT{1:3},iSave,... 
            0,sTY-1,...  
            inc(iChunk),inc(iChunk+1)-1+extend_by,...
            lastLine1,lastLine2,...
            lastLine3);
  
          fclose(reProjFile);
          system(sprintf('chmod a=wrx %s',rePrjFileName));
          
          
          
          
          [failedToRun,~] = system(sprintf('%s',rePrjFileName));
          if (failedToRun)
            % Reduce size
            switch mapBackRePrjSize
              case 512
                mapBackRePrjSize = 384;
              case 384
                mapBackRePrjSize = 256;
              case 256
                mapBackRePrjSize = 192;
              case 192
                mapBackRePrjSize = 128;
              case 128
                mapBackRePrjSize = 96;
              case 96
                % header but don't actually reproject anyth
                mapBackRePrjSize = 64;
              case 64
                mapBackRePrjSize = 32;
              case 32
                mapBackRePrjSize = 16;
              case 16
                mapBackRePrjSize = 8;
              case 8
                mapBackRePrjSize = 4;
              case 4
                system(sprintf('%s > failReProj.log',rePrjFileName))
                error(['mapBackRePrjSize = 4 is still too much for the'],...
                  ['GPU which probably means your local alignments'],...
                  ['are too large'])
            end
            % Update the stored size
            subTomoMeta.mapBackGeometry.(tiltNameList{iTiltSeries}).('tomoCprRePrjSize') = mapBackRePrjSize;
            
            system(sprintf('mv %s %s.tmp',rePrjFileName,rePrjFileName));
            system(sprintf('cp %s.tmp %s',rePrjFileName,rePrjFileName));
            system(sprintf('rm %s.tmp %s',rePrjFileName,rePrjFileName));
                          system(sprintf('%s',rePrjFileName));
            
            Break out to next iter of while loop, recalculating the
            chunk size at reduced depth.
            
            break
            
          end
          
          
          if iChunk == nChunks
            keepItRunning = 0;
          end
          
          
        end % end of for loop over chunks
      end % end of while loop
    end % loop over error and masked tomo
    
    
    fprintf(reModFile,['#!/bin/bash\n\n',...
      'tilt -StandardInput << EOF\n',...
      'input %s\n', ...
      'output %smapBack%d/%s.fid\n', ...
      'COSINTERP 0\n', ...
      'THICKNESS %d\n', ...
      'TILTFILE %smapBack%d/%s_align.rawtlt \n', ...
      'DefocusFile %smapBack%d/%s_align.defocus \n', ...
      'PixelForDefocus %f,%f\n', ...
      'AngleOutputFile %smapBack%d/%s.defAngTilt\n', ...
      'AlignTransformFile %smapBack%d/%s_align.XF\n', ...
      'ProjectModel %smapBack%d/%s.3dfid\n', ...
      '%s\n',...
      '%s\n',...
      '%s\n',...
      'EOF'], ...
      tilt_binned_filename, ... % input
      mbOUT{1:3}, ... % output fiducial model
      maxZ, ... % thickness
      mbOUT{1:3},... % tilt angle file
      mbOUT{1:3},... % defocus file
      pixel_size./10, ...
      0, ... % flgInvertTiltAngles,... % Ang --> nm
      mbOUT{1:3},...
      mbOUT{1:3},...
      mbOUT{1:3},...
      lastLine1,lastLine2,...
      lastLine3);
    
    fclose(reModFile);
    system(sprintf('chmod a=wrx %s',reModFileName));
    
    
    [failedToRun,~] = system(sprintf('%s',reModFileName));
    % Sometimes the file is busy
    if (failedToRun)
      system(sprintf('mv %s %s.tmp',reModFileName,reModFileName));
      system(sprintf('cp %s.tmp %s',reModFileName,reModFileName));
      system(sprintf('rm %s.tmp',reModFileName));
      system(sprintf('%s',reModFileName));
    end
    
    % re-write the projected coords
    system(sprintf(['model2point -float -contour -zero ', '%smapBack%d/%s.fid %smapBack%d/%s.coordPrj'], mbOUT{1:3}, mbOUT{1:3}))
    
  end
  
  for iSave = 1
    % Remove the full size tomo
    system(sprintf('rm %smapBack%d/%s.tmpRot%d',mbOUT{1:3},iSave));
  end
  
  fidList = load(sprintf('%smapBack%d/%s.coordPrj',mbOUT{1:3}));
  parList = load(sprintf('%smapBack%d/%s.coord_start',mbOUT{1:3}));
  defList = load(sprintf('%smapBack%d/%s.defAngTilt',mbOUT{1:3}));
  
  %   Need to shift again from the model coordinate system
  fidList(:,[2,3]) = fidList(:,[2,3]) + repmat(emc.prjVectorShift(1:2), size(fidList,1),1);
  foundNans = sum(isnan(fidList(:,3)));
  if (foundNans)
    fprintf('\n\t\tThere are %d NaNs in the projected fiducial list %3.3f\n\n',foundNans, foundNans/size(fidList,1)*100);
    fprintf('The only confirmed case that produced this were NaNs in the fixedStacks/tiltN.local file.\n');
    error("Exiting");
  end
  
  % Give every instance of each fiducial a unique identifier.
  fidList = [1:size(fidList,1);fidList']';

  % for center of mass
  COM = 3;
  [bx,by] = ndgrid(-COM:COM,-COM:COM);

  % add optional half radius for edge case and make the padding more
  % logical, twice the particle radius, and then CTF size using mulit_iter
  % with an optimization step
  particlePad = 2.5;
  tileRadius = floor(particlePad.*particle_radius);
  tileSize = (2.*tileRadius + 1).*[1,1];
  
  CTFSIZE = BH_multi_iterator([2.*tileSize,1], 'fourier');
  
  CTFSIZE = CTFSIZE(1:2);
  ctfOrigin = emc_get_origin_index(CTFSIZE);
  
  padCTF = BH_multi_padVal(tileSize,CTFSIZE);
  ctfMask = BH_mask3d('sphere',CTFSIZE,ctfOrigin-7,[0,0],'2d');
  
  % if (eraseMask)
  %   peakMask = BH_mask3d(eraseMaskType,CTFSIZE,3.*eraseMaskRadius,[0,0],'2d');
  % else
  %   peakMask = BH_mask3d('sphere',CTFSIZE,peak_search_radius,[0,0],'2d');
  % end
  
  % peakMask(peakMask < 0.99) = 0;

  dataMask = EMC_gaussianKernel(tileSize,particle_radius,'gpu',{});
  dataMask = dataMask ./ max(dataMask(:));

  peakMask = EMC_gaussianKernel(CTFSIZE,particle_radius./3,'gpu',{});
  peakMask = peakMask ./ max(peakMask(:));
 
  
  
  
  bandPassPrj = BH_bandpass3d([sTX,sTY,1],0,0,lowPassCutoff,'cpu',pixel_size);
  
  diagnosticCell = cell(nPrjs,1);
  evalMaskCell = cell(nPrjs,1);
  
  for iPrj = 1:nPrjs
    evalMaskCell{iPrj} = zeros(gather([sTX,sTY,1]),'uint8');
  end
  
  % Any large shifts should be obvious in the original alignment, so only
  % look around +/- this value
  % globalPeak = max(2,ceil(10/pixel_size));
  % globalPeak = globalPeak + mod(globalPeak,2);
  globalPeak = floor( min(max(sTX,sTY).*0.25, 120/pixel_size) );
  globalPeak = globalPeak + mod(globalPeak,2);
  
  globalPeakMask = zeros([sTX,sTY,1],'single');
  
  globalPeakMask(binned_specimen_origin_in_specimen_frame(1) -globalPeak : binned_specimen_origin_in_specimen_frame(1) + globalPeak,...
    binned_specimen_origin_in_specimen_frame(2) -globalPeak : binned_specimen_origin_in_specimen_frame(2) + globalPeak) = 1;
  
  globalBinary = ( globalPeakMask > 0 );
  % Zero and only changed if CTF is refined.
  defocusShifts = cell(nPrjs,1);
  
  if (calcCTF)
    nToCheck = floor(ctfRange./ctfInc);
    defShiftVect = ctfInc.*[-nToCheck:nToCheck]';
  else
    nToCheck = 1;
    defShiftVect = 0;
  end
  nDefTotal = length(defShiftVect);
  
  defocusCCC = cell(nPrjs,1);
  expectedDefocusPerFiducial=cell(nPrjs,1);
  
  
  if samplingRate > 1
    tiltSeries = sprintf('%scache/%s_ali%d_bin%d.fixed',CWD,tiltName,mapBackIter+1,samplingRate);
  else
    tiltSeries = sprintf('%saliStacks/%s_ali%d.fixed',CWD,tiltName,mapBackIter+1);
  end
  
  
  
  % Optionally restrict the search to a given number of fiducials:
  nUniqueFids = numel(unique(fidList(:,2))); % I think the max val of this column should also be okay (+1)
  nFidsTotal = nUniqueFids;
  %     nFidsTotal =  sum(fidList(:,5) == 1 );
  if emc.tomoCPR_random_subset == -1 || emc.tomoCPR_random_subset > nUniqueFids
    fprintf('Using all of the %d available fiducials\n',nUniqueFids);
  else
    fprintf('Using a random subset of %d fiducials from the %d available\n',...
      emc.tomoCPR_random_subset, nUniqueFids);
    
    keepFids = datasample(0:nUniqueFids-1,emc.tomoCPR_random_subset,'Replace',false);
    fidList(~ismember(fidList(:,2),keepFids),2) = -9999;
    nFidsTotal = emc.tomoCPR_random_subset;
  end
  
  
  for iPrj = 1:nPrjs
    % I must specify the number of fiducials somehwere else, replace the
    % unique when there is time.
    defocusShifts{iPrj} = 0;
    defocusCCC{iPrj} = zeros(nDefTotal, nFidsTotal,'single','gpuArray');
    expectedDefocusPerFiducial{iPrj} = zeros(nDefTotal,nFidsTotal,'single');
  end
  
  %Put back into a natural order
  TLT = sortrows(TLT,1);
  
  if isnan(emc.k_factor_scaling)
    emc.k_factor_scaling = 10 / sqrt(nFidsTotal);
  end
  
  
  % for iPrj = 1:nPrjs % recert
  parfor iPrj = 1:nPrjs
    
    %   	    % For some reason if these mrc objects are created before the parfor
    % loop begins, they fail to load. It is fine as a regular for loop
    % though - annoying, but very little overhead. It would be nice
    % to know what is going on here.
    bhF = fourierTransformer(randn(CTFSIZE,'single','gpuArray'));
    
    iMrcObj = MRCImage(tiltSeries,0);
    iMrcObjRef = MRCImage(sprintf('%smapBack%d/%s_1_mapBack.st',mbOUT{1:3}),0);
    iMrcObjSamplingMask = MRCImage(sprintf('%saliStacks/%s_ali%d.fixed.samplingMask',CWD,tiltName,mapBackIter+1),0);
    
    
    tic
    while toc < 300
      try
        dataPrj = OPEN_IMG('single', iMrcObj,[1,sTX],[1,sTY],iPrj,'keep');
        break
      catch
        pause(1e-1)
        if ~(mod(toc,10))
          fprintf('Waiting for %d sec for dataPrj to be available\n',toc);
        end
      end
    end
    if toc == 300
      error('failed to load dataPrj %s at %d after %d tries\n', ...
        tiltSeries,iPrj,3000);
    else
      %         fprintf('loaded dataPrj %d on try %d\n',iPrj,floor(toc./0.1));
    end
    tic
    while toc < 300
      try
        refPrj  = OPEN_IMG('single', iMrcObjRef,[1,sTX],[1,sTY],iPrj,'keep');
        break
      catch
        pause(1e-1)
        if ~(mod(toc,10))
          fprintf('Waiting for %d sec for refPrj to be available\n',toc);
        end
      end
    end
    if toc == 300
      error('failed to load refPrj %s at %d after %d tries\n', ...
        tiltSeries,iPrj,3000);
    else
      %         fprintf('loaded refPrj %d on try %d\n',iPrj,floor(toc./0.1));
    end
    tic
    while toc < 300
      try
        samplingMask  = gpuArray(OPEN_IMG('single', iMrcObjSamplingMask,[],[],iPrj,'keep'));
        break
      catch
        pause(1e-1)
        if ~(mod(toc,10))
          fprintf('Waiting for %d sec for samplingMask to be available\n',toc);
        end
      end
    end
    if toc == 300
      error('failed to load samplingMask %s at %d after %d tries\n', ...
        tiltSeries,iPrj,3000);
    else
      %         fprintf('loaded refPrj %d on try %d\n',iPrj,floor(toc./0.1));
    end
    
    if (samplingRate > 1)
      samplingMask = BH_resample2d(samplingMask,[0,0,0],[0,0],'Bah','GPU','forward',1/samplingRate,[sTX,sTY]);
    else
      
    end
    
    % stored as uint8, it comes out 128,129 instead of 0,1.FIXME
    samplingMask = ( samplingMask == max(samplingMask(:)) );
    
    % In case there is any carbon or other bright shit in the periphery
    normSize = floor([256,256]./samplingRate);
    % set to even dimension
    normSize = normSize + mod(normSize,2);
    
    eSize = floor([64,64]./samplingRate);
    eSize = eSize + mod(eSize,2);
    % subtract the global mean and global center so an unsampled area mask can be made
    dataPrj = dataPrj - mean(dataPrj(:));
    dataPrj = dataPrj ./ rms(dataPrj(:));
    
    % % %       dataPrj = dataPrj - BH_movingAverage(dataPrj,normSize);
    dataAVG = BH_movingAverage(dataPrj,normSize);
    % This needs to be calculated prior to normalizing the dataPrj
    dataRMS = BH_movingRMS(dataPrj-dataAVG,eSize);
    mRms = mean(dataRMS(:));
    sRms = rms(dataRMS(:)-mRms);
    
    % FIXME
    if (emc.whitenProjections)
      whitenBP = [2*particle_radius,lowPassCutoff,pixel_size,particle_radius];
      [dataPrj,NPS] = BH_whitenNoiseSpectrum(dataPrj,'',whitenBP,1);
      % Create a matched filter.
      refPrj = refPrj ./ NPS; NPS = [];
    else
      % Local scaling
      dataPrj = dataPrj - dataAVG;
      dataPrj = gather(dataPrj ./ BH_movingRMS(dataPrj,normSize));
    end
    
    
    % Global scaling
    dataPrj = dataPrj - mean(dataPrj(:));
    dataPrj = dataPrj ./ rms(dataPrj(:));
    refPrj = refPrj - mean(refPrj(:));
    refPrj = refPrj ./ rms(refPrj(:));
    
    
    % In principle the robust fitting in imod tilt should get rid of
    % outliers, however, with a systematic error in tracking this helps
    % to avoid those outliers in the first place.
    evalMask = ( dataRMS > (mRms - 2*sRms) );
    dataRMS = [];
    
    minEval = 0.9 * (2*tileRadius)^2;
    evalMask(~samplingMask) = false;
    samplingMask = [];
    
    % FIXME: The old convention was negative. For now, assuming no support for overfocus, so if a negative value is encountered, invert it.
    mean_defocus = abs(TLT(iPrj,15));
    half_astigmatism = TLT(iPrj,12);
    angle_astigmatism = TLT(iPrj,13);
    defVect = [mean_defocus + half_astigmatism, mean_defocus - half_astigmatism, angle_astigmatism];
    
    [Hqz, HqzUnMod] = BH_ctfCalc(TLT(iPrj,16).*samplingRate,TLT(iPrj,17), ...
      TLT(iPrj,18),defVect,size(refPrj), ...
      1.*TLT(iPrj,19),-0.15);
    
    Hqz = gather(Hqz);
    
    HqzUnMod = gather(HqzUnMod);
    
    
    
    
    cccPrj = fftshift(real(ifftn(bandPassPrj.*fftn(dataPrj).*(conj(fftn(refPrj).*Hqz)))));
    %      cccPrj = fftshift(real(ifftn(bandPassPrj.*fftn(dataPrj).*...
    %    
    %                           conj(fftn(refPrj).*Hqz))));
    if (use_PCF)
      cccPrj = cccPrj .* cccPrj ./ (abs(cccPrj) + 0.1);
    end

    cccPrj = (cccPrj-min(cccPrj(:))) .* globalPeakMask;
    [~,maxPRJ] = max(cccPrj(:));
    [mRx, mRy] = ind2sub(size(cccPrj), maxPRJ);
    try
      cccPRJ = cccPrj(mRx-COM:mRx+COM, mRy-COM:mRy+COM);
    catch
      mRx
      mRy
      COM
      SAVE_IMG(MRCImage(gather(cccPrj)),'err.mrc');
      error('failed to box out the cross-correlation for image\n%s\non Projection %d\n', ...
        tiltName, iPrj);
    end
    cccPRJ = cccPRJ - min(cccPRJ(:));
    comPRJX = sum(sum(bx.*cccPRJ))./sum(cccPRJ(:));
    comPRJY = sum(sum(by.*cccPRJ))./sum(cccPRJ(:));
    
  
    
    % max index - origin + fraxtional COM
    estimated_global_offset = [mRx, mRy] - binned_specimen_origin_in_specimen_frame(1:2) + [comPRJX, comPRJY];
    glbList = fopen(sprintf('%smapBack%d/%s_%03d.global',mbOUT{1:3},iPrj),'w');
    % Add unique indicies to prevent ambiquity when comparing with paral
    fprintf(glbList,'%f  degree tilt at %f %f\n', TLT(iPrj,4),estimated_global_offset);
    dataPrj = BH_resample2d(dataPrj,[0,0,0],[estimated_global_offset,0],'Bah','GPU','inv',1,size(dataPrj));
    %mapBack%d, imshow3D(gather(fftshift(real(ifftn(fftn(dataPrj).*conj(fftn(refPrj)))))));
    
    cccPrj = fftshift(real(ifftn(bandPassPrj.*fftn(dataPrj).* abs(HqzUnMod).*...
      conj(fftn(refPrj).*Hqz))));
    %      cccPrj = fftshift(real(ifftn(bandPassPrj.*fftn(dataPrj).*...
    %    
    %                           conj(fftn(refPrj).*Hqz))));

    
    %     figure, imshow3D(gather(cccPrj));
    
    % Sanity check for geometry
    cccPrj = cccPrj - mean(cccPrj(:));
    
    diagnosticCell{iPrj} = reshape(cccPrj(globalBinary),2*globalPeak+1,2*globalPeak+1);
    cccPrj = [];
    
    % Both are ordered by fiducial (imod contour number) but are not
    % explicitly checked to correspond. Should this be done?
    
    % fid list is produced by projection of the 3dmodel using tilt with
    % the angles supplied sorted from (-) --> (+)
    wrkPrjIDX = ( fidList(:,5) == iPrj - 1 );
    wrkFid = fidList(wrkPrjIDX,:);
    wrkPar = parList(wrkPrjIDX,:);
    
    wrkDefAngTilt = defList(wrkPrjIDX,[7,6,5]);
    
    coordOUT = fopen(sprintf('%smapBack%d/%s_%03d.coordFIT',mbOUT{1:3},iPrj),'w');
    
    tmpOut = zeros([CTFSIZE,size(wrkFid,1)],'single','gpuArray');
    dXY = [0,0];
    for iFid = 1:size(wrkFid,1)
      
      if wrkFid(iFid,2) == -9999
        fprintf(coordOUT,'%d %d %0.4f %0.4f %d\n', wrkFid(iFid,1:2), [-4,-4], -9999);
        continue
      end
      
      % pixelX = wrkFid(iFid,3) - emc.pixelShift + emc.flgPostShift(1);
      % pixelY = wrkFid(iFid,4) - emc.pixelShift + emc.flgPostShift(2);
      pixelX = wrkFid(iFid,3);
      pixelY = wrkFid(iFid,4);
      
      ox = floor(pixelX+0.5) - tileRadius;
      oy = floor(pixelY+0.5) - tileRadius;

      sx = (pixelX - floor(pixelX+0.5));
      sy = (pixelY - floor(pixelY+0.5));
      
      % sx = emc.pixelMultiplier*(pixelX - floor(pixelX));
      % sy = emc.pixelMultiplier*(pixelY - floor(pixelY));
      
      %         ox = floor(wrkFid(iFid,3)) - tileRadius;
      %         oy = floor(wrkFid(iFid,4)) - tileRadius;
      oxEval = [floor(wrkFid(iFid,3) - particle_radius),floor(wrkFid(iFid,3) + particle_radius)];
      oyEval = [floor(wrkFid(iFid,4) - particle_radius),floor(wrkFid(iFid,4) + particle_radius)];
      % it would be good to try a smaller tile.
      % First check that the data are found in this given projection
      
      try
        % If any zeros values within the particle radius, do not evaluate
        iSkipEval = any(any(evalMask(oxEval(1):oxEval(2),oyEval(1):oyEval(2)) == 0));
        %           iSkipEval = sum(evalMask(oxEval(1):oxEval(2),oyEval(1):oyEval(2)),'all') < minEval;
      catch
        % If the particle was outof bounds, do not evaluate
        iSkipEval = 1;
      end
      
      if ( iSkipEval )
        %           fprintf('\nThe current fiducial %d is not sampled in this projection %d\n',iFid,iPrj);
        fprintf(coordOUT,'%d %d %0.4f %0.4f %d\n', wrkFid(iFid,1:2), [2,-2], -9999);
        continue
      end
      
      if  (ox < 1 || oy < 1 || ox +2*tileRadius > sTX || oy +2*tileRadius > sTY )
        fprintf(coordOUT,'%d %d %0.4f %0.4f %d\n', wrkFid(iFid,1:2), [2,-2], -9999);
        continue
      end
      
      dataTile = dataMask.*dataPrj(ox:ox+2.*tileRadius,oy:oy+2.*tileRadius);
      dataTile = dataTile - mean(dataTile(:));
      
      
      refTile = dataMask.*refPrj(ox:ox+2.*tileRadius,oy:oy+2.*tileRadius);
      refTile = refTile - mean(refTile(:));
      
      
      dataTile = dataTile./rms(dataTile(:));
      refTile = refTile ./ rms(refTile(:));

      if (save_diagnostic_ccf)
        SAVE_IMG(dataTile, 'dataTile.mrc');
        SAVE_IMG(refTile, 'refTile.mrc');
      end

      
      dataTile = ctfMask.*BH_padZeros3d(dataTile,'fwd',padCTF, ...
        'GPU','singleTaper');
      
      refTile = ctfMask.*BH_padZeros3d(refTile,'fwd',padCTF, ...
        'GPU','singleTaper');
      
        if (save_diagnostic_ccf)
          SAVE_IMG(dataTile, 'dataTile_pad.mrc');
          SAVE_IMG(refTile, 'refTile_pad.mrc');
        end

      
      df1 = (wrkDefAngTilt(iFid,1) + wrkPar(iFid,5)) * 10;
      df2 = (wrkDefAngTilt(iFid,1) - wrkPar(iFid,5)) * 10;
      dfA = wrkPar(iFid,6);
      
      if (calcCTF)
        dataFT = bhF.swapPhase(bhF.fwdFFT(dataTile,1,1,[1e-5,40,min(min_res_for_ctf_fitting,sqrt(2).*pixel_size),pixel_size]),'fwd');
        refFT  = conj(bhF.fwdFFT(refTile,1,1,[1e-5,40,min(min_res_for_ctf_fitting,sqrt(2).*pixel_size),pixel_size]));
      else
        dataFT = bhF.swapPhase(bhF.fwdFFT(dataTile,1,1,[1e-5,400,lowPassCutoff,pixel_size]),'fwd');
        refFT  = conj(bhF.fwdFFT(refTile,1,1,[1e-5,400,lowPassCutoff,pixel_size]));
      end
      
      bestScore = -inf;
      bestCTF = 1;
      for deltaCTF = 1:nDefTotal
        iRefCTF = refFT .* ...
          mexCTF(true,false,int16(CTFSIZE(1)),int16(CTFSIZE(2)),single(samplingRate*TLT(iPrj,16)*10^10), ...
                single(TLT(iPrj,18)*10^10),single(TLT(iPrj,17)*10^3),...
                single(df1 + defShiftVect(deltaCTF)),single(df2 + defShiftVect(deltaCTF)),single(dfA),single(TLT(iPrj,18)));
        %           try
        iRefCTF = iRefCTF ./ sqrt(2.*sum(abs(iRefCTF(1:end-bhF.invTrim,:)).^2,'all'));
        cccMap = dataFT .* iRefCTF;


        if (use_PCF)
          cccMap = cccMap .* cccMap ./ (abs(cccMap) + 0.001);
        end
    
        cccMap = peakMask.*real(bhF.invFFT(cccMap));
    
        if (save_diagnostic_ccf)
          SAVE_IMG(cccMap, 'ccfMap.mrc');
        end

        tmpOut(:,:,iFid) = cccMap;
   
        
        [maxVal,maxMap] = max(cccMap(:));
        defocusCCC{iPrj}(deltaCTF,iFid) = maxVal;
        
        
        % It might be better to reorder the search, to reduce the
        % number of times this loop is executed (assumming the value is
        % closer to the center of the defocusVector)
        if maxVal > bestScore
          bestScore = maxVal;
          bestCTF = deltaCTF;
          if ~(calcCTF)
            [mMx, mMy] = ind2sub(size(cccMap), maxMap);
            
            try
              
              % It would be good to know why this is out of bounds
              % sometimes. FIXME
              
              cccMap = cccMap(mMx-COM:mMx+COM, mMy-COM:mMy+COM);
              cccMap = cccMap - min(cccMap(:));
              
              comMapX = sum(sum(bx.*cccMap))./sum(cccMap(:));
              comMapY = sum(sum(by.*cccMap))./sum(cccMap(:));
         
              
              % peak in Map is where query is relative to ref, dXY then is the shift
              % needed to move the predicted position to the measured.
              % Data moved from a position of estimated_global_offset, so add this to dXY
              
              dXY = [mMx,mMy]+[comMapX,comMapY] - ctfOrigin(1:2)+ estimated_global_offset + [sx,sy];
            catch
              dXY = estimated_global_offset + [sx,sy]; % TODO double check me
            end
            
          end
        end
      end % End of loop over defocus values.
      
      
      if (calcCTF)
        
        [~,imDefC] = max(defocusCCC{iPrj}(:,iFid),[],1);
        dCTF = imDefC;
        %           fprintf('New best score %3.6f for defocus shift %3.3eAng\n', bestScore, defShiftVect(dCTF));
        
        dataFT = bhF.swapPhase(bhF.fwdFFT(dataTile,1,1,[1e-5,400,lowPassCutoff,pixel_size]),'fwd');
        refFT  = conj(bhF.fwdFFT(refTile,1,1,[1e-5,400,lowPassCutoff,pixel_size]));
        
        iRefCTF = refFT .* ...
          mexCTF(true,false,int16(CTFSIZE(1)),int16(CTFSIZE(2)),single(samplingRate*TLT(iPrj,16)*10^10), ...
          single(TLT(iPrj,18)*10^10),single(TLT(iPrj,17)*10^3),...
          single(df1 + defShiftVect(dCTF)),single(df2 + defShiftVect(dCTF)),single(dfA),single(TLT(iPrj,18)));
        % Renormalize
        dataFT = dataFT ./ (sum(abs(dataFT(:)).^2)./numel(dataFT));
        iRefCTF = iRefCTF ./ (sum(abs(iRefCTF(:)).^2)./numel(iRefCTF));
        
        cccMap = dataFT .* iRefCTF;
        if (use_PCF)
          cccMap = cccMap .* cccMap ./ (abs(cccMap) + 0.1);
          
        end
        cccMap = peakMask.*real(bhF.invFFT(cccMap));
        
        if (save_diagnostic_ccf)
          SAVE_IMG(cccMap, 'ccfMap_2.mrc');
        end
        [~,maxMap] = max(cccMap(:));
        
        [mMx, mMy] = ind2sub(size(cccMap), maxMap);
        
        cccMap = cccMap(mMx-COM:mMx+COM, mMy-COM:mMy+COM);
        
        cccMap = cccMap - min(cccMap(:));
        
        comMapX = sum(sum(bx.*cccMap))./sum(cccMap(:));
        comMapY = sum(sum(by.*cccMap))./sum(cccMap(:));
    
        % peak in Map is where query is relative to ref, dXY then is the shift
        % needed to move the predicted position to the measured.
        % Data moved from a position of estimated_global_offset, so add this to dXY
        
        dXY = [mMx,mMy]+[comMapX,comMapY] - ctfOrigin(1:2)+ estimated_global_offset + [sx,sy];
      end
      fprintf(coordOUT,'%d %d %0.4f %0.4f %d\n', wrkFid(iFid,1:2), dXY, wrkFid(iFid,5));
      if (save_diagnostic_ccf)
        error('save_diagnostic_ccf');
      end
    end % end of loop over fiducials

    if (calcCTF)
      [~,imDefC] = max(defocusCCC{iPrj},[],1);
      expectedDefocus = mean(defShiftVect(imDefC));
      defocusShifts{iPrj} = expectedDefocus;
      fprintf('prj %d delDef %3.3e\n',expectedDefocus);
    end
    evalMaskCell{iPrj} = uint8(evalMask); evalMask = [];
    fclose(coordOUT);
  end % end of the parfor loop

  if ( calcCTF )
    
    save(sprintf('%smapBack%d/%s%s.defShiftsMat',mbOUT{1:3},outCTF),'defocusShifts');
    defShifts = fopen(sprintf('%smapBack%d/%s%s.defShifts',mbOUT{1:3},outCTF),'w');
    defCCC = sprintf('%smapBack%d/%s%s_defCCC.mat',mbOUT{1:3},outCTF);
    save(defCCC,'defocusCCC','expectedDefocusPerFiducial');
    for iPrj = 1:nPrjs
      fprintf(defShifts,'%6.6e\n',defocusShifts{iPrj});
    end
  end

  
  clear diagnosticCell evalMaskCell

  system(sprintf('cat %smapBack%d/%s_???.coordFIT | sort -k 1 -g > %smapBack%d/%s.coordFIT',mbOUT{1:3},mbOUT{1:3}));
  system(sprintf('rm %smapBack%d/%s_???.coordFIT',mbOUT{1:3}));
  
  system(sprintf('cat %smapBack%d/%s_???.global | sort -k 1 -g > %smapBack%d/%s.global',mbOUT{1:3},mbOUT{1:3}));
  system(sprintf('rm %smapBack%d/%s_???.global ',mbOUT{1:3}));
  
  % create model tomogram for cross correlation
  fidShifts = load(sprintf('%smapBack%d/%s.coordFIT',mbOUT{1:3}));
  fidShifts = fidShifts(:,2:end);
  fidCombine = fopen(sprintf('%smapBack%d/%s.coordCombine',mbOUT{1:3}),'w');
  % Output at full sampling for tiltalign
  fidFull = fopen(sprintf('%smapBack%d/%s.coordFull',mbOUT{1:3}),'w');
  if (calcCTF)
    fidDefFull = fopen(sprintf('%smapBack%d/%s%s.defFidFull',mbOUT{1:3},outCTF),'w');
  end
  fidBin  = fopen(sprintf('%smapBack%d/%s.coordBin%d',mbOUT{1:3},samplingRate),'w');
  fidList = fidList(:,2:end);
  
  
  fCombine = [fidShifts(:,1),fidList(:,2:3)+fidShifts(:,2:3),fidShifts(:,4)];
  fprintf('\n\n%d/%d pts ignored\n\n',sum(fCombine(:,4)==-9999),size(fCombine,1));
  
  fFull = fCombine;
  fDefFull = [fCombine,zeros(size(fCombine,1),1)];
  % % % % %     fDefFull(:,2:3) = fDefFull(:,2:3).*pixel_size;
  fDefFull(:,2:3) = fDefFull(:,2:3).*samplingRate;
  for iPrj = 1:nPrjs
    % Create a file that has the X,Y,defocus positions for all fiducials
    % in each tilt to use in ctf correction.
    wrkFidIDX = ( fidList(:,4) == iPrj - 1 );
    fDefFull(wrkFidIDX,5) = defList(wrkFidIDX,7) + defocusShifts{iPrj};
  end
  
  fDefFull = fDefFull(fDefFull(:,4)~=-9999,:);
  fFull = fFull(fFull(:,4)~=-9999,:);
  fCombine = fCombine(fCombine(:,4)~=-9999,:);
  
  if (calcCTF)
    fprintf(fidDefFull,'%d %4.4f %4.4f %d %2.6e\n',fDefFull');
    fclose(fidDefFull);
  end
  
  fprintf(fidBin,'%d %4.4f %4.4f %d\n',fFull');
  fclose(fidBin);
  
  % % % % %     fFull(:,2:3) = fFull(:,2:3).*samplingRate;
    % The model ends up seeing the pixel size as 1, so even though it loads
  % properly on the full aligned stack, these coords need to be scaled by
  % the pixel size since this is the input to tiltalign.
  % fFull(:,2:3) = fFull(:,2:3).*pixel_size;

  fFull(:,2:3) = (samplingRate.*(fFull(:,2:3) - emc_get_origin_index(fullTiltSizeXandY./samplingRate)) + emc_get_origin_index(fullTiltSizeXandY));
  
  fprintf(fidCombine,'%d %4.4f %4.4f %d\n',fCombine');
  fclose(fidCombine);
  
  fprintf(fidFull,'%d %4.4f %4.4f %d\n',fFull');
  fclose(fidFull);
  % convert to model
  system(sprintf('point2model -zero -circle 3 -color 0,0,255 -ImageForCoordinates %smapBack%d/%s_1_mapBack.st %smapBack%d/%s.coordCombine %smapBack%d/%s_fit-comb.fid',mbOUT{1:3},mbOUT{1:3},mbOUT{1:3}));
  system(sprintf('point2model -zero -circle 3 -color 0,0,255 -ImageForCoordinates %s %smapBack%d/%s.coordFull %smapBack%d/%s_fit-full.fid',tilt_filename,mbOUT{1:3},mbOUT{1:3}));
  system(sprintf('point2model -zero -circle 3 -color 0,0,255 -ImageForCoordinates %smapBack%d/%s_1_mapBack.st %smapBack%d/%s.coordBin%d %smapBack%d/%s_fit-bin%d.fid',mbOUT{1:3},mbOUT{1:3},samplingRate,mbOUT{1:3},samplingRate));
  % write the com script for running tiltalign
  RotDef = 5;
  TltDef = 4;
  
  if (emc.shift_z_to_to_centroid)
    final_line1 =  'ShiftZFromOriginal';
    final_line2 =  'AxisZShift 0.0';
    final_line3 =  'LocalOutputOptions 1,1,1';

    
  else
    final_line3 =  '';
    final_line2 =  '';
    final_line1 =  'LocalOutputOptions 1,1,1';
  end
  
  mbOutAlt = mbOUT;
  tilt_script_name = iRawTltName;
  if (multi_node_run)
    mbOutAlt{1} = 'cache/';
    [~,tn2,tn3] = fileparts(iRawTltName);
    tilt_script_name = sprintf('cache/mapBack%d/%s%s',mbOUT{2},tn2,tn3);
  end

  aliCom_name = sprintf('%smapBack%d/%s.align',mbOutAlt{1:3});
  aliCom = fopen(aliCom_name,'w');

  % Testing local alignment with optimized parameters using the new imod options for leave out
  fprintf(aliCom,['%smapBack%d/%s_fit-full.fid\n',... %1
                  '%smapBack%d/%s%s.3dmod\n',... %2
                  '%smapBack%d/%s%s.resid\n',... %3
                  '%smapBack%d/%s%s.xyz\n',... %4
                  '%smapBack%d/%s%s.tlt\n',... %5
                  '%smapBack%d/%s%s.xtilt\n',... %6
                  '%smapBack%d/%s%s.tltxf\n',... %7
                  '%s\n',... %8 input tilt file
                  '%3.3f\n',... KFactorScaling %9
                  '%smapBack%d/%s%s.local\n',... OutputLocalFile %10
                  '%d\n%d\n', ...TargetPatchSizeXandY %11 12
                  '%d\n%d\n',... MinFidsTotalAndEachSurface %13 14
                  '%f\n%f\n',... MinSizeOrOverlapXandY %15 16
                  '%smapBack%d/%s.align_ta.log\n',... output log file %17 
                  '%smapBack%d/%s.align_ta_optimizer.log\n'],... output log file for optimizer %18
                  mbOutAlt{1:3},... % for ModelFile
                  mbOutAlt{1:3},outCTF,... %2
                  mbOutAlt{1:3},outCTF,... %3
                  mbOutAlt{1:3},outCTF,... %4
                  mbOutAlt{1:3},outCTF,... %5
                  mbOutAlt{1:3},outCTF,... %6
                  mbOutAlt{1:3},outCTF, ... %7
                  tilt_script_name,...
                  emc.k_factor_scaling, ...
                  mbOutAlt{1:3},outCTF, ...
                  targetPatchSize(1), ...
                  targetPatchSize(2),...
                  nFiducialsPerPatch, ...
                  floor(nFiducialsPerPatch/3),...
                  emc.min_overlap, ...
                  emc.min_overlap,...
                  mbOutAlt{1:3},...
                  mbOutAlt{1:3});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % % TODO: It looks like the output model file (3dmod) is the solved positions,
  % % but is saved at a pixel size of 1. Scaling by the sampling rate in all dimensions
  % % and then adding the origin (only for Z) places the coordinates back into the bin6 model
  % % I think we could get shifts from this
  % % TODO: could use ImageOriginXandY to accound for a diffence in origin due to binning
  % % as ImageSizeXandY given as binned size*sampling rate, which may not equal full size
  % % fullTiltSizeXandY,...
  % %   unsampled_pixel_size,unsampled_pixel_size,...
  %   % 'ImageSizeXandY %d,%d\n',...
  %   %   'ImagePixelSizeXandY %f,%f\n',...
  % fprintf(aliCom,['#!/bin/bash\n\n',...
  %   '#iTiltSeries %d\n',...
  %   'tiltalign -StandardInput << EOF\n',...
  %   'ModelFile %smapBack%d/%s_fit-full.fid\n',...
  %   'ImagesAreBinned 1\n',...
  %   'OutputModelFile %smapBack%d/%s%s.3dmod\n',...
  %   'OutputResidualFile %smapBack%d/%s%s.resid\n',...
  %   'OutputFidXYZFile	%smapBack%d/%s%s.xyz\n',...
  %   'OutputTiltFile	%smapBack%d/%s%s.tlt\n',...
  %   'OutputXAxisTiltFile	%smapBack%d/%s%s.xtilt\n',...
  %   'OutputTransformFile	%smapBack%d/%s%s.tltxf\n',...
  %   'RotationAngle	0.00\n',... % assumed to be rotated already
  %   'TiltFile	%s\n',...
  %   'SurfacesToAnalyze	%d\n',...
  %   'RotOption	%d\n',... % def solve all rotations
  %   'RotDefaultGrouping	3\n',... % if rot option --> 5 use def group size
  %   'TiltOption	%d\n',... % Tilts are harder use automapping
  %   'TiltDefaultGrouping	%d\n',...
  %   'MagOption	%d\n',... % def solve all mags
  %   'MagDefaultGrouping	%d\n',...
  %   'XStretchOption	0\n',...
  %   'SkewOption	0\n',...
  %   'BeamTiltOption	0\n',...
  %   'XTiltOption	0\n',...
  %   'ResidualReportCriterion	0.001\n',...
  %   'RobustFitting\n',...
  %   'KFactorScaling %3.3f\n',...
  %   'LocalAlignments\n',...
  %   'LocalRotOption %d\n',...
  %   'LocalRotDefaultGrouping %d\n',...
  %   'LocalTiltOption %d\n',...
  %   'LocalTiltDefaultGrouping %d\n',...
  %   'LocalMagOption %d\n',...
  %   'LocalMagDefaultGrouping %d\n',...
  %   'OutputLocalFile %smapBack%d/%s%s.local\n',...
  %   'TargetPatchSizeXandY %d,%d\n', ...
  %   'MinFidsTotalAndEachSurface %d,%d\n',...
  %   'MinSizeOrOverlapXandY %f,%f\n',...
  %   '%s\n',...
  %   '%s\n',...
  %   '%s\n',...
  %   'EOF'],...
  %   iTiltSeries,...
  %   mbOutAlt{1:3},... % for ModelFile
  %   mbOutAlt{1:3},...
  %   outCTF,...
  %   mbOutAlt{1:3},outCTF,mbOutAlt{1:3},outCTF,...
  %   mbOutAlt{1:3},outCTF,mbOutAlt{1:3},outCTF,mbOutAlt{1:3},outCTF, ...
  %   tilt_script_name,...
  %   n_surfaces, ...
  %   emc.rot_option_global, ...
  %   emc.tilt_option_global, ...
  %   emc.tilt_default_grouping_global, ...
  %   emc.mag_option_global, ...
  %   emc.mag_default_grouping_global, ...
  %   emc.k_factor_scaling, ...
  %   emc.rot_option_local, ...
  %   emc.rot_default_grouping_local, ...
  %   emc.tilt_option_local, ...
  %   emc.tilt_default_grouping_local, ...
  %   emc.mag_option_local, ...
  %   emc.mag_default_grouping_local, ...
  %   mbOutAlt{1:3},outCTF,targetPatchSize, ...
  %   targetPatchSize,...
  %   nFiducialsPerPatch, ...
  %   floor(nFiducialsPerPatch/3),...
  %   emc.min_overlap,emc.min_overlap,...
  %   final_line1,final_line2,final_line3);
  % % % % Assume that any backlash was solved well enough that there are no major
  % % % % discontinuities in the coarse alignment. Mag and rot are solved/ tilt in
  % % % % the global solution anyhow, so this shouldn't be a bit deal.
  % % % % 'SeparateGroup	1-%d\n',...
  % % % % iViewGroup,
  % %        fprintf(aliCom,'\n\ngrep -A %d  " At minimum tilt" ./mapBack%d/%s_ta.log >  ./mapBack%d/tmp.log',nPrjs+2,mbOUT{1:3},mbOUT{1:3});
  % %        fprintf(aliCom,'\nawk ''{if(NR >3) print $5}'' ./mapBack%d/tmp.log > mapBack%d/%s.mag',mbOUT{1:3},mbOUT{1:3});
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  fclose(aliCom);
  % system(sprintf('chmod a=wrx %smapBack%d/%s.align',mbOUT{1:3}));
  
  if (is_first_run)
    % if ( multi_node_run )
    %   fOUT = fopen(sprintf('%smapBack%d/runAlignments_%d_%d.sh',mbOUT{1:2},tiltStart,nTiltSeries),'w');
    %   fprintf(fOUT,['%smapBack%d/%s.align > ',...
    %     '%smapBack%d/%s.align_ta.log &\n'],mbOutAlt{1:3},mbOutAlt{1:3});
    % else
    %   fOUT = fopen(sprintf('%smapBack%d/runAlignments.sh',mbOUT{1:2}),'w');
    %   fprintf(fOUT,['#!/bin/bash\n\n%smapBack%d/%s.align > ',...
    %     '%smapBack%d/%s.align_ta.log &\n'],mbOutAlt{1:3},mbOutAlt{1:3});
    % end

    if ( multi_node_run )
      fOUT = fopen(sprintf('%smapBack%d/runAlignments_%d_%d.sh',mbOUT{1:2},tiltStart,nTiltSeries),'w');
    else
      fOUT = fopen(sprintf('%smapBack%d/runAlignments.sh',mbOUT{1:2}),'w');
      fprintf(fOUT,'#!/bin/bash\n\n');
    end
    % fprintf(fOUT,'cat %s | /scratch/etna/master_align.sh `xargs` &\n',aliCom_name);
    fprintf(fOUT,'%s\n',aliCom_name);

    % Since we send to the background in a shell, makes sure the
    % function waits on children.
    %if (iTiltSeries == nTiltSeries)
    %  fprintf(fOUT,'\nwait\n');
    %end
    fclose(fOUT);
    is_first_run = false;
  else
    if ( multi_node_run )
      fOUT = fopen(sprintf('%smapBack%d/runAlignments_%d_%d.sh',mbOUT{1:2},tiltStart,nTiltSeries),'a');
    else
      fOUT = fopen(sprintf('%smapBack%d/runAlignments.sh',mbOUT{1:2}), 'a');
    end
    % fprintf(fOUT,['%smapBack%d/%s.align > ',...
    %   '%smapBack%d/%s.align_ta.log &\n'], ...
    %   mbOutAlt{1:3},mbOutAlt{1:3});
    fprintf(fOUT,'%s\n',aliCom_name);
  
    % Since we send to the background in a shell, makes sure the
    % function waits on children.
    %if (iTiltSeries == nTiltSeries)
    %  fprintf(fOUT,'\nwait\n');
    %end
    fclose(fOUT);
  end
  
  %system(sprintf('./mapBack/%s.align > ./mapBack/%s_ta.log',TN,TN));
  
  %%%%%%%%% There is still sometimes a shift in Z, fit slope of the X shifts
  %%%%%%%%% and make cutoff compared to zshift to re-run the alignment with
  %%%%%%%%% this value entered for AxisZShift.
  
  % Until I hear from DAVID get the mag from the tilt log
  %%%system(sprintf('grep -A %d  " At minimum tilt" ./mapBack/%s_ta.log >  tmp.log',nPrjs+2,TN));
  %%%system(sprintf('awk ''{if(NR >3) print $5}'' tmp.log > mapBack/%s.mag',TN));
  %%%end %uf cibdutuib

end % loop over tilts


if ( flgRunAlignments )

  mainFile = sprintf('%smapBack%d/runAlignments.sh',mbOUT{1:2});
  altFiles = sprintf('%smapBack%d/runAlignments_*.sh',mbOUT{1:2});
  
  if (multi_node_run)
    % fOUT = fopen(mainFile,'w');
    % fprintf(fOUT,'#!/bin/bash\n\n');
    % fclose(fOUT);
    system(sprintf('rm %s && touch %s',mainFile,mainFile));
    fprintf('Combining Results from alt and main\n');
    system(sprintf('cat %s >> %s',altFiles,mainFile));
  end
  
  % fOUT = fopen(mainFile,'a');
  % fprintf(fOUT,'\nwait\n');
  % fclose(fOUT);
  
  % system(sprintf('chmod a=wrx %smapBack%d/runAlignments.sh', mbOUT{1:2}));
  % system(sprintf('%smapBack%d/runAlignments.sh', mbOUT{1:2}));
  system(sprintf('cat %smapBack%d/runAlignments.sh | parallel -j%d "cat {} | /scratch/etna/master_align.sh `xargs`"', mbOUT{1:2}, emc.nCpuCores))
end

if ( conserveDiskSpace )
  system(sprintf('rm %smapBack%d/%s_mapBack.st', mbOUT{1:3}));
end



if (tmpCache)
  if (flgRunAlignments)
    system(sprintf('mv %smapBack%d mapBack%d', mbOUT{1:2}, mbOUT{2}));
  else
    % This is an partial run
    system(sprintf('mkdir -p  cache/mapBack%d', mbOUT{2}));
    system(sprintf('mv %smapBack%d/* cache/mapBack%d', mbOUT{1:2}, mbOUT{2}));
  end
end

% if (flgCleanCache)
%   % Double check that this exists to avoid data loss.
%   checkDir = dir(tmpCache);
%   if isempty(checkDir)
%     fprintf('not removing the temp cache because it did not eval with dir\n');
%   else
%     cleanItUp = sprintf('rm  %s/*',tmpCache);
%     system(cleanItUp);
%   end
% end
% Since we've updated (potentially) mapBackRePrjSize, save the new metaData.
if (flgRunAlignments)
  subTomoMeta.currentTomoCPR =  subTomoMeta.currentTomoCPR + 1;
  
  if isfield(subTomoMeta,'tomoCPR_run_in_cycle')
    subTomoMeta.('tomoCPR_run_in_cycle') = cat(1,subTomoMeta.('tomoCPR_run_in_cycle'),...
      [subTomoMeta.currentTomoCPR,CYCLE]);
  else
    subTomoMeta.('tomoCPR_run_in_cycle') = [subTomoMeta.currentTomoCPR,CYCLE];
  end
  
  save(emc.('subTomoMeta'), 'subTomoMeta');
end

end

