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
buildTomo=1;% % % % % % %
METHOD = 'GPU';
flgRunAlignments = true;
COLOR_MAP= '0';

pBH = BH_parseParameterFile(PARAMETER_FILE);

try
  flgColorMap = pBH.('flgColorMap');
catch 
  flgColorMap = 0;
end

try
  preShift = pBH.('preShift');
catch
  preShift = [-0.5,-0.5,0.5];
end

try
  postShift = pBH.('postShift');
catch
  postShift = [-0.5,-0.5];
end

try
  prjVectorShift = pBH.('prjVectorShift')';
catch
  prjVectorShift = [0.5,0.5,1.0]';
end

try
  pixelShift = pBH.('pixelShift');
catch
  pixelShift = -1;
end

try
  pixelMultiplier = pBH.('pixelMultiplier');
catch
  pixelMultiplier = 1;
end


	
CYCLE = EMC_str2double(CYCLE);
cycle_numerator = '';
cycle_denominator ='';
skip_to_the_end_and_run = false;

if numel(CYCLE) == 3
  
  % After splitting, run the alignments while skipping everything else
  if CYCLE(2) == 0 && CYCLE(3) == 0
    skip_to_the_end_and_run = true;
    flgRunAlignments = true;
  else
    flgRunAlignments = false;
  end
  cycle_numerator = CYCLE(2);
  cycle_denominator = CYCLE(3);
  CYCLE = CYCLE(1);
  flgAltRun = 1;

else
  flgAltRun = 0; % Could just use one flag for RunAlignments and ALt Run
end


cycleNumber = sprintf('cycle%0.3u', CYCLE);



reconScaling = 1;
samplingRate = pBH.('Ali_samplingRate');
% used to determine the number of fiducials/patch for local area. 
MOL_MASS = pBH.('particleMass');
molMass = MOL_MASS.*(25/samplingRate); 

try
  tomoCPR_random_subset = pBH.('tomoCPR_randomSubset')
catch
  tomoCPR_random_subset = -1
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Parameters I am currently experimenting with as of Jan 2018

 try
   rmsScale = pBH.('rmsScale');
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
  probabilityPeakiness = pBH.('probPeakiness');
catch
  probabilityPeakiness = 0;
end

try
  useAverageDefocus = pBH.('useAverageDefocus');
catch
  % While it seems like using the per fiducial defocus max for the
  % refinement makes sense, have the default be the average of the
  % projection 
  useAverageDefocus = 0;
end

try
  whitenProjections = pBH.('whitenProjections');
catch
  whitenProjections = 0;
end

% Used to calc defocus values using tilt instead of manually. Convention
% diff.
  flgInvertTiltAngles = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Playing around with the model
n_surfaces=2;
try 
  rot_option_global = pBH.('rot_option_global'); 
catch
  rot_option_global = 1; 
end
try 
  rot_option_local = pBH.('rot_option_local'); 
catch
  rot_option_local = 1;
end
try 
  rot_default_grouping_global = pBH.('rot_default_grouping_global'); 
catch
  rot_default_grouping_global = 3;
end
try 
  rot_default_grouping_local = pBH.('rot_default_grouping_local');
catch
  rot_default_grouping_local = 3;
end
try
  mag_option_global = pBH.('mag_option_global');
catch
  mag_option_global = 1;
end
try
  mag_option_local = pBH.('mag_option_local');
catch
  mag_option_local = 1;
end
try
  mag_default_grouping_global = pBH.('mag_default_grouping_global');
catch
  mag_default_grouping_global = 5;
end
try
  mag_default_grouping_local = pBH.('mag_default_grouping_local');
catch
  mag_default_grouping_local = 5;
end
try
  tilt_option_global = pBH.('tilt_option_global');
catch
  tilt_option_global = 5;
end
try
  tilt_option_local = pBH.('tilt_option_local');
catch
  tilt_option_local = 5;
end
try
  tilt_default_grouping_global = pBH.('tilt_default_grouping_global');
catch
  tilt_default_grouping_global = 5;
end
try
  tilt_default_grouping_local = pBH.('tilt_default_grouping_local');
catch
  tilt_default_grouping_local = 5;
end
try
  nFiducialsPerPatch = pBH.('n_fiducials_per_patch');
catch
  % TODO how smooth should the solutions really be - should multiple
  % results be run and compared?
  nFiducialsPerPatch = ceil(100./sqrt(molMass));
end
try
  target_patch_size = pBH.('target_patch_size');
catch
  target_patch_size = 500;
end
try
  peak_mask_fraction = pBH.('peak_mask_fraction');
catch
  peak_mask_fraction = 0.4;
end
try
  min_overlap = pBH.('min_overlap');
catch
  min_overlap = 0.5;
end
try
  k_factor_scaling = pBH.('k_factor_scaling');
catch
  k_factor_scaling = nan;
end
try
  shift_z_to_to_centroid = pBH.('shift_z_to_to_centroid');
catch
  shift_z_to_to_centroid = true;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
  use_PCF =  pBH.('use_PCF')
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
 tmpCache= pBH.('fastScratchDisk');
end

if strcmpi(tmpCache, 'ram') 
  if isempty(getenv('EMC_CACHE_MEM'))
    fprintf('Did not find a variable for EMC_CACHE_MEM\nSkipping ram\n');
    tmpCache= '';
  else
    % I have no ideah how much is needed
    if EMC_str2double(getenv('EMC_CACHE_MEM')) < 64
      fprintf('There is only 64 Gb of cache on ramdisk, not using');
      tmpCache = '';
    else
      tmpCache=getenv('MCR_CACHE_ROOT');
      fprintf('Using the tmp EMC cache in ram at %s\n',tmpCache);
    end
  end
end

% % % nWorkers = EMC_str2double(nWORKERS)
nGPUs = pBH.('nGPUs');
pInfo = parcluster();
gpuScale=3*samplingRate;
nWorkers = min(nGPUs*gpuScale,pBH.('nCpuCores')); % 18
fprintf('Using %d workers as max of %d %d*nGPUs and %d nWorkers visible\n', ...
        nWorkers,gpuScale,nGPUs*gpuScale,pInfo.NumWorkers);
      
% Check to make sure it even exists
if isempty(dir(tmpCache))
    fprintf('\n\nIt appears your fastScratchDisk\n\t%s\ndoes not exist!\n\n',tmpCache);
    tmpCache = '';
end
if isempty(tmpCache) 
  tmpCache='cache/';
  flgCleanCache = 0;
  CWD = '';
else
  flgCleanCache = 1;
  CWD = sprintf('%s/',pwd);
  % Check for a trailing slash
  slashCheck = strsplit(tmpCache,'/');
  if isempty(slashCheck{end})
    tmpCache = sprintf('%scache/',tmpCache); % prefix for mapBack
  else
    tmpCache = sprintf('%s/cache/',tmpCache); % prefix for mapBack
  end
  
end

system(sprintf('mkdir -p %s',tmpCache));


load(sprintf('%s.mat', pBH.('subTomoMeta')), 'subTomoMeta');
mapBackIter = subTomoMeta.currentTomoCPR;




% Add error check onrange for reasonable values.
ctfRange = pBH.('tomoCprDefocusRange')*10^10; 
ctfInc = pBH.('tomoCprDefocusStep')*10^10;

calcCTF = pBH.('tomoCprDefocusRefine');


[tiltNameList, nTiltSeries] = BH_returnIncludedTilts( subTomoMeta.mapBackGeometry );
  
if (flgAltRun && ~skip_to_the_end_and_run)
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

  
refName = pBH.('Raw_className');

try
  symmetry = pBH.('symmetry');
catch
  error('You must now specify a symmetry=X parameter, where symmetry E (C1,C2..CX,O,I)');
end

classVector{1}  = pBH.('Raw_classes_odd')(1,:);
classSymmetry{1}= pBH.('Raw_classes_odd')(2,:);
classVector{2}  = pBH.('Raw_classes_eve')(1,:);
classSymmetry{2}= pBH.('Raw_classes_eve')(2,:);

nRefs = length(classVector{1})


particleMask = cell(nRefs,1);
for iGold = 1:2

  if iGold == 1
    halfSet = 'ODD';
  else
    halfSet = 'EVE';
  end
 

  imgNAME = sprintf('class_%d_Locations_REF_%s', refName, halfSet)   

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

% try
%   refNameODD = sprintf('%s_%s_class0_REF_ODD.mrc', ...
%                                              cycleNumber,pBH.('subTomoMeta'));
%   refNameEVE = sprintf('%s_%s_class0_REF_EVE.mrc', ...
%                                              cycleNumber,pBH.('subTomoMeta'));     
%   refODD = getVolume(MRCImage(refNameODD));
%   refEVE = getVolume(MRCImage(refNameEVE));                                           
% catch
%   fprintf('\nDid not find either %s or %s, trying Raw prefix\n',refNameODD,refNameEVE);
%   try
%     refNameODD = sprintf('%s_%s_class0_Raw_ODD.mrc', ...
%                                              cycleNumber,pBH.('subTomoMeta'));
%     refNameEVE = sprintf('%s_%s_class0_Raw_EVE.mrc', ...
%                                              cycleNumber,pBH.('subTomoMeta'));                                             
% 
%     refODD = getVolume(MRCImage(refNameODD));
%     refEVE = getVolume(MRCImage(refNameEVE));
%   catch
%     error('\nDid not find either %s or %s\n',refNameODD,refNameEVE)
%   end
% end

try
  conserveDiskSpace = pBH.('conserveDiskSpace');
catch
  conserveDiskSpace = 0;
end
tiltGeometry = subTomoMeta.tiltGeometry;

outCTF = '';
is_first_run = true;

mbOUT = {[tmpCache],[mapBackIter+1],'dummy'};
fprintf('\nmBOUT name is %smapBack%d/%s\n',mbOUT{1:3});
  
for iTiltSeries = tiltStart:nTiltSeries
  if (skip_to_the_end_and_run)
    continue;
  end

  tiltNameList
  mapBackRePrjSize = subTomoMeta.mapBackGeometry.(tiltNameList{iTiltSeries}).('tomoCprRePrjSize');
% % %   iViewGroup = subTomoMeta.mapBackGeometry.viewGroups.(tiltNameList{iTiltSeries});
  nTomograms = subTomoMeta.mapBackGeometry.(tiltNameList{iTiltSeries}).nTomos
  if nTomograms == 0
    % No points were saved after template matching so skip this tilt series
    % altogether.
    continue
  end
  
  
    
  tiltList = cell(nTomograms,1);
 % tomoList = fieldnames(subTomoMeta.mapBackGeometry.tomoName);
  tomoList = {};
  tomoIDX = 1;
  for iTomo = 1:size(subTomoMeta.mapBackGeometry.(tiltNameList{iTiltSeries}).coords,1)
    % This is dumb, fix it to be explicit.
    if any(subTomoMeta.mapBackGeometry.(tiltNameList{iTiltSeries}).coords(iTomo,:))
      tomoList{tomoIDX} = sprintf('%s_%d',tiltNameList{iTiltSeries},iTomo);
     
    

        tiltList{tomoIDX} = sprintf('%saliStacks/%s_ali%d.fixed',...
                                 CWD,tiltNameList{iTiltSeries},mapBackIter+1);
        outCTF='_ctf';                       
     
    % Only increment if values found.
    tomoIDX = tomoIDX + 1;
    end
    
  end
  
  [~,tiltBaseName,~] = fileparts(tiltList{1});
  mbOUT{3} = tiltBaseName;


  if (mapBackIter)
    localFile = sprintf('%smapBack%d/%s_ali%d_ctf.local', ...
                        CWD,mapBackIter,tiltNameList{iTiltSeries},mapBackIter)
  else
    localFile = sprintf('%sfixedStacks/%s.local',CWD,tiltNameList{iTiltSeries})
  end

  if exist(localFile,'file')
    fprintf('Found local file\n.');
  else
    fprintf('No local transforms requested.\n');
    localFile = 0;
  end

  % For now assume that all of these are the same - this is a shitty way to
  % handle it, but keeps things general and simple.
  reconRotation = zeros(nTomograms,3);
  % Remove this option (always assume -rx)
  for iTomo = 1:nTomograms
    reconRotation(iTomo,:) = [0,0,-90];
  end
  % The model is scaled to full sampling prior to passing to tiltalign,
  % make sure the header in the synthetic stack is set appropriately.
  fullPixelSize = pBH.('PIXEL_SIZE').*10^10;
  if pBH.('SuperResolution')
    fullPixelSize = fullPixelSize * 2;
  end
  pixelSize = fullPixelSize.*samplingRate;
  
try 
  eraseMaskType = pBH.('Peak_mType');
	eraseMaskRadius = pBH.('Peak_mRadius')./pixelSize;
  fprintf('Further restricting peak search to radius %f %f %f\n',...
          eraseMaskRadius);
  eraseMask = 1;
catch
  eraseMask = 0;
  fprintf('\n');
end


  [ ~,~,maskRadius,~ ] = BH_multi_maskCheck(pBH,'Ali',pixelSize)
   PARTICLE_RADIUS = floor(max(pBH.('particleRadius')./pixelSize));
  
  %PARTICLE_RADIUS = floor(mean(pBH.('particleRadius')./pixelSize));
  % TODO, is this too restricted?
  peakSearchRad = floor(peak_mask_fraction*PARTICLE_RADIUS.*[1,1]);
  try
    lowPassCutoff = pBH.('tomoCprLowPass');
    fprintf('Using a user supplied lowpass cutoff of %3.3f Ang\n.',...
            lowPassCutoff);
  catch
    % TODO are these range limits okay?
    lowPassCutoff = 1.5.*mean(subTomoMeta.currentResForDefocusError);
    if (lowPassCutoff < 10)
      lowPassCutoff = 10;
    elseif (lowPassCutoff > 24)
      lowPassCutoff = 24;
    end
    fprintf('Using an internatlly determined lowpass cutoff of %3.3f Ang\n.',...
            lowPassCutoff);
  end
  if lowPassCutoff < 2* pixelSize
    fprintf('Psych, the cutoff is being set to Nyquist');
    lowPassCutoff = 2*pixelSize;
  end
  
 min_res_for_ctf_fitting = 10.0; 
  if (calcCTF)
    try 
      min_res_for_ctf_fitting = pBH.('min_res_for_ctf_fitting');
    catch
    end
    
    if sqrt(2)*pixelSize > min_res_for_ctf_fitting
      fprintf('Warning the current resolution is too low to refine the defocus. Turning off this feature');
      calcCTF = false;
    end
  end

% % % % %   targetPatchSize = max(500, ceil(2.*(PARTICLE_RADIUS).*sqrt(nFiducialsPerPatch)))

  nFiducialsPerPatch = ceil(100./sqrt(molMass))
  targetPatchSize = max(500, ceil(2.*(PARTICLE_RADIUS).*sqrt(nFiducialsPerPatch)))

  
% % % 
% % %   % Check to see if this tilt has already been worked on, if so skip
% % %   aliCmdFileCheck = sprintf('%smapBack%d/%s.align',mbOUT{1:3});
% % %   if exist(aliCmdFileCheck,'file')
% % %     fprintf('\n\nFound aliCmdFileCheck, skipping rather than overwrite.\n');
% % %     continue
% % %   end



        
  if (samplingRate > 1)

    for iTomo = 1:nTomograms
          [~, tltName, tltExt] = fileparts(tiltList{iTomo});
                                
      % Resample the tilt if necessary, then modify the tilt list

        BH_multi_loadOrBin(tiltList{iTomo},-1.*samplingRate, 2);
        tiltList{iTomo} = sprintf('%scache/%s_bin%d%s', ...
                                     CWD,tltName, samplingRate,tltExt); 

          
  
  
                       
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




  
  % re-initialize the parpool for each tilt series to free up mem.
  if ~isempty(gcp('nocreate'))
    delete(gcp('nocreate'))
    EMC_parpool(nWorkers);
  else
    EMC_parpool(nWorkers);
  end
  fprintf('init with %d workers\n',nWorkers);

  outputPrefix = sprintf('%s_%s', cycleNumber, pBH.('subTomoMeta'));





  

  % Get the thickest for recon
  maxZ = 0;
  overSampleZforProjection = 1.0;
  tiltHeader = getHeader(MRCImage(tiltList{1},0));

  for iTomo = 1:nTomograms

    tomoNumber = subTomoMeta.mapBackGeometry.tomoName.(tomoList{iTomo}).tomoNumber;
    tiltName = subTomoMeta.mapBackGeometry.tomoName.(tomoList{iTomo}).tiltName;
    nZdZ = subTomoMeta.mapBackGeometry.(tiltName).coords(tomoNumber,[4,6])./samplingRate
    
    % half the size in z plus the shift back to the microscope coords.
    sZneeded = 2.*ceil(overSampleZforProjection * (nZdZ(1)/2+abs(nZdZ(2))));
    if sZneeded > maxZ
      maxZ = sZneeded;
    end

    clear tomoNumber nZdZ
  end
  maxZ = maxZ + (samplingRate*2);
  fprintf('combining thickness and shift, found a maxZ of %d\n',maxZ);

  % xyzproj assumes centered in Z, so add extra height for z offsets to create
  % the true "in microsope" dimension

    reconstructionSize = [tiltHeader.nX,tiltHeader.nY,maxZ]
    originRec = ceil((reconstructionSize+1)./2)
    avgTomo = cell(3,1);

    
    avgSampling = zeros(reconstructionSize,'uint8');
    % These two are mutually exclusive for now, but not enforced.
    if (flgClassAvg)
      avgColor = zeros(reconstructionSize, 'int16');
    end

    if (flgColorMap)
      avgColor = zeros(reconstructionSize, 'int16');
    end

  % as the projection of the 3dModel with tilt will use this file and it
  % must match the zCoords in the defAng file.
  tomoList{1}
  pause(3)
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
      
      sprintf('[%d,%d]',maxZ,samplingRate)
    tiltNameList{iTiltSeries}

    
  backgroundName = sprintf('%scache/%s_%d_bin%d_backgroundEst.rec',CWD,tiltNameList{iTiltSeries},1, samplingRate);
%   emClarity('internal','ctf','3d',PARAMETER_FILE,sprintf('[%d,%d]',maxZ,samplingRate),tiltNameList{iTiltSeries},'dummy');
  BH_ctf_Correct3d(PARAMETER_FILE,sprintf('[%d,%d]',maxZ,samplingRate),tiltNameList{iTiltSeries},'dummy');
  
   % re-initialize the parpool for each tilt series to free up mem.
   delete(gcp('nocreate'))
   EMC_parpool(nWorkers);

  
  avgTomo{1} = getVolume(MRCImage(backgroundName));
  
   system(sprintf('rm %s',backgroundName));

      for iRef = 1:nRefs
        refVol{1}{iRef} = gpuArray(refVol{1}{iRef});
        refVol{2}{iRef} = gpuArray(refVol{2}{iRef});      
        particleMask{iRef} = gpuArray(particleMask{iRef});
      end
          avgTomo{1} = avgTomo{1} ./ (overSampleZforProjection.*rmsScale*rms(avgTomo{1}(:)));
    
% % %   % Now reset the binned tilt to the non-ctf corrected. Could probably
% % %   % just temporarily rename, but for testing do this.
% % %   if (samplingRate > 1)
% % %      rmTiltName = sprintf('%scache/%s_ali%d_bin%d.fixed', ...
% % %                                    CWD,tiltNameList{iTiltSeries}, mapBackIter+1, samplingRate);
% % %     % Force removal so that a binned version of the ctf stack will be
% % %     % created
% % %     system(sprintf('rm %s',rmTiltName));
% % % 
% % %     % Resample the tilt if necessary, then modify the tilt list
% % %       BH_multi_loadOrBin(sprintf('aliStacks/%s_ali%d.fixed',tiltNameList{iTiltSeries}, mapBackIter+1),-1.*samplingRate, 2);
% % %   end
    
%     avgTomo{1} = zeros(reconstructionSize,'single');



  if (flgColorMap)
    avgColor = zeros(reconstructionSize, 'int16');
  end

  if (buildTomo)
    coordOUT = fopen(sprintf('%smapBack%d/%s.coord',mbOUT{1:3}),'w');
    coordSTART = fopen(sprintf('%smapBack%d/%s.coord_start',mbOUT{1:3}),'w');

    defOUT   = fopen(sprintf('%smapBack%d/%s.defAng',mbOUT{1:3}),'w');
  end
  
  % Track the number of fiducials in order to scale the K-factor to more or less
  % aggressivley downweight outliers in the alignment
  nFidsTotal = 0;
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
    fprintf(iDefocusFile,'%f\n',defTLT(:,2)'.*(-1*10^9));
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
    tomoNumber = subTomoMeta.mapBackGeometry.tomoName.(tomoList{iTomo}).tomoNumber;
    tiltName   = subTomoMeta.mapBackGeometry.tomoName.(tomoList{iTomo}).tiltName;
    coords = subTomoMeta.mapBackGeometry.(tiltName).coords(tomoNumber,1:4);
    
%     [ binShift, ~ ] = BH_multi_calcBinShift( coords, samplingRate);
    binShift = [0,0,0];
    positionList = positionList(positionList(:,26) ~= -9999,:);
    nFidsTotal = nFidsTotal + size(positionList,1);

    % Need to store tilt name/path explicity in meta deta
    tiltName    = tiltList{iTomo};
    

    tiltHeader = getHeader(MRCImage(tiltName,0));

    fullTiltSizeXandY = [tiltHeader.nX,tiltHeader.nY].*samplingRate;
 

    sTX = floor(tiltHeader.nX );
    sTY = floor(tiltHeader.nY );
    iTLT = floor(tiltHeader.nZ);
	% FIXME the z-dimension should be 1 right?
    originPrj = ceil(([sTX,sTY,0]+1)./2);

    tomoNumber = subTomoMeta.mapBackGeometry.tomoName.(tomoList{iTomo}).tomoNumber;
    tiltName = subTomoMeta.mapBackGeometry.tomoName.(tomoList{iTomo}).tiltName;
    reconCoords = subTomoMeta.mapBackGeometry.(tiltName).coords(tomoNumber,:);
          

%     iGPU=1;

    if (buildTomo)
%        [tomo,tomoReconCoords] = BH_multi_loadOrBuild(tomoList{iTomo}, ...
%                                               reconCoords, mapBackIter, ...
%                                               samplingRate, iGPU,reconScaling,1);
                                           
           doRecon = 0;
        doLoad = false;
        reconCoords
       [~,tomoReconCoords] = BH_multi_loadOrBuild(tomoList{iTomo}, ...
                                              reconCoords, mapBackIter, ...
                                              samplingRate, doRecon,reconScaling,...
                                              doLoad, 'tomoCPR');
         


      originVol = ceil((tomoReconCoords(1,1:3)+1)./2);

      reconShift = tomoReconCoords(2,1:3);
 
    end



    nPrjs = size(TLT,1)
    nSubTomos = size(positionList,1);

    
    % TODO need to update this.
    if (flgColorMap)
      colorMap = single(getVolume(MRCImage(COLOR_MAP)));
      % should be the same size as the average

      if any(size(refVol{1})-size(colorMap))
        error('Color map and average vol must be the same size.\n');
      end
        colorMap = colorMap(avgOrigin(1)-maxRad:avgOrigin(1)+maxRad,...
                    avgOrigin(2)-maxRad:avgOrigin(2)+maxRad,...
                    avgOrigin(3)-maxRad:avgOrigin(3)+maxRad);
    end              



    % Switch from maskRadius to particleRadius 20180129
    sizeAvgVol = size(refVol{1}{1});

    



    for iRef = 1:nRefs

	% FIXME change to EMC_maskreference
        refVol{1}{iRef} = gpuArray(refVol{1}{iRef});
        refVol{2}{iRef} = gpuArray(refVol{2}{iRef});
        particleMask{iRef} = BH_mask3d('sphere',sizeAvgVol,PARTICLE_RADIUS.*[1,1,1],[0,0,0]).* ...
                       BH_mask3d(refVol{1}{iRef} + refVol{2}{iRef} ,pixelSize,'','');

%         binaryMask = particleMask{iRef} > 0.01;
%         for rV = 1:2
%           refVol{rV}{iRef} = refVol{rV}{iRef} - mean(refVol{rV}{iRef}(binaryMask));
%           refVol{rV}{iRef} = refVol{rV}{iRef} ./ (0.5.*rms(refVol{rV}{iRef}(binaryMask)));
%           refVol{rV}{iRef} = refVol{rV}{iRef} .* particleMask{iRef};
%         end
    end

  


    if (iTomo == 1)
      fidIDX = 0;
    end

    %%%%%%%%%%
    if (buildTomo)

      modelRot = BH_defineMatrix([0,90,0],'Bah','forwardVector');

    for iSubTomo = 1:nSubTomos


% 
%       prjVector = (positionList(iSubTomo,11:13)./samplingRate + binShift) - ...
%                                           originVol + reconShift;
%                                         
      rSubTomo = reshape(positionList(iSubTomo,17:25),3,3);
      prjVector = (positionList(iSubTomo,11:13)./samplingRate) - originVol + reconShift;
                                       
      iRefIDX = 1;
      iClassIDX = 1;
      if (nRefs > 1)
        % Assuming generally there are fewer classes seleceted as references than there are total classes
        % For those that aren't on of the select ones, we could try to track the best matched reference from the most recent
        % alignment
        iClassIDX = positionList(iSubTomo,26);
        if ~(ismember(iClassIDX,classVector{1}) || ismember(iClassIDX,classVector{2}))
          iClassIDX = datasample(classVector{1},1);          
        end
        iRefIDX = find(classVector{1} == iClassIDX);
      end
      


% % %       prjVector = prjVector + [0.5,0.0,-0.5];
%       prjVector = prjVector + [0.0,0.0,1.0];
      prjVector = prjVector - preShift;
      recVector = (originPrj + [0,0,ceil((reconstructionSize(3)+1)/2)] + prjVector); % subTomo origin relative to reconLowerLeft

      %Resample a copy of the average to match the position in the tomogram
      % The third entry is a dummy, normally used to make sure at least the
      % particle was being extracted even if the surrounding density (where
      % some delocalized values may be located) are not.
      [ indVAL, padVAL, shiftVAL ] = ...
                          BH_isWindowValid(reconstructionSize, sizeAvgVol, sizeAvgVol./5, recVector);
                        
       

      if ~ischar(indVAL)                  


       if positionList(iSubTomo,7) == 1
        iAvgResamp = BH_resample3d(refVol{1}{iRefIDX},rSubTomo',shiftVAL,'Bah',METHOD,'forward');  
       elseif positionList(iSubTomo,7) ==2
         iAvgResamp = BH_resample3d(refVol{2}{iRefIDX},rSubTomo',shiftVAL,'Bah',METHOD,'forward');
       else
         error('positionList iSubtomo %d col 7 is %d',iSubTomo,positionList(iSubTomo,7));
       end
       iMaskResamp = BH_resample3d(particleMask{iRefIDX},rSubTomo',shiftVAL,'Bah',METHOD,'forward');
         


  
        iAvgResamp = gather(iMaskResamp.*iAvgResamp);

        if (flgColorMap || flgClassAvg)
          if (flgColorMap)
            iColorMap = gather(int16(iMaskResamp.* BH_resample3d(colorMap, ...
                                   rSubTomo',shiftVAL,'Bah',METHOD,'forward')));  
          else
                        % Set value to class average number
            iColorMap = iMaskResamp;
            iColorMap(iColorMap < 0.05) = 0;
            iColorMap(iColorMap >= 0.05) = iRefIDX;
            iColorMap = gather(int16(iColorMap));  
          end


          if ( flgClassAvg )

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
          fprintf('Warning, subTomo %d appears to be out of bounds in mapBack?\n');
          continue
        end
          
          
        

       
        % Reproject using tilt, so just save the 3d coords. 
       fprintf(coordOUT,'%0.4f %0.4f %0.4f %d\n', modelRot*prjVector' + [originRec(1),originRec(3),originRec(2)]'- prjVectorShift([1,3,2]), fidIDX);

        for iPrj = 1:nPrjs
          
          iPrj_nat = find(TLT(:,1) == iPrj);
          % imod is indexing from zero
          % imod is indexing from zero
          zCoord = iPrj_nat; 

          rTilt = BH_defineMatrix([90,1.*TLT(iPrj_nat,4),-90],'Bah','forwardVector');
          
          prjCoords = rTilt*prjVector';

          fprintf(defOUT,'%d %d %6.6e\n', fidIDX, zCoord, samplingRate.*prjCoords(3).*fullPixelSize.*10^-10+TLT(iPrj_nat,15));
%           d1 = -1.*((samplingRate.*prjCoords(3).*fullPixelSize.*10^-10+TLT(iPrj_nat,15)) - TLT(iPrj_nat,12))*10^10;
%           d2 = -1.*((samplingRate.*prjCoords(3).*fullPixelSize.*10^-10+TLT(iPrj_nat,15)) + TLT(iPrj_nat,12))*10^10;
          
          d1 = -1.*(samplingRate.*prjVector(3).*fullPixelSize.*10^-10+TLT(iPrj_nat,15))*10^9; % Defocus value adjusted for Z coordinate in the tomogram. nm
          d2 = TLT(iPrj_nat,12)*10^9; % half astigmatism value

          fprintf(coordSTART,'%d %d %d %3.3f %3.3f %3.3f %3.3f %3.3f %3.3f %3.3f %3.3f %3.3f %3.3f %3.3f %3.3f %3.3f %3.3f %d\n',fidIDX, tomoNumber,positionList(iSubTomo,4),d1,d2,180./pi.*TLT(iPrj_nat,13),reshape(rSubTomo,1,9) , preExposure(iPrj_nat), postExposure(iPrj_nat),positionList(iSubTomo,7));

          % These shifts are a record of transformation from the raw data, but here
          % we are comparing with [CTF] corrected data, from which the
          % reconstructino was made directly

        end % loop over tilt projections      


       fidIDX = fidIDX + 1;
      else
        fprintf('ignoring subTomo %d for out of bounds conditions.\n', iSubTomo);
        reconstructionSize
        sizeAvgVol
        sizeAvgVol./5
        recVector
      end
    end % loop over subtomos
    end %%%% temp condition to skip building full tomo
    

  end


  if (buildTomo)
   fclose(coordOUT);
   fclose(coordSTART);
   

    p2m = sprintf(['point2model -zero -circle 3 -color 0,0,255 -values -1 ',...
                   '%smapBack%d/%s.coord %smapBack%d/%s.3dfid'], ...
                   mbOUT{1:3},mbOUT{1:3})
    system(p2m);

    for iSave = 1
      SAVE_IMG(MRCImage(gather(avgTomo{iSave})),sprintf('%smapBack%d/%s.tmpTomo%d', mbOUT{1:3},iSave),pixelSize);
      avgTomo{iSave} = [];      
    end
    clear avgTomo
    
    if (flgColorMap || flgClassAvg)
      SAVE_IMG(MRCImage(gather(avgColor)),sprintf('%smapBack%d/%s.tmpTomoColor', mbOUT{1:3}),4.0);
      clear avgColor    
    end
    % If not planning on visualization, save only a binned copy of the synthetic
    % tomo.
    
%       tmpTomoBin = floor(1/samplingRate*6);
% TODO make this an adjustable parameter
%       tmpTomoBin = ceil(6/pixelSize);
%       
%       for iSave = 1:1+(3*testSubtraction)      
%         system(sprintf(['binvol -bin %d %smapBack%d/%s.tmpTomo%d ',...
%                        '%smapBack%d/%s.bin%dTomo%d.mrc'], ...
%                        tmpTomoBin,mbOUT{1:3},iSave,mbOUT{1:3},tmpTomoBin,iSave));
%       end
      tmpTomoBin = 2;
      if (flgColorMap || flgClassAvg)                 
        system(sprintf(['binvol -bin %d %smapBack%d/%s.tmpTomoColor ',...
                       '%smapBack%d/%s.bin%dTomoColor.mrc'], ...
                       tmpTomoBin,mbOUT{1:3},mbOUT{1:3},tmpTomoBin)); 
        system(sprintf('rm %smapBack%d/%s.tmpTomoColor ', mbOUT{1:3}));
          
      end
      
   


      rotSize = [tiltHeader.nX,maxZ,tiltHeader.nY]

      for iSave = 1
        rotCMD = sprintf(['rotatevol -angles 0,0,90 -size %d,%d,%d ',...
               '%smapBack%d/%s.tmpTomo%d %smapBack%d/%s.tmpRot%d'], ...
               rotSize, mbOUT{1:3},iSave,mbOUT{1:3},iSave);

        system(rotCMD);
        
        system(sprintf('rm %smapBack%d/%s.tmpTomo%d',  mbOUT{1:3},iSave));


      end

 


  % % % %   if (flgColorMap)
  % % % %     SAVE_IMG(MRCImage(avgColor),sprintf('mapBack/%s_colorMap.mrc',tiltBaseName));
  % % % %     % -90 is assumed for trim vol, so if rotate vol is used add 90
  % % % %     if (rotateVol)
  % % % %       system(sprintf('rotatevol -angles 0,0,%d mapBack/%s_colorMap.mrc mapBack/%s_colorMap.rot',reconRotation(iTomo,3)+90.0,tiltBaseName,tiltBaseName));
  % % % %       system(sprintf('mv mapBack/%s_colorMap.rot mapBack/%s_colorMap.mrc',tiltBaseName,tiltBaseName));
  % % % %       system(sprintf('rm mapBack/%s_colorMap.rot',tiltBaseName));
  % % % %     end 
  % % % %   end


    clear avgTomo{1}  wgt
  end
  % % %   % It may be faster to work with a rotated vol since the reading in may cause
  % % %   % problems, but the projection is so slow, that this isn't worth dealing with
  % % %   % now.
    if (nWorkers > 1)
      chunkSize = ceil(sTY./nWorkers);
      chunkInc = zeros(nWorkers,3);
      for iWorker = 1:nWorkers-1
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
        cpuLastLine = lastLine1;  
      else
        lastLine1 = '';
        cpuLastLine = '';
      end      

      if strcmpi(METHOD, 'GPU')
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
%         refPrj = zeros(sTX,sTY,iTLT, 'single');
        
     for iSave = 1
        keepItRunning = 1;
        outputStackName = sprintf('%smapBack%d/%s_%d_mapBack.st',mbOUT{1:3},iSave);
        
        while (keepItRunning)
     

          inc = 0:mapBackRePrjSize:sTY-1;
          if inc(end) < sTY-1
            inc = [inc,sTY-1];
          end
          inc
          nChunks = length(inc)-1;

          for iChunk = 1:nChunks
            
            if iChunk == 1
              
              if exist(outputStackName,'file')
                fprintf('removing %s\n',outputStackName);
                system(sprintf('rm %s',outputStackName));
              end
                % Special case, initialize the full sized volume and the
                % header but don't actually reproject anything.
              fprintf('Initializing volume %d/%d with size %d\n',...
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
                                  'EOF'],tiltList{1} ,outputStackName, maxZ, ...
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
            end
              
                % Special case, initialize the full sized volume and the
                % header but don't actually reproject anything.
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
                                  'EOF'],tiltList{1} ,outputStackName, maxZ, ...
                                         mbOUT{1:3},...
                                         taStr, mbOUT{1:3},iSave,...
                                         0,sTY-1,...
                                         inc(iChunk),inc(iChunk+1),...
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
%               system(sprintf('%s',rePrjFileName));                            
              
              % Break out to next iter of while loop, recalculating the
              % chunk size at reduced depth.
              
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
                          'EOF'],tiltList{1}, mbOUT{1:3}, maxZ, ...
                                 mbOUT{1:3},...
                                 mbOUT{1:3},...
                                 pixelSize./10, flgInvertTiltAngles,... % Ang --> nm
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
      system(sprintf(['model2point -float -contour -zero ',...
                      '%smapBack%d/%s.fid %smapBack%d/%s.coordPrj'],...
                      mbOUT{1:3}, mbOUT{1:3}))




  
  end 
  

    for iSave = 1
    % Remove the full size tomo
      system(sprintf('rm %smapBack%d/%s.tmpRot%d',mbOUT{1:3},iSave));  
    end

    
  fidList = load(sprintf('%smapBack%d/%s.coordPrj',mbOUT{1:3}));
  parList = load(sprintf('%smapBack%d/%s.coord_start',mbOUT{1:3}));
  defList = load(sprintf('%smapBack%d/%s.defAngTilt',mbOUT{1:3}));

  
  %   Need to shift again from the model coordinate system
  fidList(:,[2,3]) = fidList(:,[2,3]) + repmat(prjVectorShift(1:2)', size(fidList,1),1);  
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
    [bx,by] = ndgrid(-COM:COM,-COM:COM)
    % add optional half radius for edge case and make the padding more
    % logical, twice the particle radius, and then CTF size using mulit_iter
    % with an optimization step
    particlePad = 1.5;
    tileRadius = floor(particlePad.*PARTICLE_RADIUS);
    tileSize = (2.*tileRadius + 1).*[1,1];
    
    CTFSIZE = BH_multi_iterator([2.*tileSize,1], 'fourier');

    CTFSIZE = CTFSIZE(1:2);
    ctfOrigin = floor(CTFSIZE./2) + 1;

    padCTF = BH_multi_padVal(tileSize,CTFSIZE);
    ctfMask = BH_mask3d('sphere',CTFSIZE,ctfOrigin-7,[0,0],'2d');  
    
    if (eraseMask)
      peakMask = BH_mask3d(eraseMaskType,CTFSIZE,eraseMaskRadius,[0,0],'2d');
    else
      peakMask = BH_mask3d('sphere',CTFSIZE,peakSearchRad,[0,0],'2d');      
    end
    
    peakMask(peakMask < 0.99) = 0;
    
    
    
    bandPassPrj = BH_bandpass3d([sTX,sTY,1],0,0,lowPassCutoff,'cpu',pixelSize);

    diagnosticCell = cell(nPrjs,1);
    evalMaskCell = cell(nPrjs,1);

    for iPrj = 1:nPrjs
      evalMaskCell{iPrj} = zeros(gather([sTX,sTY,1]),'uint8');
    end

    % Any large shifts should be obvious in the original alignment, so only
    % look around +/- this value
    globalPeak = max(2,ceil(10/pixelSize));
    globalPeak = globalPeak + mod(globalPeak,2);
    
    globalPeakMask = zeros([sTX,sTY,1],'single');

    globalPeakMask(originPrj(1) -globalPeak : originPrj(1) + globalPeak,...
                   originPrj(2) -globalPeak : originPrj(2) + globalPeak) = 1;

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
    if tomoCPR_random_subset == -1 || tomoCPR_random_subset > nUniqueFids
      fprintf('Using all of the %d available fiducials\n',nUniqueFids);
    else
      fprintf('Using a random subset of %d fiducials from the %d available\n',...
              tomoCPR_random_subset, nUniqueFids);
      
      keepFids = datasample(0:nUniqueFids-1,tomoCPR_random_subset,'Replace',false);
      fidList(~ismember(fidList(:,2),keepFids),2) = -9999;
      nFidsTotal = tomoCPR_random_subset;
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

   if isnan(k_factor_scaling)
      k_factor_scaling = 10 / sqrt(nFidsTotal);
   end

parfor iPrj = 1:nPrjs  

%   	    % For some reason if these mrc objects are created before the parfor
	    % loop begins, they fail to load. It is fine as a regular for loop
	    % though - annoying, but very little overhead. It would be nice
	    % to know what is going on here.
      bhF = fourierTransformer(randn(CTFSIZE,'single','gpuArray'),'OddSizeOversampled');

      iMrcObj = MRCImage(tiltSeries,0);
      iMrcObjRef = MRCImage(sprintf('%smapBack%d/%s_1_mapBack.st',mbOUT{1:3}),0);
      iMrcObjSamplingMask = MRCImage(sprintf('%saliStacks/%s_ali%d.fixed.samplingMask',CWD,tiltName,mapBackIter+1),0);

   
      tic
      while toc < 300
        try
          dataPrj = single(getVolume(iMrcObj,[1,sTX],[1,sTY],iPrj,'keep'));
          break
        catch
          pause(1e-1)
          if ~(mod(toc,10))
            fprintf('Waiting for %d sec for dataPrj to be available\n',toc);
          end
        end
      end
      if toc == 300
        error('failed to load dataPrj %s at %d after %d tries\n.', ...
                                                tiltSeries,iPrj,3000);
      else
%         fprintf('loaded dataPrj %d on try %d\n',iPrj,floor(toc./0.1));
      end
      tic
      while toc < 300
        try
          refPrj  = single(getVolume(iMrcObjRef,[1,sTX],[1,sTY],iPrj,'keep'));
          break
        catch
          pause(1e-1)
          if ~(mod(toc,10))
            fprintf('Waiting for %d sec for refPrj to be available\n',toc);
          end
        end
      end
      if toc == 300
        error('failed to load refPrj %s at %d after %d tries\n.', ...
                                                tiltSeries,iPrj,3000);
      else
%         fprintf('loaded refPrj %d on try %d\n',iPrj,floor(toc./0.1));
      end   
            tic
      while toc < 300
        try
          samplingMask  = gpuArray(single(getVolume(iMrcObjSamplingMask,[],[],iPrj,'keep')));
          break
        catch
          pause(1e-1)
          if ~(mod(toc,10))
            fprintf('Waiting for %d sec for samplingMask to be available\n',toc);
          end
        end
      end
      if toc == 300
        error('failed to load samplingMask %s at %d after %d tries\n.', ...
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
      if (whitenProjections)
        whitenBP = [2*PARTICLE_RADIUS,lowPassCutoff,pixelSize,PARTICLE_RADIUS];
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
                                  
      mean_defocus = TLT(iPrj,15);          
      half_astigmatism = TLT(iPrj,12);
      angle_astigmatism = TLT(iPrj,13);
      defVect = [mean_defocus - half_astigmatism, mean_defocus + half_astigmatism, angle_astigmatism];
      
      [Hqz, HqzUnMod] = BH_ctfCalc(TLT(iPrj,16).*samplingRate,TLT(iPrj,17), ...
                                   TLT(iPrj,18),defVect,size(refPrj), ...
                                   1.*TLT(iPrj,19),-0.15);

      Hqz = gather(Hqz);
      
	HqzUnMod = gather(HqzUnMod);
      



       cccPrj = fftshift(real(ifftn(bandPassPrj.*fftn(dataPrj).* abs(HqzUnMod).*...
                              conj(fftn(refPrj).*Hqz))));
%      cccPrj = fftshift(real(ifftn(bandPassPrj.*fftn(dataPrj).*...
%                              conj(fftn(refPrj).*Hqz))));      

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




      estPeak = [mRx, mRy] - originPrj(1:2) + [comPRJX, comPRJY];
      glbList = fopen(sprintf('%smapBack%d/%s_%03d.global',mbOUT{1:3},iPrj),'w');
    % Add unique indicies to prevent ambiquity when comparing with paral
      fprintf(glbList,'%f  degree tilt at %f %f\n', TLT(iPrj,4),estPeak);
      dataPrj = BH_resample2d(dataPrj,[0,0,0],[estPeak,0],'Bah',METHOD,'inv',1,size(dataPrj));
      %mapBack%d, imshow3D(gather(fftshift(real(ifftn(fftn(dataPrj).*conj(fftn(refPrj)))))));


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

 
              
      for iFid = 1:size(wrkFid,1)
        
        if wrkFid(iFid,2) == -9999
          fprintf(coordOUT,'%d %d %0.4f %0.4f %d\n', wrkFid(iFid,1:2), [-4,-4], -9999);
          continue
        end

        pixelX = wrkFid(iFid,3) - pixelShift + postShift(1);
        pixelY = wrkFid(iFid,4) - pixelShift + postShift(2);
      
        ox = floor(pixelX) - tileRadius;
        oy = floor(pixelY) - tileRadius;     
   
        sx = pixelMultiplier*(pixelX - floor(pixelX));
        sy = pixelMultiplier*(pixelY - floor(pixelY));
      
%         ox = floor(wrkFid(iFid,3)) - tileRadius;
%         oy = floor(wrkFid(iFid,4)) - tileRadius;
        oxEval = [floor(wrkFid(iFid,3) - PARTICLE_RADIUS),floor(wrkFid(iFid,3) + PARTICLE_RADIUS)];
        oyEval = [floor(wrkFid(iFid,4) - PARTICLE_RADIUS),floor(wrkFid(iFid,4) + PARTICLE_RADIUS)];
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

        dataTile = dataPrj(ox:ox+2.*tileRadius,oy:oy+2.*tileRadius);
        dataTile = dataTile - mean(dataTile(:));
       
        
        refTile = refPrj(ox:ox+2.*tileRadius,oy:oy+2.*tileRadius);         
        refTile = refTile - mean(refTile(:));
       
        
        dataTile = dataTile./rms(dataTile(:));
        refTile = refTile ./ rms(refTile(:));

        dataTile = ctfMask.*BH_padZeros3d(dataTile,'fwd',padCTF, ...
                                                 'GPU','singleTaper'); 

        refTile = ctfMask.*BH_padZeros3d(refTile,'fwd',padCTF, ...
                                                 'GPU','singleTaper');  
                                               
                                                       
                        

        df1 = (wrkDefAngTilt(iFid,1) + wrkPar(iFid,5)) * 10;
        df2 = (wrkDefAngTilt(iFid,1) - wrkPar(iFid,5)) * 10;
        dfA = wrkPar(iFid,6);
          
        if (calcCTF)
          dataFT = bhF.swapPhase(bhF.fwdFFT(dataTile,1,1,[1e-5,40,min(min_res_for_ctf_fitting,sqrt(2).*pixelSize),pixelSize]),'fwd');
          refFT  = conj(bhF.fwdFFT(refTile,1,1,[1e-5,40,min(min_res_for_ctf_fitting,sqrt(2).*pixelSize),pixelSize]));
        else
          dataFT = bhF.swapPhase(bhF.fwdFFT(dataTile,1,1,[1e-5,400,lowPassCutoff,pixelSize]),'fwd');
          refFT  = conj(bhF.fwdFFT(refTile,1,1,[1e-5,400,lowPassCutoff,pixelSize]));         
        end
     
        bestScore = -1000000;
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
              cccMap = cccMap .* cccMap ./ (abs(cccMap) + 0.1);
            else
% % % % %               cccMap = peakMask.*real(bhF.invFFT(bhF.swapPhase(bhF.fwdFFT(dataTile,1,0,[0,300,lowPassCutoff,pixelSize]).*conj(bhF.fwdFFT(refTile,1,0) .* iCTF),'fwd')));
        
            end          
            cccMap = peakMask.*real(bhF.invFFT(cccMap));

                                    

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
                  % Data moved from a position of estPeak, so add this to dXY

                  dXY = [mMx,mMy]+[comMapX,comMapY] - ctfOrigin(1:2)+ estPeak - [sx,sy]; 
                catch
                  dXY = estPeak - [sx,sy]; % TODO double check me
                end

              end
            end
        end % End of loop over defocus values.
        
   
        if (calcCTF)

          [~,imDefC] = max(defocusCCC{iPrj}(:,iFid),[],1);
          dCTF = imDefC;
%           fprintf('New best score %3.6f for defocus shift %3.3eAng\n', bestScore, defShiftVect(dCTF));

          dataFT = bhF.swapPhase(bhF.fwdFFT(dataTile,1,1,[1e-5,400,lowPassCutoff,pixelSize]),'fwd');
          refFT  = conj(bhF.fwdFFT(refTile,1,1,[1e-5,400,lowPassCutoff,pixelSize]));         
      
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

          
          [~,maxMap] = max(cccMap(:));

          [mMx, mMy] = ind2sub(size(cccMap), maxMap);

          cccMap = cccMap(mMx-COM:mMx+COM, mMy-COM:mMy+COM);

          cccMap = cccMap - min(cccMap(:));

          comMapX = sum(sum(bx.*cccMap))./sum(cccMap(:));
          comMapY = sum(sum(by.*cccMap))./sum(cccMap(:));

          % peak in Map is where query is relative to ref, dXY then is the shift
          % needed to move the predicted position to the measured.
          % Data moved from a position of estPeak, so add this to dXY

          dXY = [mMx,mMy]+[comMapX,comMapY] - ctfOrigin(1:2)+ estPeak - [sx,sy];  
        end

              
        fprintf(coordOUT,'%d %d %0.4f %0.4f %d\n', wrkFid(iFid,1:2), dXY, wrkFid(iFid,5));


            
      
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
%     evalMaskStack = zeros(sTX,sTY,nPrjs);
%     diagnosticStack = zeros([(globalPeak.*2+1).*[1,1],nPrjs],'single');
%     for iPrj = 1:nPrjs
%       diagnosticStack(:,:,iPrj) = gather(diagnosticCell{iPrj});
%       evalMaskStack(:,:,iPrj) = int16(gather(evalMaskCell{iPrj}));
%     end
%     if (bh_global_save_tomoCPR_diagnostics)
%       diagnosticStack = zeros([(globalPeak.*2+1).*[1,1],nPrjs],'single');
%       for iPrj = 1:nPrjs
%         diagnosticStack(:,:,iPrj) = gather(diagnosticCell{iPrj});
%         evalMaskStack(:,:,iPrj) = int16(gather(evalMaskCell{iPrj}));
%       end
%     end
    
    clear diagnosticCell evalMaskCell
%     if ~(conserveDiskSpace) && bh_global_save_tomoCPR_diagnostics
%       SAVE_IMG(MRCImage(gather(diagnosticStack)),sprintf('%smapBack%d/%s_diagnostic.mrc',mbOUT{1:3}));
%       SAVE_IMG(MRCImage(gather(evalMaskStack)),sprintf('%smapBack%d/%s_evalMask.mrc',mbOUT{1:3}));
%     end
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


    % shifts/List col 1/4 should match - maybe add a check to be safe
    size(fidShifts)
    size(fidList)
   
    fCombine = [fidShifts(:,1),fidList(:,2:3)+fidShifts(:,2:3),fidShifts(:,4)];
    fprintf('\n\n%d/%d pts ignored\n\n',sum(fCombine(:,4)==-9999),size(fCombine,1));

    
    
    fFull = fCombine;
    fDefFull = [fCombine,zeros(size(fCombine,1),1)];
% % % % %     fDefFull(:,2:3) = fDefFull(:,2:3).*pixelSize;
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
    fFull(:,2:3) = fFull(:,2:3).*pixelSize;

    fprintf(fidCombine,'%d %4.4f %4.4f %d\n',fCombine');
    fclose(fidCombine);

    fprintf(fidFull,'%d %4.4f %4.4f %d\n',fFull');
    fclose(fidFull);
    % convert to model 
    system(sprintf(['point2model -zero -circle 3 -color 0,0,255 ',...
                    '%smapBack%d/%s.coordCombine %smapBack%d/%s_fit-comb.fid'],mbOUT{1:3},mbOUT{1:3}));
    system(sprintf(['point2model -zero -circle 3 -color 0,0,255 ',...
                    '%smapBack%d/%s.coordFull %smapBack%d/%s_fit-full.fid'],mbOUT{1:3},mbOUT{1:3}));
    system(sprintf(['point2model -zero -circle 3 -color 0,0,255 ',...
                    '%smapBack%d/%s.coordBin%d %smapBack%d/%s_fit-bin%d.fid'],mbOUT{1:3},samplingRate,mbOUT{1:3},samplingRate));
    % write the com script for running tiltalign
    TN =tiltBaseName
    RotDef = 5;
    TltDef = 4;
    aliCom = fopen(sprintf('%smapBack%d/%s.align',mbOUT{1:3}),'w');
% % % % %     fprintf(aliCom,['#!/bin/bash\n\n',...
% % % % %                     'tiltalign -StandardInput << EOF\n',...
% % % % %                     'ModelFile %smapBack%d/%s_fit-full.fid\n',...
% % % % %                     'ImageSizeXandY %d,%d\n',...
% % % % %                     'ImagePixelSizeXandY %f,%f\n',...
% % % % %                     'ImagesAreBinned 1\n',...
% % % % %                     'OutputModelFile %smapBack%d/%s%s.3dmod\n',...
% % % % %                     'OutputResidualFile %smapBack%d/%s%s.resid\n',...
% % % % %                     'OutputFidXYZFile	%smapBack%d/%s%s.xyz\n',...
% % % % %                     'OutputTiltFile	%smapBack%d/%s%s.tlt\n',...
% % % % %                     'OutputXAxisTiltFile	%smapBack%d/%s%s.xtilt\n',...
% % % % %                     'OutputTransformFile	%smapBack%d/%s%s.tltxf\n',...
% % % % %                     'RotationAngle	0.00\n',... % assumed to be rotated already
% % % % %                     'TiltFile	%s\n',...
% % % % %                     'SurfacesToAnalyze	2\n',...
% % % % %                     'RotOption	1\n',... % def solve all rotations
% % % % %                     'RotDefaultGrouping	3\n',... % if rot option --> 5 use def group size
% % % % %                     'TiltOption	%d\n',... % Tilts are harder use automapping
% % % % %                     'TiltDefaultGrouping	%d\n',...                 
% % % % %                     'MagOption	1\n',... % def solve all mags
% % % % %                     'MagDefaultGrouping	3\n',...
% % % % %                     'XStretchOption	0\n',...
% % % % %                     'SkewOption	0\n',...          
% % % % %                     'BeamTiltOption	0\n',...  
% % % % %                     'XTiltOption	0\n',...
% % % % %                     'ResidualReportCriterion	0.001\n',...
% % % % %                     'ShiftZFromOriginal\n',...
% % % % %                     'AxisZShift 0.0\n',...
% % % % %                     'RobustFitting\n',...
% % % % %                     'KFactorScaling %3.3f\n',...
% % % % %                     'LocalAlignments\n',...
% % % % %                     'LocalRotOption 1\n',...
% % % % %                     'LocalRotDefaultGrouping 3\n',...
% % % % %                     'LocalTiltOption %d\n',...
% % % % %                     'LocalTiltDefaultGrouping %d\n',...
% % % % %                     'LocalMagOption %d\n',...
% % % % %                     'LocalMagDefaultGrouping 5\n',...
% % % % %                     'OutputLocalFile %smapBack%d/%s%s.local\n',...
% % % % %                     'TargetPatchSizeXandY %d,%d\n', ...
% % % % %                     'MinFidsTotalAndEachSurface %d,%d\n',...
% % % % %                     'MinSizeOrOverlapXandY 0.5,0.5\n',...
% % % % %                     'LocalOutputOptions 1,1,1\n', ...
% % % % %                     'EOF'],mbOUT{1:3},fullTiltSizeXandY,...
% % % % %                            fullPixelSize,fullPixelSize,...
% % % % %                            mbOUT{1:3},outCTF,mbOUT{1:3},outCTF,mbOUT{1:3},outCTF,...
% % % % %                            mbOUT{1:3},outCTF,mbOUT{1:3},outCTF,mbOUT{1:3},outCTF, ...
% % % % %                            iRawTltName,tiltAliOption(1:2),...
% % % % %                            10 / sqrt(nFidsTotal),tiltAliOption(3:4),flgLocalMag, ...
% % % % %                            mbOUT{1:3},outCTF,targetPatchSize,targetPatchSize,...
% % % % %                            nFiducialsPerPatch,floor(nFiducialsPerPatch/3));

    if (shift_z_to_to_centroid)
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
    if (flgAltRun)
      mbOutAlt{1} = 'cache/';
      [~,tn2,tn3] = fileparts(iRawTltName);
      tilt_script_name = sprintf('cache/mapBack%d/%s%s',mbOUT{2},tn2,tn3);
    end
    
    fprintf(aliCom,['#!/bin/bash\n\n',...
                    '#iTiltSeries %d\n',...
                    'tiltalign -StandardInput << EOF\n',...
                    'ModelFile %smapBack%d/%s_fit-full.fid\n',...
                    'ImageSizeXandY %d,%d\n',...
                    'ImagePixelSizeXandY %f,%f\n',...
                    'ImagesAreBinned 1\n',...
                    'OutputModelFile %smapBack%d/%s%s.3dmod\n',...
                    'OutputResidualFile %smapBack%d/%s%s.resid\n',...
                    'OutputFidXYZFile	%smapBack%d/%s%s.xyz\n',...
                    'OutputTiltFile	%smapBack%d/%s%s.tlt\n',...
                    'OutputXAxisTiltFile	%smapBack%d/%s%s.xtilt\n',...
                    'OutputTransformFile	%smapBack%d/%s%s.tltxf\n',...
                    'RotationAngle	0.00\n',... % assumed to be rotated already
                    'TiltFile	%s\n',...
                    'SurfacesToAnalyze	%d\n',...
                    'RotOption	%d\n',... % def solve all rotations
                    'RotDefaultGrouping	3\n',... % if rot option --> 5 use def group size
                    'TiltOption	%d\n',... % Tilts are harder use automapping
                    'TiltDefaultGrouping	%d\n',...                 
                    'MagOption	%d\n',... % def solve all mags
                    'MagDefaultGrouping	%d\n',...
                    'XStretchOption	0\n',...
                    'SkewOption	0\n',...          
                    'BeamTiltOption	0\n',...  
                    'XTiltOption	0\n',...
                    'ResidualReportCriterion	0.001\n',...
                    'RobustFitting\n',...
                    'KFactorScaling %3.3f\n',...
                    'LocalAlignments\n',...
                    'LocalRotOption %d\n',...
                    'LocalRotDefaultGrouping %d\n',...
                    'LocalTiltOption %d\n',...
                    'LocalTiltDefaultGrouping %d\n',...
                    'LocalMagOption %d\n',...
                    'LocalMagDefaultGrouping %d\n',...
                    'OutputLocalFile %smapBack%d/%s%s.local\n',...
                    'TargetPatchSizeXandY %d,%d\n', ...
                    'MinFidsTotalAndEachSurface %d,%d\n',...
                    'MinSizeOrOverlapXandY %f,%f\n',...
                    '%s\n',...
                    '%s\n',...
                    '%s\n',...
                    'EOF'],iTiltSeries,mbOutAlt{1:3},fullTiltSizeXandY,...
                           fullPixelSize,fullPixelSize,...
                           mbOutAlt{1:3},outCTF,mbOutAlt{1:3},outCTF,mbOutAlt{1:3},outCTF,...
                           mbOutAlt{1:3},outCTF,mbOutAlt{1:3},outCTF,mbOutAlt{1:3},outCTF, ...
                           tilt_script_name,n_surfaces,rot_option_global,...
                           tilt_option_global,tilt_default_grouping_global,...
                           mag_option_global,mag_default_grouping_global,...
                           k_factor_scaling,...
                           rot_option_local,rot_default_grouping_local,...
                           tilt_option_local,tilt_default_grouping_local,...
                           mag_option_local,mag_default_grouping_local, ...
                           mbOutAlt{1:3},outCTF,targetPatchSize,targetPatchSize,...
                           nFiducialsPerPatch,floor(nFiducialsPerPatch/3),...
                           min_overlap,min_overlap,...
                           final_line1,final_line2,final_line3);
                         
% % % Assume that any backlash was solved well enough that there are no major
% % % discontinuities in the coarse alignment. Mag and rot are solved/ tilt in
% % % the global solution anyhow, so this shouldn't be a bit deal.                         
% % % 'SeparateGroup	1-%d\n',... 
% % % iViewGroup,
%        fprintf(aliCom,'\n\ngrep -A %d  " At minimum tilt" ./mapBack%d/%s_ta.log >  ./mapBack%d/tmp.log',nPrjs+2,mbOUT{1:3},mbOUT{1:3});
%        fprintf(aliCom,'\nawk ''{if(NR >3) print $5}'' ./mapBack%d/tmp.log > mapBack%d/%s.mag',mbOUT{1:3},mbOUT{1:3});

        fclose(aliCom);
        system(sprintf('chmod a=wrx %smapBack%d/%s.align',mbOUT{1:3}));
        
        if (is_first_run)
          if ( flgAltRun )
            fOUT = fopen(sprintf('%smapBack%d/runAlignments_%d_%d.sh',mbOUT{1:2},tiltStart,nTiltSeries),'w');
            fprintf(fOUT,['%smapBack%d/%s.align > ',...
                        '%smapBack%d/%s.align_ta.log &\n'],mbOutAlt{1:3},mbOutAlt{1:3});
          else
            fOUT = fopen(sprintf('%smapBack%d/runAlignments.sh',mbOUT{1:2}),'w');
            fprintf(fOUT,['#!/bin/bash\n\n%smapBack%d/%s.align > ',...
                        '%smapBack%d/%s.align_ta.log &\n'],mbOutAlt{1:3},mbOutAlt{1:3});
          end

          % Since we send to the background in a shell, makes sure the
          % function waits on children.
          %if (iTiltSeries == nTiltSeries)
          %  fprintf(fOUT,'\nwait\n');
          %end
          fclose(fOUT);
          is_first_run = false;
        else
          if ( flgAltRun )
            fOUT = fopen(sprintf('%smapBack%d/runAlignments_%d_%d.sh',mbOUT{1:2},tiltStart,nTiltSeries),'a');
          else
            fOUT = fopen(sprintf('%smapBack%d/runAlignments.sh',mbOUT{1:2}), 'a');
          end
          fprintf(fOUT,['%smapBack%d/%s.align > ',...
                        '%smapBack%d/%s.align_ta.log &\n'], ...
                        mbOutAlt{1:3},mbOutAlt{1:3});  
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
  
  
  if (flgAltRun)
    fOUT = fopen(mainFile,'w');
    fprintf(fOUT,'#!/bin/bash\n\n');
    fclose(fOUT);
    
    fprintf('Combining Results from alt and main\n');
    system(sprintf('cat %s >> %s',altFiles,mainFile));
  end
 
  
  fOUT = fopen(mainFile,'a');
  fprintf(fOUT,'\nwait\n');
  fclose(fOUT);
  
  system(sprintf('chmod a=wrx %smapBack%d/runAlignments.sh', mbOUT{1:2}));
  
  system(sprintf('%smapBack%d/runAlignments.sh', mbOUT{1:2}));
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

  save(pBH.('subTomoMeta'), 'subTomoMeta');
end

end

