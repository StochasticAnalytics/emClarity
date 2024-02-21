function [  ] = BH_ctf_Correct3d( PARAMETER_FILE, varargin )
%Apply a full 3d CTF ala Jensen and Kornberg approach
%   Run in a few new modes
%   emClarity ctf 3d paramN.m = Normal, reconstruction for subTomo

%   emClarity ctf 3d paramN.m [THICKNESS, BINNING] = use this for picking
%   subRegions
%
%   emClarity ctf 3d paramN.m 'templateMatching' = use this to make recs
%   for higher res template matching (in the works)

% Read in 2dCtf stacks to trouble shoot
PosControl2d=0;
emc = BH_parseParameterFile(PARAMETER_FILE);

% Apply a Wiener filter with this many zeros during Ctf multiplication
global bh_global_turn_on_phase_plate
masterTM = struct();
resTarget = 15;

% TODO remove thise params
tiltWeight = [0.2,0];
shiftDefocusOrigin = 1;
tiltStart = 1;

try
  flgEraseBeads_aferCTF = emc.('erase_beads_after_ctf');
catch
  flgEraseBeads_aferCTF = false; % If false they SHOULD be erased in ctf estimate/update, but since the user could change parameter, include here.
end

% Test David's new super sampling in reconstruction. No check that this
% version (currently 4.10.40) is properly sourced.
try
  super_sample = emc.('super_sample');
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
end


try
  expand_lines = emc.('expand_lines');
  if isempty(super_sample) || expand_lines == false
    expand_lines = '';
  else
    expand_lines = ' -ExpandInputLines';
  end
catch
  expand_lines = '';
end

fprintf('\n Superampling in imod is [%s] with expandLines [%s]\n',super_sample ,expand_lines);
%default to cycle number zero for
%determining mean z height of particles
recWithoutMat = false;
reconstructionParameters = 0;
filterProjectionsForTomoCPRBackground=0;
loadSubTomoMeta = true;
flgWhitenPS = [0,0,0.0];
use_existing_tmpCache='';
if nargin > 2
  if isempty(EMC_str2double(varargin{1}))
    error('Extra argument to ctf 3d should be a vector [THICKNESS, BINNING] tiltN, or a string templateSearch');
  else
    reconstructionParameters = EMC_str2double(varargin{1});
    recWithoutMat = true;
    if length(varargin) > 2
      % Full recon for tomoCPR
      bh_global_turn_on_phase_plate = varargin{3};
      filterProjectionsForTomoCPRBackground = varargin{4};
      if length(varargin) > 4
        use_existing_tmpCache = varargin{5};
      end
    else
      loadSubTomoMeta = false;
      % Default to on for subregion picking
      % If user has specified phakePhasePlate, don;t use ...otherwise
      if isempty(bh_global_turn_on_phase_plate(1)) || bh_global_turn_on_phase_plate(1) == 0
        bh_global_turn_on_phase_plate = [1,2];
      end
    end
  end
elseif nargin > 1
  if strcmpi(varargin{1},'templateSearch')
    recWithoutMat = true;
    loadSubTomoMeta = false;
    if (bh_global_turn_on_phase_plate(1))
      fprintf('WARNING: the filtered tomogram should only be used for viz, not template matching.');
    end
  else
    error('Extra argument to ctf 3d should be a vector [THICKNESS, BINNING] tiltN, or a string templateSearch');
  end
else
  % Default to zero for normal use
  if isempty(bh_global_turn_on_phase_plate)
    bh_global_turn_on_phase_plate = 0;
  end
end

try
  % -1, whiten before ctf, 1 whiten after - test both.
  usr_flgWhitenPS = emc.('whitenPS');
  if (numel(usr_flgWhitenPS) == 3)
    flgWhitenPS = usr_flgWhitenPS;
  else
    error('flgWhitenPS should be a 3 element vector');
  end
catch
end

if (bh_global_turn_on_phase_plate(1) && any(emc.whitenPS))
  fprintf('WARNING: phakePhasePlate and whitening are conflicting preocesses. Turning off whitening.\n');
  emc.whitenPS = [0,0,0];
end


try
  applyExposureFilter = emc.('applyExposureFilter')
catch
  applyExposureFilter = 1;
end

% This will be set false if the reconstruction is for template matching or
% for tomoCPR
try
  useSurfaceFit = emc.('useSurfaceFit')
catch
  useSurfaceFit = 1;
end

try
  % Not for normal use, pass the total dose less first frame to flip values.
  invertDose = emc.('invertDose')
catch
  invertDose = 0;
end

%cycleNumber = sprintf('cycle%0.3d',CYCLE);
%fprintf('cycle is %d\n',CYCLE);


fprintf('tiltweight is %f %f\n',tiltWeight);



reconScaling = 1;

[tmpCache, flgCleanCache, CWD] = EMC_setup_tmp_cache(emc.fastScratchDisk, use_existing_tmpCache, 'ctf3d', false);


if (recWithoutMat)
  useSurfaceFit = false;
  if (loadSubTomoMeta)
    load(sprintf('%s.mat', emc.('subTomoMeta')), 'subTomoMeta');
    mapBackIter = subTomoMeta.currentTomoCPR;
    masterTM = subTomoMeta; clear subTomoMeta
    CYCLE = masterTM.currentCycle;
  else
    mapBackIter = 0;
    CYCLE = 0;
  end
else
  load(sprintf('%s.mat', emc.('subTomoMeta')), 'subTomoMeta');
  mapBackIter = subTomoMeta.currentTomoCPR;
  masterTM = subTomoMeta; clear subTomoMeta
  CYCLE = masterTM.currentCycle;
end

cycleNumber = sprintf('cycle%0.3d',CYCLE);
% This should be run after raw alignment and after cycle 0
if (CYCLE)
  if isfield(masterTM.(cycleNumber),'RawAlign')
    fprintf(' %s \n',cycleNumber);
  elseif isfield(masterTM.(sprintf('cycle%0.3d',CYCLE-1)),'RawAlign')
    cycleNumber = sprintf('cycle%0.3d',CYCLE-1);
    fprintf('Falling back the previous alignment cycle\n');
  else
    error('Did not find the geometry from RawAlign for %s!',cycleNumber);
  end
end


try
  flgDampenAliasedFrequencies = emc.('flgDampenAliasedFrequencies');
catch
  flgDampenAliasedFrequencies = 0;
end

try
  flg2dCTF = emc.('flg2dCTF');
catch
  flg2dCTF = 0;
end

try
  % Part of the experiment with template matching using higher res info, also
  % allow for a median filter post CTF correction, pre reconstruction to
  % further denoise prior to template matching.
  flgMedianFilter = emc.('ctfMedianFilter');
catch
  flgMedianFilter = 0;
end



%%%%% Take these from param file later.
if (reconstructionParameters(1))
  samplingRate = reconstructionParameters(2);
else
  if (loadSubTomoMeta)
    samplingRate = emc.('Ali_samplingRate');
    % This number is used to roughly balance the trade off between
    % achievable resolution, and run time during reconstruction as
    % determined by the thickness of each 3d slab reconstructed. Given that
    % we expect the resolution to improve beyond our current value, we
    % multiply by 1/2, which gives a (only loosely optimized) resTarget.
    resTarget = mean(masterTM.('currentResForDefocusError')*0.5);
    if (emc.whitenPS(1))
      emc.whitenPS(2) = resTarget;
    end
  else
    % For template search
    samplingRate = emc.('Tmp_samplingRate');
    try
      resTarget = emc.('lowResCut');
    catch
      resTarget = 12;
    end
    
  end
  
end

try
  max_ctf3dDepth = emc.('max_ctf3dDepth');
catch
  max_ctf3dDepth = 500*10^-9;
end

if (max_ctf3dDepth < 1 * 10^-9 || max_ctf3dDepth > 1000 * 10^-9)
  error('max_ctf3dDepth should be between 1 and 1000 nm');
else
  fprintf('Using a max_ctfDepth of %2.2f nm\n',max_ctf3dDepth*10^9);
end

fprintf('Using a target resolution of %2.2f Angstroms\n',resTarget);

nGPUs = emc.('nGPUs');
% Optionally specify gpu idxs
if numel(nGPUs) == 1
  gpuList = 1:nGPUs;
else
  gpuList = nGPUs;
  nGPUs = length(gpuList);
end

emc.pixel_size_angstroms = emc.pixel_size_angstroms .* samplingRate;


eraseRadius = ceil(1.5.*(emc.('beadDiameter')./emc.pixel_size_si.*0.5) / samplingRate);

nTomosPerTilt = 0;
recGeom = 0;

if (recWithoutMat)
  if reconstructionParameters(1) && loadSubTomoMeta
    tiltList{1} = varargin{2};
    nTilts = 1;
    % We just need one valid subtomot
    iTry = 1;
    tomoList{1} = '';
    while iTry < 25
      if (isfield(masterTM.mapBackGeometry.tomoName,sprintf('%s_%d',tiltList{1},iTry)))
        tomoList{1} = sprintf('%s_%d',tiltList{1},iTry);
        break;
      end
      iTry = iTry + 1;
    end
    if isempty(tomoList{1})
      error('Did not find a valid tomogram in the searchspace ->25');
    end
  else
    % TODO set up a check on the recon folder to get what is needed for
    % templateSearch
    getCoords = dir('recon/*.coords');
    nTilts = length(getCoords);
    if (nTilts == 0)
      error('Did not find any tomogram coordinates in recon/TS*.coords');
    end
    tiltList = cell(nTilts,1);
    nTomosTotal = 0;
    nTomosPerTilt = cell(nTilts,1);
    recGeom = cell(nTilts,1);
    for iStack = 1:nTilts
      [ recGeom{iStack}, tiltName, nTomosPossible] = BH_multi_recGeom( sprintf('recon/%s',getCoords(iStack).name) );
      nTomosTotal = nTomosTotal + nTomosPossible;
      nTomosPerTilt{iStack} = nTomosPossible;
      tiltList{iStack} = tiltName;
    end
    
    tomoList = cell(nTomosTotal,1);
    nTomosAdd = 0;
    for iStack = 1:nTilts
      for iTomo = 1:nTomosPerTilt{iStack}
        tomoList{nTomosAdd+1} = sprintf('%s_%d',tiltName,iTomo);
        nTomosAdd = nTomosAdd +1;
      end
    end
    
  end
else
  [tiltList,nTilts] = BH_returnIncludedTilts(masterTM.mapBackGeometry);
  tomoList = fieldnames(masterTM.mapBackGeometry.tomoName);
end


% Divide the tilt series up over each gpu
iterList = cell(nGPUs,1);
% If there is only one tilt, things break in a weird way
nGPUs = min(nGPUs, nTilts);
for iGPU = 1:nGPUs
  iterList{gpuList(iGPU)} = iGPU+(tiltStart-1):nGPUs:nTilts;
  iterList{gpuList(iGPU)};
end
try
  EMC_parpool(nGPUs)
catch
  delete(gcp('nocreate'))
  EMC_parpool(nGPUs)
end

% FIXME 
parfor iGPU = 1:nGPUs
% for iGPU = 1:nGPUs

  for iTilt = iterList{gpuList(iGPU)}
    
    iTomoList = {};
    
    % For now, since the tilt geometry is not necessarily updated (it is manual)
    % in the masterTM, check that newer (possible perTilt refined) data is
    % not present.
    TLTNAME = sprintf('fixedStacks/ctf/%s_ali%d_ctf.tlt',tiltList{iTilt},mapBackIter+1);
    TLT = load(TLTNAME);
    fprintf('iGPU %d and iTilt %d using TLT %s\n', iGPU, iTilt, TLTNAME);

    
    % Get all the tomogram names that belong to a given tilt-series.
    if (~recWithoutMat)
      nTomos = 0;
      alreadyMade = 0;
      for iTomo = 1:length(tomoList)
        if strcmp(tiltList{iTilt},masterTM.mapBackGeometry.tomoName.(tomoList{iTomo}).tiltName)
          iTomoList{nTomos+1} = tomoList{iTomo};
          nTomos = nTomos + 1;
        end
        % The order of tomo num could be off but only if all are present do we
        % skip.
        if (bh_global_turn_on_phase_plate(1))
          filtered = '_filtered';
        else
          filtered = '';
        end
        checkRecon = sprintf('cache/%s_%d_bin%d%s.rec', ...
          tiltList{iTilt},iTomo,samplingRate,filtered);
        if exist(checkRecon, 'file')
          fprintf('found %s to already exits\n',checkRecon);
          alreadyMade = alreadyMade +1;
        end
        
      end
      
      if alreadyMade == nTomos
        fprintf('All tomos 1-%d found to exist for tilt-series %s\n',nTomos,tiltList{iTilt});
        continue
      end
    end
    
    
    
    preBinStacks(TLT, tiltList{iTilt}, mapBackIter,1,...
      samplingRate,...
      PosControl2d,...
      tiltWeight,...
      flgMedianFilter);
    
  end
end

% All data is handled through disk i/o so everything unique created in the
% parfor is also destroyed there as well.
parfor iGPU = 1:nGPUs
  % for iGPU = 1:nGPUs
  gpuDevice(gpuList(iGPU));
  % Loop over each tilt
  for iTilt = iterList{gpuList(iGPU)}

    
    if (recWithoutMat)
      if (loadSubTomoMeta)
        nTomos = 1;
      else
        % templaterch
        nTomos  = nTomosPerTilt{iTilt};
        iCoords = recGeom{iTilt};
      end
    else
      nTomos = masterTM.mapBackGeometry.(tiltList{iTilt}).nTomos;
      iCoords = masterTM.mapBackGeometry.(tiltList{iTilt}).coords;
      % FIXME
    end
    
    if (recWithoutMat && ~loadSubTomoMeta) || ~recWithoutMat
      targetSizeY = diff(floor(iCoords(:,2:3)),1,2)+1;
      iCoords = iCoords ./ samplingRate;
      iCoords(:,1:4) = floor(iCoords(:,1:4));
      iCoords(:,3) = iCoords(:,3) - (diff(floor(iCoords(:,2:3)),1,2)+1 - floor(targetSizeY./samplingRate));
    end
    iTomoList = cell(nTomos,1);
    
    
    
    TLTNAME = sprintf('fixedStacks/ctf/%s_ali%d_ctf.tlt',tiltList{iTilt},mapBackIter+1);
    TLT = load(TLTNAME);
    fprintf('iGPU %d and iTilt %d using TLT %s\n', iGPU, iTilt, TLTNAME);

  
    
    if (~recWithoutMat)
      % Get all the tomogram names that belong to a given tilt-series.
      nTomos = 0;
      alreadyMade = 0;
      for iTomo = 1:length(tomoList)
        if strcmp(tiltList{iTilt},masterTM.mapBackGeometry.tomoName.(tomoList{iTomo}).tiltName)
          iTomoList{nTomos+1} = tomoList{iTomo};
          nTomos = nTomos + 1;
        end
        % The order of tomo num could be off but only if all are present do we
        % skip.
        if (bh_global_turn_on_phase_plate(1))
          filtered = '_filtered';
        else
          filtered = '';
        end
        checkRecon = sprintf('cache/%s_%d_bin%d%s.rec', ...
          tiltList{iTilt},iTomo,samplingRate,filtered);
        if exist(checkRecon, 'file')
          fprintf('found %s to already exits\n',checkRecon);
          alreadyMade = alreadyMade +1;
        end
        
      end
      
      if alreadyMade == nTomos
        fprintf('All tomos 1-%d found to exist for tilt-series %s\n',nTomos,tiltList{iTilt});
        continue
      end
    end
    
    
    
    if samplingRate > 1
      fullStack = sprintf('%aliStacks/%s_ali%d.fixed', ...
        tiltList{iTilt},mapBackIter+1);
      inputStack = sprintf('cache/%s_ali%d_bin%d.fixed',...
        tiltList{iTilt},mapBackIter+1,samplingRate);
      if ~exist(inputStack, 'file')
        %     binCMD = sprintf('newstack -bin %d -antialias 6 %s %s > /dev/null',samplingRate,fullStack,inputStack);
        % %     binCMD = sprintf('newstack -bin %d -antialias 6 %s %s ',samplingRate,fullStack,inputStack);
        %
        %     system(binCMD);
        BH_multi_loadOrBin(fullStack, samplingRate, 2, false);
        
      end
    else
      inputStack = sprintf('aliStacks/%s_ali%d.fixed',...
        tiltList{iTilt},mapBackIter+1);
    end
    
    %       system(sprintf('header %s',inputStack));
    % iHeader = MRCImage(inputStack,0);
    % STACK = gpuArray(single(getVolume(iHeader)));
    
    maskedStack = single(getVolume(MRCImage(inputStack)));
    
    if (recWithoutMat)
      if (reconstructionParameters(1) && loadSubTomoMeta)
        NX = size(maskedStack,1);
        NY = size(maskedStack,2);
        
        %           NY = size(maskedStack,2)-1;
        NZ = floor(reconstructionParameters(1))
        maxZ = NZ;
        iCoords = [NX,0,NY-1,NZ,0,0];
        tomoNumber = 1;
      else
        [ ~, maxZ, tomoNumber, ~ ] = calcAvgZ('dummy',iCoords,tiltList{iTilt}, ...
          iTomoList,nTomos, emc.pixel_size_angstroms, ...
          samplingRate, cycleNumber,...
          0,1);
      end
    else
      [ ~, maxZ, tomoNumber, ~ ] = calcAvgZ(masterTM,iCoords,tiltList{iTilt}, ...
        iTomoList,nTomos, emc.pixel_size_angstroms, ...
        samplingRate, cycleNumber,...
        0,1);
    end
    
    if ( flg2dCTF || recWithoutMat && loadSubTomoMeta)
      nSections = 1;
      ctf3dDepth = maxZ * 10 ^ -9;
    else
      dampeningMax = 0.90;
      [ ctf3dDepth ] = BH_ctfCalcError( samplingRate*mean(TLT(:,16)), ...
        TLT(1,17),TLT(1,18),abs(TLT(1,15)), ...
        2048, TLT(1,19), ...
        resTarget,maxZ*10, ...
        dampeningMax,CYCLE);
      fprintf('\n\nCalculated a ctfDepth of %2.2f nm for %s\n\n',ctf3dDepth*10^9,tiltList{iTilt});
      if (ctf3dDepth > max_ctf3dDepth)
        ctf3dDepth = max_ctf3dDepth;
        fprintf('Calculated ctfDepth exceeds user specified max, so actually using a max_ctfDepth of %2.2f nm\n',ctf3dDepth*10^9);
      end
      % sections centered at 0, which for now is also supposed to coincide with
      % the mean defocus determination, although this could be corrected using
      % knowledge of particle positions given assurance that particles are the
      % primary source of signal (and not carbon for example).
      nSections = ceil(maxZ/(ctf3dDepth*10^9));
      % max odd number
      nSections = nSections + ~mod(nSections,2);
    end
    fprintf('with %3.3f nm sections, correcting %d tilt-series\n', ctf3dDepth*10^9, nSections);
    
    % For each tomo create a list of slices that are to be reconstructed
    % for every section section.
    
    [ sectionList ] = calcTomoSections(iCoords, tomoNumber, emc.pixel_size_angstroms, nSections, tiltList{iTilt}, ctf3dDepth);
    

    if (recWithoutMat)
      avgZ = 0;
      surfaceFit = 0;
    else
      [ avgZ, maxZ, tomoNumber, surfaceFit ] = calcAvgZ(masterTM,iCoords,tiltList{iTilt}, ...
                                                        iTomoList,nTomos, emc.pixel_size_angstroms, ...
                                                        samplingRate, cycleNumber,...
                                                        sectionList, 0);
    end
    
    if ( shiftDefocusOrigin )
      fprintf('Using avgZ %3.3e nm as the defocus origin\n',avgZ*10^9);
    else
      avgZ = 0;
      fprintf('Using sample origin as the defocus origin\n');
      fprintf('If you want to use the COM of subTomos, set shiftDefocusToSubTomoCOM=1\n');
    end
    
    
    
    for iSection = 1:nSections
      
      defFitFull = '';
      preCombDefocus = 0;
      if (mapBackIter)
        defFitFull = sprintf('mapBack%d/%s_ali%d_ctf.defFidFull',mapBackIter, tiltList{iTilt},mapBackIter);
        if exist(defFitFull,'file')
          preCombDefocus = load(defFitFull);
          fprintf('3dCTF using pre calc combined per tilt defocus %s\n',defFitFull);
        end
      end
      
      if (PosControl2d)
        correctedStack = maskedStack;
      else
        
        % I would have thought the global would be recognized, but it looks
        % like there is something odd about its use with a parfor loop
        % FIXME, when setting up the iterator, make clean copies for each
        % worker that are local in scope.e
        
        
        [ correctedStack ] = ctfMultiply_tilt(nSections,iSection,ctf3dDepth, ...
          avgZ,TLT,emc.pixel_size_angstroms,maskedStack,...
          maxZ*10/emc.pixel_size_angstroms,flgDampenAliasedFrequencies,...
          preCombDefocus,samplingRate,...
          applyExposureFilter,surfaceFit,...
          useSurfaceFit,invertDose,...
          bh_global_turn_on_phase_plate,...
          filterProjectionsForTomoCPRBackground,...
          emc.whitenPS);
      end
      % Write out the stack to the cache directory as a tmp file
      
      if (flgEraseBeads_aferCTF)
        scalePixelsBy = samplingRate;
        correctedStack = BH_eraseBeads(correctedStack,eraseRadius, tiltList{iTilt}, scalePixelsBy,mapBackIter,sortrows(TLT,1));
      end
      
      
      outputStack = sprintf('%s/%s_ali%d_%d.fixed', tmpCache,tiltList{iTilt},mapBackIter+1,iSection);
      SAVE_IMG(correctedStack,outputStack,emc.pixel_size_angstroms);
      correctedStack = [];
      
      % Loop over tomos reconstructing section and appending a file to
      for iT = 1:nTomos
        
        thisTomo = tomoNumber(iT);
        
        if any(sectionList{iT}(iSection,:)+9999)
          
          
          
          reconName = sprintf('%s/%s_ali%d_%d_%d.rec', ...
            tmpCache,tiltList{iTilt},mapBackIter+1,thisTomo,iSection);
          
          
          if (loadSubTomoMeta)
            if (recWithoutMat)
              TA = sortrows(masterTM.tiltGeometry.(tomoList{1}),1);
            else
              TA = sortrows(masterTM.tiltGeometry.(sprintf('%s_%d',tiltList{iTilt},thisTomo)),1);
            end
            TA = TA(:,4);
          else
            if (mapBackIter)
              % FIXME: I don't think this block should work, it should only be the tilt angles!
              error('THis block should not be reached.')
              TA = load(sprintf('%smapBack%d/%s_ali%d_ctf.tlt',CWD,mapBackIter,tiltList{iTilt},...
                mapBackIter));
            else
              TA = load(sprintf('%sfixedStacks/%s.tlt',CWD,tiltList{iTilt}));
            end
          end
          
          rawTLT = sprintf('cache/%s_%d.rawtlt',tiltList{iTilt},thisTomo);
          rawTLT_file = fopen(rawTLT, 'w');
          fprintf(rawTLT_file,'%f\n', TA');
          fclose(rawTLT_file);
          
          if (mapBackIter)
            
            LOCAL = sprintf('%smapBack%d/%s_ali%d_ctf.local',CWD,mapBackIter,tiltList{iTilt}, ...
              mapBackIter);
          else
            LOCAL = sprintf('%sfixedStacks/%s.local',CWD,tiltList{iTilt});
          end
          
          
          
          if exist(LOCAL,'file')
            flgLocal = 1;
          else
            fprintf('Did not find local alignment information at %s\n',LOCAL);
            flgLocal = 0;
          end
          
          % hangover from slab padding, remove later.
          padRec = 0;
          
          nTiltWorkers = 2;
          nTotalSlices = (iCoords(thisTomo,3)-iCoords(thisTomo,2)+1);
          tiltChunkSize = ceil(nTotalSlices/nTiltWorkers);
          tiltChunks = iCoords(thisTomo,2):tiltChunkSize:iCoords(thisTomo,3);
          tiltChunks(end) = iCoords(thisTomo,3);
          totalSlices = [tiltChunks(1),tiltChunks(end)];
          
          
          
          rCMD = sprintf(['tilt %s %s -input %s -output %s.TMPPAD -TILTFILE %s -UseGPU %d ', ...
            '-WIDTH %d -COSINTERP 0 -THICKNESS %d -SHIFT %f,%f '],...
            super_sample, expand_lines, ...
            outputStack, reconName, rawTLT, gpuList(iGPU), ...
            iCoords(thisTomo,1),floor(sectionList{iT}(iSection,5))+2*padRec,...
            iCoords(thisTomo,5),sectionList{iT}(iSection,6));
          
          
          
          % Explicitly set Radial to Nyquist
          if (flgLocal)
            rCMD = [rCMD sprintf('-LOCALFILE %s -RADIAL 0.5,.05 -MODE 2 -SCALE 0,%d',LOCAL,reconScaling)];
          else
            rCMD = [rCMD sprintf('-RADIAL 0.5,.05 -MODE 2 -SCALE 0,%d',reconScaling)];
          end
          
          system(sprintf('rm -f %s.sh',reconName));
          
          recScript = fopen(sprintf('%s.sh',reconName),'w');
          fprintf(recScript,'#!/bin/bash\n\n');
          fprintf(recScript,'%s -SLICE -1,-1 -TOTALSLICES %d,%d\n',rCMD,totalSlices);
          for iRecSec = 1:nTiltWorkers-1
            if iRecSec < nTiltWorkers -1
              iShift = 1;
            else
              iShift = 0;
            end % /dev/null
            fprintf(recScript,'%s -SLICE %d,%d -TOTALSLICES %d,%d > /dev/null  &\n',rCMD, ...
              tiltChunks(iRecSec),...
              tiltChunks(iRecSec+1)-iShift,...
              totalSlices);
          end
          fprintf(recScript,'\n\nwait\n\n');
          fclose(recScript);
          system(sprintf('chmod a=wrx %s.sh',reconName));
          
          [recError,~] = system(sprintf('%s.sh > /dev/null ',reconName)); % /dev/null
          if (recError)
            system(sprintf('%s.sh',reconName));
            error('\n\nerror during reconstruction %s\n\n',reconName);
          end
          % Z coords (y in this orientation) are decreasing into the
          % monitor. For symmetrical padding this doesn't matter, but keep
          % in mind. /dev/null
          trimCMD = sprintf('trimvol -rx -y %d,%d %s.TMPPAD %s > /dev/null  ' , ...
            padRec+1,floor(sectionList{iT}(iSection,5))+padRec,reconName,reconName);
          % % %           trimCMD = sprintf('newstack -fromone -secs %d-%d %s.TMPPAD %s > /dev/null', ...
          % % %                              padRec+1,floor(sectionList{iT}(iSection,5))+padRec,reconName,reconName)
          [msg,~]= system(trimCMD);
          if (msg)
            fprintf('%d from trimCMD\n',msg)
            trimCMDPrintError = sprintf('trimvol -rx -y %d,%d %s.TMPPAD %s', ...
              padRec+1,floor(sectionList{iT}(iSection,5))+padRec,reconName,reconName)
            % % %             trimCMDPrintError = sprintf('newstack -fromone -secs %d-%d %s.TMPPAD %s', ...
            % % %                              padRec+1,floor(sectionList{iT}(iSection,5))+padRec,reconName,reconName)
            system(trimCMDPrintError);
          end
          system(sprintf('rm %s.TMPPAD', reconName));
          %  fprintf([trimCMD ' \n'])
          
          
        end
        
      end % end loop over tomos for this section
      
      system(sprintf('rm %s',outputStack));
      
    end % end loop over sections
    
    deltaZ = [];
    evalMask = [];
    maskedStack = [];
    
    for iT = 1:nTomos
      thisTomo = tomoNumber(iT);
      
      if (bh_global_turn_on_phase_plate(1))
        reconNameFull = sprintf('cache/%s_%d_bin%d_filtered.rec', ...
          tiltList{iTilt},thisTomo,samplingRate);
      elseif reconstructionParameters(1)
        reconNameFull = sprintf('%scache/%s_%d_bin%d_backgroundEst.rec', ...
          CWD,tiltList{iTilt},thisTomo,samplingRate);
      else
        reconNameFull = sprintf('cache/%s_%d_bin%d.rec', ...
          tiltList{iTilt},thisTomo,samplingRate);
      end
      fprintf('in ctf3d reconNameFull is %s\n\n',reconNameFull);
      
      recCMD = 'newstack -fromone';
      for iSection = 1:nSections
        reconName = sprintf('%s/%s_ali%d_%d_%d.rec', ...
          tmpCache, tiltList{iTilt},mapBackIter+1,thisTomo,iSection);
        
        if any(sectionList{iT}(iSection,:)+9999)
          recCMD = [recCMD,sprintf(' -secs 1-%d %s', ...
            floor(sectionList{iT}(iSection,5)), ...
            reconName)];
        else
          
          fprintf('no info for section %d for tomo %d\n',iSection,thisTomo);
        end
      end
      
      system([recCMD, sprintf(' %s > /dev/null ',reconNameFull)]); %/dev/null
      
      
      for iSection = 1:nSections
        cleanUp3 = sprintf('rm %s/%s_ali%d_%d_%d.rec', ...
          tmpCache,tiltList{iTilt},mapBackIter+1,thisTomo,iSection);
        system(cleanUp3);
        cleanUp4 = sprintf('rm %s/%s_ali%d_%d_%d.rec.sh', ...
          tmpCache,tiltList{iTilt},mapBackIter+1,thisTomo,iSection);
        system(cleanUp4);
        
      end
      
      
    end % end of recombination loop
    
    maskedStack = [];
    
  end % end of loop over tilt-series
end % end of parfor over gpus

if (flgCleanCache)
  % Double check that this exists to avoid data loss.
  checkDir = dir(tmpCache);
  if isempty(checkDir)
    fprintf('not removing the temp cache because it did not eval with dir\n');
  else
    cleanItUp = sprintf('rm  %s/*',tmpCache);
    system(cleanItUp);
  end
end

end

function [] = preBinStacks(TLT, STACK_PRFX, mapBackIter,usableArea,...
  samplingRate,...
  PosControl2d,...
  tiltWeight,...
  flgMedianFilter)



if (PosControl2d)
  prefix = 'ctf';
  suffix = '_ctf';
else
  prefix = 'ali';
  suffix = '';
end


fullStack = sprintf('%sStacks/%s_ali%d%s.fixed', ...
  prefix,STACK_PRFX,mapBackIter+1,suffix);
inputStack = sprintf('cache/%s_ali%d%s_bin%d.fixed',...
  STACK_PRFX,mapBackIter+1,suffix,samplingRate);
if ~exist(inputStack, 'file')
  BH_multi_loadOrBin(fullStack, samplingRate, 2, false);
end


end


function  [ sectionList ] = calcTomoSections(iCoords, tomoNumber, pixel_size_angstroms, nSections,tiltName, ctf3Depth)

nTomos = length(tomoNumber);
sectionList = cell(nTomos,1);
for iTomo = 1:nTomos
  % min and max in absolute pixels min and max from 1:nZrecon
  sectionList{iTomo} = zeros(nSections,6);
end

% With rounding this could end up a bit short except the top and bottom are both
% half a section larger than minimally needed.
nSec = floor(ctf3Depth*10^10/pixel_size_angstroms)      ;
nSec = nSec + ~mod(nSec,2);
halfSec = (nSec-1)/2;

for iT = 1:length(tomoNumber)
  iTomo = tomoNumber(iT);
  % Origin + originshift
  
  -1.*(ceil((iCoords(iTomo,4)+1)/2)-1) + iCoords(iTomo,6),
  reconRange = floor([-1.*(ceil((iCoords(iTomo,4)+1)/2)-1) + iCoords(iTomo,6),0]);
  reconRange(2) = reconRange(1) + iCoords(iTomo,4) - 1;
  nZ = 1;
  flgFirstSec = 1;
  
    for iSection = 1:nSections
    
    sectionCenter = ((nSections-1)/-2+(iSection-1))*(nSec-1);
    
    % Check that sectionCenter is within range
    if sectionCenter + halfSec < reconRange(1) || ...
        sectionCenter - halfSec > reconRange(2)
      sectionList{iT}(iSection,:) = -9999;
    else
      
      if (sectionCenter - halfSec > 0)
        sectionList{iT}(iSection,1) = max(sectionCenter - halfSec ,reconRange(1));
        if sectionList{iT}(iSection,1) ~= reconRange(1)
          sectionList{iT}(iSection,1) = sectionList{iT}(iSection,1) +1;
        end
      elseif (sectionCenter - halfSec < 0)
        sectionList{iT}(iSection,1) = max(sectionCenter - halfSec,reconRange(1));
        if sectionList{iT}(iSection,1) ~= reconRange(1)
          sectionList{iT}(iSection,1) = sectionList{iT}(iSection,1) +1;
        end
      else
        sectionList{iT}(iSection,1) = max(-halfSec,reconRange(1));
      end
      
      if (sectionCenter - halfSec > 0)
        sectionList{iT}(iSection,2) = min(sectionCenter + halfSec,reconRange(2));
      elseif (sectionCenter - halfSec < 0)
        sectionList{iT}(iSection,2) = min(sectionCenter + halfSec ,reconRange(2));
      else
        sectionList{iT}(iSection,2) = min(halfSec,reconRange(2));
      end
      
      
      % Check that first section starts in the correct place. If a small error,
      % just shift the results, otherwise complain.
      if (flgFirstSec)
        if sectionList{iT}(iSection,1) ~= reconRange(1)
          if abs(sectionList{iT}(iSection,1) - reconRange(1)) < 10
            sectionList{iT}(iSection,1) = reconRange(1);
            fprintf('\n\nShifting first section %s_n%d\n\n',tiltName,iTomo);
          else
            error('section start %d is too far off from expected %d\n', ...
              sectionList{iT}(iSection,1), reconRange(1))
          end
        end
        flgFirstSec = 0;
      end
      
      secZ = sectionList{iT}(iSection,2) -  sectionList{iT}(iSection,1);
      sectionList{iT}(iSection,3) = nZ;
      sectionList{iT}(iSection,4) = nZ + secZ;
      sectionList{iT}(iSection,5) = secZ +1;
      sectionList{iT}(iSection,6) = (secZ+1)./2 + sectionList{iT}(iSection,1);
      nZ = nZ + secZ + 1;
    end
    
    
    % not a good solution, but not sure just yet why I'm getting some occasionally
    % weird results.
    if sectionList{iT}(iSection,5) < 3
      sectionList{iT}(iSection,:) = -9999;
    end
  end % end loop over sections
  
  % TroubleShoot
  tSHT = fopen(sprintf('.tblSht_%s_i%d.txt',tiltName,iTomo),'w');
  fprintf(tSHT,'%2.2f %2.2f %2.2f %2.2f %2.2f %2.2f\n', iCoords(iTomo,:)');
  fprintf(tSHT,'%2.2f %2.2f %2.2f %2.2f %2.2f %2.2f\n', sectionList{iT}');
  fclose(tSHT);
end % end loop over tomos



end




function [correctedStack] = ctfMultiply_tilt(nSections,iSection,ctf3dDepth, ...
  avgZ,TLT,pixel_size_angstroms,maskedStack,...
  maxZ,flgDampenAliasedFrequencies,...
  preCombDefocus,samplingRate,...
  applyExposureFilter,surfaceFit,...
  useSurfaceFit,invertDose, ...
  phakePhasePlate, ...
  filterProjectionsForTomoCPRBackground,...
  flgWhitenPS)
% Correct in strips which is more expensive but (hopefully) more accurate.



% For sections with too few subTomos to fit, fall back

if isa(surfaceFit,'cell')
  surfaceFit = surfaceFit{iSection};
  if ~isa(surfaceFit,'sfit')
    fprintf('Warning, surfaceFit is not an sfit object\n');
    useSurfaceFit = false;
  end
else
  useSurfaceFit = false;
end

[d1,d2,nPrjs] = size(maskedStack);

% This is just going to be written out to disk so keep in main memory.
correctedStack = zeros(d1,d2,nPrjs,'single');

% This should probably come from the calcSections. If avgZ has not been
% forced to zero, then assume the defocus determined is the distance from
% the center of mass of subtomograms in Z to the focal plane, rather than
% the center of mass of the tomograms (specimen)

if (useSurfaceFit)
  defocusOffset = 0;
else
  defocusOffset = (((nSections-1)/-2+(iSection-1))*ctf3dDepth);
  fprintf('Not using surface fit, so using offset %3.3e for section %d with COM offset %3.3e\n', defocusOffset, iSection, avgZ);
  % Assuming the majority of the fit defocus came from the subtomograms, then the estimated defocus value needs to be moved from
  % the origin of the specimen to the origin of the subtomograms.
  defocusOffset = (defocusOffset + avgZ); % The average height of the particles is factored into the surface fit
end



if ( flgDampenAliasedFrequencies )
  % Experiment with dampning higher frequencies where aliasing is going to
  % result in nonsens.
  flgDampenAlias = 1;
  fprintf('\n\nExperimental dampening of aliased CTF terms\n');
else
  flgDampenAlias = 0;
end

apoSize = 6;
fastFTSize = BH_multi_iterator([d1,d2],'fourier2d');


% These should be constant for a given tiltseries
Cs = TLT(1,17);
WAVELENGTH = TLT(1,18);
AMPCONT = TLT(1,19);

%if d2 < padTileSize
%  padYdim = padTileSize;
%else
%  padYdim = d2;
%end
if ( flgDampenAlias )
  % Calculate a centered grid b/c real space convolution
  [radialGrid,phi,~,~,~,~] = BH_multi_gridCoordinates(fastFTSize, ...
    'Cylindrical','GPU', ...
    {'none'},1,1,0);
else
  [radialGrid,phi,~,~,~,~] = BH_multi_gridCoordinates(fastFTSize, ...
    'Cylindrical','GPU', ...
    {'none'},1,0,0);
end
radialGrid = {radialGrid./(pixel_size_angstroms*10^-10),0,phi};
phi = [];

fprintf('%f %f\n',filterProjectionsForTomoCPRBackground,pixel_size_angstroms);
if (filterProjectionsForTomoCPRBackground ~= 0)
  bpFilter = BH_bandpass3d(fastFTSize,0, 0, filterProjectionsForTomoCPRBackground, 'GPU',pixel_size_angstroms);
else
  bpFilter = 1;
end

for iPrj = 1:nPrjs
  
  maxEval = cosd(TLT(iPrj,4)).*(d1/2) + maxZ./2*abs(sind(TLT(iPrj,4)));
  oX = emc_get_origin_index(d1);
  oY = emc_get_origin_index(d2);
  iEvalMask = floor(oX-maxEval):ceil(oX+maxEval);
  
  if ( applyExposureFilter )
    iExposureFilter = BH_exposureFilter(fastFTSize,TLT(iPrj,:),'GPU',samplingRate,0);
  else
    iExposureFilter = 1;
  end
  
  iExposureFilter = iExposureFilter .* bpFilter;

  
  ddF = TLT(iPrj,12);
  dPhi = TLT(iPrj,13);
  D0 = abs(TLT(iPrj,15));
  %TLT(iPrj,16);
  
  
  padVal = BH_multi_padVal([d1,d2],fastFTSize);
  trimVal = BH_multi_padVal(fastFTSize,[d1,d2]);
  
  
  iProjection = BH_padZeros3d(maskedStack(:,:,TLT(iPrj,1)),padVal(1,:),padVal(2,:),'GPU','singleTaper');
  iProjectionFT = fftn(iProjection).*iExposureFilter; clear iExposureFilter
  correctedPrj = zeros([d1,d2],'single','gpuArray');
  
  % Gridvectors for the specimen plane
  [rX,rY,~] = BH_multi_gridCoordinates([d1,d2],'Cartesian','GPU',{'none'},0,1,0);
  % Assuming the plane fit is from the origin as it is.
  
  if (useSurfaceFit)
    %rZ = (surfaceFit.p00 + surfaceFit.p10.*(rX+oX)) + surfaceFit.p01.*(rY+oY);
    try
      rZ = feval(surfaceFit,rX,rY);
    catch
      d1
      d2
      rX
      rY
      surfaceFit
      error('surface fit failed');
    end
  else
    rZ = zeros([d1,d2],'single','gpuArray');
  end

  
  defocus_adj = D0 - (defocusOffset.*cosd(TLT(iPrj,4)));
  
  % For a positive angle, this will rotate the positive X axis farther from the focal plane (more underfocus)
  rA = BH_defineMatrix([0,TLT(iPrj,4),0],'SPIDER','inv');

  % Transform the specimen plane
  tX = round(rA(1).*rX + rA(4).*rY + rA(7).*rZ + oX);
  tY = round(rA(2).*rX + rA(5).*rY + rA(8).*rZ + oY);
  % undefocus is positive, but we have stored the negative value (so we just add this to the positional offset)
  tZ = defocus_adj - (pixel_size_angstroms*10^-10).*(rA(3).*rX + rA(6).*rY + rA(9).*rZ);

  % Some edge pixels can be out of bounds depending on the orientation of
  % the plan fit. Setting to zero will will ignore them (assuming defocus
  % is always < 0)
  tZ( tX < 1 | tY < 1 | tX > d1 | tY > d2) = 1;
  
  
  minDefocus = min(tZ(:));
  maxDefocus = max(tZ(tZ < 1));
  % Spit out some info
  
  % To track sampling in case I put in overlap
  samplingMask = zeros([d1,d2],'single','gpuArray');
  
  % TODO: confirm the astigmatism is correct, check - and +/- 90
  for iDefocus = minDefocus-ctf3dDepth/1:ctf3dDepth/1:maxDefocus+ctf3dDepth/1
    defVect = [iDefocus + ddF, iDefocus - ddF, dPhi];
    
    if (phakePhasePlate(1) > 0)
      
      if numel(phakePhasePlate) == 2
        modPower = floor(phakePhasePlate(2));
        SNR = rem(phakePhasePlate(2),1);
      else
        modPower = 1;
        SNR = 1;
      end
      
      
      [Hqz, ~] = BH_ctfCalc(radialGrid,Cs,WAVELENGTH,defVect,fastFTSize,AMPCONT,-1,1,SNR);
      
      Hqz = (-1).^modPower.*(phakePhasePlate(1).*Hqz).^1;
      
      
      modHqz = [];
    else
      if (pixel_size_angstroms < 2.0)
        % use double precision - this is not enabled, but needs to be -
        % requires changes to radial grid as well.
        Hqz = BH_ctfCalc(radialGrid,Cs,WAVELENGTH,defVect,fastFTSize,AMPCONT,-1,-1);
      else
        Hqz = BH_ctfCalc(radialGrid,Cs,WAVELENGTH,defVect,fastFTSize,AMPCONT,-1);
      end
    end
    
 
    if (flgWhitenPS(3))
      tmpCorrection = BH_padZeros3d(real(ifftn(iProjectionFT.*Hqz./(abs(Hqz).^2+flgWhitenPS(3)))),trimVal(1,:),trimVal(2,:),'GPU','single');
    else
      tmpCorrection = BH_padZeros3d(real(ifftn(iProjectionFT.*Hqz)),trimVal(1,:),trimVal(2,:),'GPU','single');
    end
    
    tmpMask = (tZ > iDefocus - ctf3dDepth/2 & tZ <= iDefocus + ctf3dDepth/2);
    
    linearIDX =  unique(sub2ind([d1,d2],tX(tmpMask),tY(tmpMask)));
    
    correctedPrj(linearIDX) = correctedPrj(linearIDX) + tmpCorrection(linearIDX);
    samplingMask(linearIDX) = samplingMask(linearIDX) + 1;
    
  end % end loop over defocus values
  
  samplingMask(samplingMask == 0) = 1;
  
  if (flgWhitenPS(1))
    correctedStack(:,:,TLT(iPrj,1)) =gather(BH_whitenNoiseSpectrum(correctedPrj./samplingMask,'',pixel_size_angstroms,1));
  else
    
    correctedStack(:,:,TLT(iPrj,1)) = gather(correctedPrj./samplingMask);
  end
  clear correctedPrj samplingMask tmpMask tmpCorrection
  
  clear iProjection iProjectionFT
end % end loop over projections
clear tile Hqz
end

function [avgZ, maxZ, tomoNumber,surfaceFit] = calcAvgZ(masterTM,iCoords, ...
                                                        tiltName,tomoList,...
                                                        nTomos, pixel_size_angstroms,...
                                                        samplingRate,cycleNumber,...
                                                        sectionList,calcMaxZ)

% Calculate the maximum extensions in Z and then how many separate sections
% need to be corrected.

surfaceFit = '';
avgZ = 0;


if isa(masterTM,'struct')
  val_to_pass =  masterTM.mapBackGeometry.tomoName;
else
  val_to_pass = 'dummy';
end

[ maxZ, tomoNumber ] = emc_get_max_specimen_NZ(val_to_pass, ...
                                              iCoords * samplingRate, ...
                                              tomoList, ...
                                              nTomos, ...
                                              samplingRate);

maxZ = maxZ .* pixel_size_angstroms ./ 10;
fprintf('combining the thickness and shift on tilt %s, found a maxZ  %3.3f nm\n',tiltName,maxZ);

if (calcMaxZ)
  return;
end

% For now use cycle000, if adding a refinment focused on a specific set of
% particles, then consider that later.
try
  initGeom = masterTM.(cycleNumber).RawAlign;
  fprintf('Loaded the geometry for RawAlign %s\n',cycleNumber);
catch
  fprintf('Failed to load the geometry for RawAlign %s\nTrying cycle000\n',cycleNumber);
  try
    initGeom = masterTM.cycle000.geometry;
  catch
    error(['Could not load the initial geometry subTomoMeta.%s.geometry\n or--',...
      'subTomoMeta.cycle000.geometry\n'],cycleNumber);
  end
end

nSubTomos = 0;
totalZ = 0;
% for each tomogram get the size and origin in Z then find mean subTomo
% position.

nSections = size(sectionList{1},1);
xFull = cell(nSections,1);
yFull = cell(nSections,1);
zFull = cell(nSections,1);
surfaceFit = cell(nSections,1);

% Initialize with empty arrays
for iSection = 1:nSections
  xFull{iSection} = [];
  yFull{iSection} = [];
  zFull{iSection} = [];
end


for iT = 1:nTomos

  iTomo = tomoNumber(iT);

  % X in the Y frame means a vector from the Y lower left to the X origin
  % X origin wrt Y origin is a vector from the origin of Y to the X origin
  tomoReconCoords = masterTM.reconGeometry.(tomoList{iT}) ./ samplingRate;
  tomo_origin_in_tomo_frame = emc_get_origin_index(tomoReconCoords(1,1:3));
  tomo_origin_wrt_specimen_origin = -1.*tomoReconCoords(2,1:3);
 
  iTomoName = sprintf('%s_%d',tiltName,iTomo);
  
  % Get the z-coordinates of the origin for all included subtomograms relative to the lower left of the tomogram
  % shouldn't be any removed particles at this stage but later there would be.
  subtomo_origin_z_in_tomo_frame = initGeom.(iTomoName)(initGeom.(iTomoName)(:,26)~=-9999,13) ./ samplingRate;

  % We should not get to this point if all subtomograms have been removed.
  if isempty(subtomo_origin_z_in_tomo_frame)
    error('No subtomograms found for %s',iTomoName);
  end

  % shift from lower left to centered and include the tomos offset from the
  subtomo_origin_wrt_specimen_origin = subtomo_origin_z_in_tomo_frame - tomo_origin_in_tomo_frame(3) + tomo_origin_wrt_specimen_origin(3);

  totalZ = totalZ + sum(subtomo_origin_wrt_specimen_origin);
  fprintf('%s tomo has %d subTomos with mean Z %3.3f nm\n', ...
    iTomoName, length(subtomo_origin_wrt_specimen_origin), ...
    mean(subtomo_origin_wrt_specimen_origin) * pixel_size_angstroms ./ 10);

  nSubTomos = nSubTomos + length(subtomo_origin_wrt_specimen_origin);

  for iSection = 1:nSections
    
    iSecOrigin = sectionList{iT}(iSection,6);
    iSecRadius = sectionList{iT}(iSection,5)/2;
    inSectionIDX = subtomo_origin_wrt_specimen_origin >  iSecOrigin - iSecRadius & subtomo_origin_wrt_specimen_origin <= iSecOrigin + iSecRadius;
    
    
    x = initGeom.(iTomoName)(initGeom.(iTomoName)(:,26)~=-9999,11)./samplingRate;
    x = x - tomo_origin_in_tomo_frame(1) + tomo_origin_wrt_specimen_origin(1);
    y = initGeom.(iTomoName)(initGeom.(iTomoName)(:,26)~=-9999,12)./samplingRate;
    y = y - tomo_origin_in_tomo_frame(2) + tomo_origin_wrt_specimen_origin(2);
    
    
    xFull{iSection} = [xFull{iSection} ; x(inSectionIDX)];
    yFull{iSection} = [yFull{iSection} ; y(inSectionIDX)];
    zFull{iSection} = [zFull{iSection} ; subtomo_origin_wrt_specimen_origin(inSectionIDX)];
    
  end % loop over sections
  clear subtomo_origin_wrt_specimen_origin
end % loop over tomos


avgZ = totalZ / nSubTomos*pixel_size_angstroms / 10*10^-9;


for iSection = 1:nSections
  
  if length(xFull{iSection}) >= 6
    surfaceFit{iSection} = fit([xFull{iSection}, yFull{iSection}],zFull{iSection},'poly22','Robust','on');
  else
    surfaceFit{iSection} = 0;
  end
end

fprintf('%s tilt-series has %d subTomos with mean Z %3.3f nm\n', tiltName, nSubTomos,avgZ*10^9);

end
