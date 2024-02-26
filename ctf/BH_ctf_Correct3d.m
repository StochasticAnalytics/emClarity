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
subTomoMeta = struct();
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
reconstructionParameters = 0;
filterProjectionsForTomoCPRBackground=0;
flgWhitenPS = [0,0,0.0];
use_existing_tmpCache='';
recon_for_tomoCPR = false;
recon_for_templateMatching = false;
recon_for_subTomo = false;
if nargin > 2
  if isempty(EMC_str2double(varargin{1}))
    error('Extra argument to ctf 3d should be a vector [THICKNESS, BINNING] tiltN, or a string templateSearch');
  else
    reconstructionParameters = EMC_str2double(varargin{1});
    if length(varargin) > 2
      recon_for_tomoCPR = true;
      % Full recon for tomoCPR
      bh_global_turn_on_phase_plate = varargin{3};
      filterProjectionsForTomoCPRBackground = varargin{4};
      if length(varargin) > 4
        use_existing_tmpCache = varargin{5};
      end
    else
      error('This block should ont be reached.');
    end
  end
elseif nargin > 1
  if strcmpi(varargin{1},'templateSearch')
    recon_for_templateMatching = true;
    if (bh_global_turn_on_phase_plate(1))
      fprintf('WARNING: the filtered tomogram should only be used for viz, not template matching.');
    end
  else
    error('Extra argument to ctf 3d should be a vector [THICKNESS, BINNING] tiltN, or a string templateSearch');
  end
else
  % Default to zero for normal use
  recon_for_subTomo = true;
  if isempty(bh_global_turn_on_phase_plate)
    bh_global_turn_on_phase_plate = 0;
  end
end

if (recon_for_tomoCPR + recon_for_templateMatching + recon_for_subTomo ~= 1)
  error('Only one of the three modes can be used at a time');
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

[tmpCache, flgCleanCache, CWD] = EMC_setup_tmp_cache(emc.fastScratchDisk, use_existing_tmpCache, 'ctf3d', false);


if (recon_for_tomoCPR || recon_for_templateMatching)
  useSurfaceFit = false;
end

if (recon_for_templateMatching)
  mapBackIter = 0;
  CYCLE = 0;
else
  load(sprintf('%s.mat', emc.('subTomoMeta')), 'subTomoMeta');
  mapBackIter = subTomoMeta.currentTomoCPR;
  CYCLE = subTomoMeta.currentCycle;
end

cycleNumber = sprintf('cycle%0.3d',CYCLE);
% This should be run after raw alignment and after cycle 0
if (CYCLE)
  if isfield(subTomoMeta.(cycleNumber),'RawAlign')
    fprintf(' %s \n',cycleNumber);
  elseif isfield(subTomoMeta.(sprintf('cycle%0.3d',CYCLE-1)),'RawAlign')
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
if (recon_for_tomoCPR)
  samplingRate = reconstructionParameters(2);
else
  if (recon_for_subTomo)
    samplingRate = emc.('Ali_samplingRate');
    % This number is used to roughly balance the trade off between
    % achievable resolution, and run time during reconstruction as
    % determined by the thickness of each 3d slab reconstructed. Given that
    % we expect the resolution to improve beyond our current value, we
    % multiply by 1/2, which gives a (only loosely optimized) resTarget.
    resTarget = mean(subTomoMeta.('currentResForDefocusError')*0.5);
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
tiltRecGeom = 0;
if (recon_for_subTomo)
  [tiltList,nTilts] = BH_returnIncludedTilts(subTomoMeta.mapBackGeometry);
  tomoList = fieldnames(subTomoMeta.mapBackGeometry.tomoName);
  
else  
  if (recon_for_tomoCPR)
    tiltList{1} = varargin{2};
    nTilts = 1;
    % We just need on
    %e valid subtomot
    iTry = 1;
    tomoList{1} = '';
    while iTry < 25
      if (isfield(subTomoMeta.mapBackGeometry.tomoName,sprintf('%s_%d',tiltList{1},iTry)))
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
    tiltRecGeom = cell(nTilts,1);
    for iStack = 1:nTilts
      [ tiltRecGeom{iStack}, tiltName, nTomosPossible, tilt_geometry] = BH_multi_recGeom( sprintf('recon/%s',getCoords(iStack).name), mapBackIter);
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
end


% Divide the tilt series up over each gpu
iterList = cell(nGPUs,1);
% If there is only one tilt, things break in a weird way
nGPUs = min(nGPUs, nTilts);
for iGPU = 1:nGPUs
  iterList{gpuList(iGPU)} = iGPU+(tiltStart-1):nGPUs:nTilts;
  iterList{gpuList(iGPU)};
end

% try
% EMC_parpool(nGPUs)
% catch
% delete(gcp('nocreate'))
% EMC_parpool(nGPUs)
% end

% parfor iGPU = 1:nGPUs
for iGPU = 1:nGPUs %%revert
 
  for iTilt = iterList{gpuList(iGPU)}
    nTomos = 0;
    alreadyMade = 0;
    iTomoList = {};
    
    % For now, since the tilt geometry is not necessarily updated (it is manual)
    % in the subTomoMeta, check that newer (possible perTilt refined) data is
    % not present.
    TLTNAME = sprintf('fixedStacks/ctf/%s_ali%d_ctf.tlt',tiltList{iTilt},mapBackIter+1);
    TLT = load(TLTNAME);
    % fprintf('iGPU %d and iTilt %d using TLT %s\n', iGPU, iTilt, TLTNAME);


    
    % Get all the tomogram names that belong to a given tilt-series.
    % FIXME: I'm not sure it makes sense to restrict this block to for_subTomo
    if (recon_for_subTomo)
      for iTomo = 1:length(tomoList)
        if strcmp(tiltList{iTilt}, subTomoMeta.mapBackGeometry.tomoName.(tomoList{iTomo}).tiltName)
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
% parfor iGPU = 1:nGPUs 
for iGPU = 1:nGPUs %%revert

  % for iGPU = 1:nGPUs
  gpuDevice(gpuList(iGPU));
  % Loop over each tilt
  for iTilt = iterList{gpuList(iGPU)}
    slab_list = {};
    
    % if (recon_for_subTomo)
    %   nTomos = subTomoMeta.mapBackGeometry.(tiltList{iTilt}).nTomos;
    %   iCoords = subTomoMeta.mapBackGeometry.(tiltList{iTilt}).coords;
    %   % FIXME
    % else
    %   if (recon_for_tomoCPR)
    %     nTomos = 1;
    %   else
    %     % templaterch
    %     nTomos  = nTomosPerTilt{iTilt};
    %     iCoords = tiltRecGeom{iTilt};
    %   end
    % end
    
    
    % if (recon_for_subTomo || recon_for_templateMatching))
    %   targetSizeY = diff(floor(iCoords(:,2:3)),1,2)+1;
    %   iCoords = iCoords ./ samplingRate;
    %   iCoords(:,1:4) = floor(iCoords(:,1:4));
    %   iCoords(:,3) = iCoords(:,3) - (diff(floor(iCoords(:,2:3)),1,2)+1 - floor(targetSizeY./samplingRate));
    % end
    % iTomoList = cell(nTomos,1);
    
    

    TLTNAME = sprintf('fixedStacks/ctf/%s_ali%d_ctf.tlt',tiltList{iTilt}, mapBackIter + 1 );
    TLT = load(TLTNAME);
    fprintf('iGPU %d and iTilt %d using TLT %s\n', iGPU, iTilt, TLTNAME);

  
    iTomoList = {};
    if (recon_for_subTomo)
      % Get all the tomogram names that belong to a given tilt-series.
      nTomos = 0;
      alreadyMade = 0;
      for iTomo = 1:length(tomoList)
        if strcmp(tiltList{iTilt},subTomoMeta.mapBackGeometry.tomoName.(tomoList{iTomo}).tiltName)
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

    if (recon_for_templateMatching)
      % templaterch
      nTomos  = nTomosPerTilt{iTilt};
      % tiltRecGeom is a cell with each value being a cell returned by multi_recGeom 
      iCoords = tiltRecGeom{iTilt};
      % iCoords will be a cell indexed by each tomo wwith a struct .tomoCoords
    else
      nTomos = subTomoMeta.mapBackGeometry.(tiltList{iTilt}).nTomos;
      % subTomoMeta.('mapBackGeometry').('tomoCoords').(tomoName).('dX_specimen_to_tomo') = recGeom{tomoIdx}.tomoCoords.dX_specimen_to_tomo;
      iCoords = cell(nTomos,1);
      for iCoordIdx = 1:nTomos
        iCoords{iCoordIdx} = subTomoMeta.mapBackGeometry.tomoCoords.(iTomoList{iCoordIdx});
      % FIXME
      end
    % else
    %   if (recon_for_tomoCPR)
    %     nTomos = 1;
    %   else

    %   end
    end

    
    if samplingRate > 1
      fullStack = sprintf('%aliStacks/%s_ali%d.fixed', ...
        tiltList{iTilt},mapBackIter+1);
      inputStack = sprintf('cache/%s_ali%d_bin%d.fixed',...
        tiltList{iTilt},mapBackIter+1,samplingRate);
      if ~exist(inputStack, 'file')
        BH_multi_loadOrBin(fullStack, samplingRate, 2, false);
      end
    else
      inputStack = sprintf('aliStacks/%s_ali%d.fixed',...
        tiltList{iTilt},mapBackIter+1);
    end
    
    maskedStack = single(getVolume(MRCImage(inputStack)));
    
    if (recon_for_subTomo)
      [ ~, specimen_NX_nm, tomoIdx, ~ ] = calcAvgZ(subTomoMeta,iCoords,tiltList{iTilt}, ...
                                                  iTomoList,nTomos, emc.pixel_size_angstroms, ...
                                                  samplingRate, cycleNumber,...
                                                  0,1);
    else
      if (recon_for_tomoCPR)
        NX = size(maskedStack,1);
        NY = size(maskedStack,2);
        
        %           NY = size(maskedStack,2)-1;
        NZ = floor(reconstructionParameters(1))
        specimen_NX_nm = NZ;
        % iCoords = [NX,0,NY-1,NZ,0,0];
        tomoIdx = 1;
        % FIXME
        for iCoordIdx = 1:nTomos
          iCoords{iCoordIdx}.tomoCoords.dX_specimen_to_tomo = 0;
          iCoords{iCoordIdx}.tomoCoords.dY_specimen_to_tomo = 0;
          iCoords{iCoordIdx}.tomoCoords.dZ_specimen_to_tomo = 0;
          iCoords{iCoordIdx}.tomoCoords.NX = NX * samplingRate;
          iCoords{iCoordIdx}.tomoCoords.NY = NY * samplingRate;
          iCoords{iCoordIdx}.tomoCoords.NZ = NZ * samplingRate;
        end
      else
        [ ~, specimen_NX_nm, tomoIdx, ~ ] = calcAvgZ('dummy',iCoords,tiltList{iTilt}, ...
          iTomoList,nTomos, emc.pixel_size_angstroms, ...
          samplingRate, cycleNumber,...
          0,1);
      end
    end
    
    if ( flg2dCTF || recon_for_tomoCPR)
      n_slabs_to_reconstruct = 1;
      ctf3dDepth = specimen_NX_nm * 10 ^ -9;
    else
      % TODO: for very thick specimen, this may be preventing the avg from getting to high enough
      % resolution to be useful. So far, this is only optimized on in vitro samples.
      dampeningMax = 0.90;

      [ ctf3dDepth ] = BH_ctfCalcError( samplingRate*mean(TLT(:,16)), ...
                                        TLT(1,17),TLT(1,18),abs(TLT(1,15)), ...
                                        2048, TLT(1,19), ...
                                        resTarget,specimen_NX_nm*10, ...
                                        dampeningMax,CYCLE);
      fprintf('\n\nCalculated a ctfDepth of %2.2f nm for %s\n\n',ctf3dDepth*10^9,tiltList{iTilt});
      if (ctf3dDepth > emc.max_ctf3dDepth)
        ctf3dDepth = emc.max_ctf3dDepth;
        fprintf('Calculated ctfDepth exceeds user specified max, so actually using a max_ctfDepth of %2.2f nm\n',ctf3dDepth*10^9);
      end
      % sections centered at 0, which for now is also supposed to coincide with
      % the mean defocus determination, although this could be corrected using
      % knowledge of particle positions given assurance that particles are the
      % primary source of signal (and not carbon for example).
      n_slabs_to_reconstruct = ceil(specimen_NX_nm/(ctf3dDepth*10^9));
      % max odd number
      n_slabs_to_reconstruct = n_slabs_to_reconstruct + ~mod(n_slabs_to_reconstruct,2);
    end
    fprintf('with %3.3f nm sections, correcting %d tilt-series\n', ctf3dDepth*10^9, n_slabs_to_reconstruct);
    
    % For each tomo create a list of slices that are to be reconstructed
    % for every section section.
    
    [ slab_list ] = calc_slab_boundaries(iCoords, tomoIdx, emc.pixel_size_angstroms, n_slabs_to_reconstruct, tiltList{iTilt}, ctf3dDepth, samplingRate);
    

    if (recon_for_subTomo)
      [ avgZ, specimen_NX_nm, tomoIdx, surfaceFit ] = calcAvgZ(subTomoMeta,iCoords,tiltList{iTilt}, ...
        iTomoList,nTomos, emc.pixel_size_angstroms, ...
        samplingRate, cycleNumber,...
        slab_list, 0);
    else
      avgZ = 0;
      surfaceFit = 0;
    end
    
    if ( shiftDefocusOrigin )
      fprintf('Using avgZ %3.3e nm as the defocus origin\n',avgZ*10^9);
    else
      avgZ = 0;
      fprintf('Using sample origin as the defocus origin\n');
      fprintf('If you want to use the COM of subTomos, set shiftDefocusToSubTomoCOM=1\n');
    end
    
    
    first_slab = true(nTomos,1);
    for iSection = 1:n_slabs_to_reconstruct
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
        
        [ correctedStack ] = ctfMultiply_tilt(n_slabs_to_reconstruct,iSection,ctf3dDepth, ...
                                              avgZ,TLT,emc.pixel_size_angstroms,maskedStack,...
                                              specimen_NX_nm*10/emc.pixel_size_angstroms,flgDampenAliasedFrequencies,...
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
        
        thisTomo = tomoIdx(iT);
        
        if (slab_list{iT}(iSection,1))
          reconName = sprintf('%s/%s_ali%d_%d_%d.rec', ...
            tmpCache,tiltList{iTilt},mapBackIter+1,thisTomo,iSection);
            
          if (recon_for_tomoCPR)
            TA = sortrows(subTomoMeta.tiltGeometry.(tomoList{1}),1);
            TA = TA(:,4);
          end

          if (recon_for_subTomo)
            TA = sortrows(subTomoMeta.tiltGeometry.(sprintf('%s_%d',tiltList{iTilt},thisTomo)),1);
            TA = TA(:,4);
          end
          

          if (recon_for_templateMatching)
            if (mapBackIter)
              % FIXME: I don't think this block should work, it should only be the tilt angles!
              error('THis block should not be reached.')
              TA = load(sprintf('%smapBack%d/%s_ali%d_ctf.tlt',CWD,mapBackIter,tiltList{iTilt}, mapBackIter));
            else
              TA = load(sprintf('%sfixedStacks/%s.tlt',CWD,tiltList{iTilt}));
            end
          end
          
          rawTLT = sprintf('cache/%s_%d.rawtlt',tiltList{iTilt},thisTomo);
          rawTLT_file = fopen(rawTLT, 'w');
          fprintf(rawTLT_file,'%f\n', TA');
          fclose(rawTLT_file);
          
          if (mapBackIter)
            LOCAL = sprintf('%smapBack%d/%s_ali%d_ctf.local',CWD,mapBackIter,tiltList{iTilt}, mapBackIter);
          else
            LOCAL = sprintf('%sfixedStacks/%s.local',CWD,tiltList{iTilt});
          end

          if exist(LOCAL,'file')
            flgLocal = 1;
          else
            fprintf('Did not find local alignment information at %s\n',LOCAL);
            flgLocal = 0;
          end
         
          % round down and then we'll add any extra needed to the final chunk
          tiltChunkSize = floor(iCoords{thisTomo}.NY ./ samplingRate ./ emc.n_tilt_workers);
          % This shoulid never happen, but to be safe
          if (emc.n_tilt_workers > floor(iCoords{thisTomo}.NY ./ samplingRate))
            error('n_tilt_workers is greater than the number of slices in the tilt series');
          end
          y_i = floor(iCoords{thisTomo}.y_i ./ samplingRate);
          y_f = floor(iCoords{thisTomo}.y_f ./ samplingRate);


          tiltChunks = y_i:tiltChunkSize:y_f;
          tiltChunks(end) = y_f;
          % Imod expects zero indexed slices
          tiltChunks = tiltChunks - 1;
          totalSlices = [tiltChunks(1),tiltChunks(end)];

          % Make sure we didn't go OOB on the first chunk
          if (tiltChunks(1) < 0)
            tiltChunks(1) = 0;
            n_slices_in_Y = tiltChunks(end) - tiltChunks(1) + 1;
          end
          
         
          rCMD = sprintf(['tilt %s %s -input %s -output %s.TMPPAD -TILTFILE %s -UseGPU %d ', ...
            '-WIDTH %d -COSINTERP 0 -THICKNESS %d -SHIFT %f,%f '],...
            super_sample, ... 
            expand_lines, ...
            outputStack, ...
            reconName, ...
            rawTLT, ...
            gpuList(iGPU), ...
            floor(iCoords{thisTomo}.NX ./ samplingRate),... % WIDTH = NX
            floor(round(slab_list{iT}(iSection,5))), ... % THICKNESS = NZ
            iCoords{thisTomo}.dX_specimen_to_tomo ./ samplingRate, ... % SHIFT X
            slab_list{iT}(iSection,6));
          
          
          reconScaling = 1;
          % Explicitly set Radial to Nyquist
          if (flgLocal)
            rCMD = [rCMD sprintf('-LOCALFILE %s -RADIAL 0.5,.05 -MODE 2 -SCALE 0,%d',LOCAL,reconScaling)];
          else
            rCMD = [rCMD sprintf('-RADIAL 0.5,.05 -MODE 2 -SCALE 0,%d',reconScaling)];
          end
          
          if isfile(sprintf('%s.sh',reconName))
            system(sprintf('rm %s.sh',reconName)); 
          end
          
          recScript = fopen(sprintf('%s.sh',reconName),'w');
          fprintf(recScript,'#!/bin/bash\n\n');
          fprintf(recScript,'%s -SLICE -1,-1 -TOTALSLICES %d,%d\n',rCMD,totalSlices);

          iShift = 1;
          for iChunk = 1:length(tiltChunks)-1
            if (iChunk == length(tiltChunks)-1)
              iShift = 0;
            end
            fprintf(recScript,'%s -SLICE %d,%d -TOTALSLICES %d,%d > /dev/null  &\n', ...
                              rCMD, ...
                              tiltChunks(iChunk),...
                              tiltChunks(iChunk+1) - iShift,...
                              totalSlices);
          end
          fprintf(recScript,'\n\nwait\n\n');
          fclose(recScript);
          pause(1);
          system(sprintf('chmod a=wrx %s.sh',reconName));
          
          [recError,~] = system(sprintf('%s.sh > /dev/null',reconName)); % /dev/null 
          if (recError)
            system(sprintf('%s.sh',reconName));
            error('\n\nerror during reconstruction %s\n\n',reconName);
          end
          % Z coords (y in this orientation) are decreasing into the
          % monitor. For symmetrical padding this doesn't matter, but keep
          % in mind. /dev/null
          trimCMD = sprintf('trimvol -rx -y %d,%d %s.TMPPAD %s > /dev/null  ' , ...
                            1,floor(round(slab_list{iT}(iSection,5))),reconName,reconName);
          [msg,~]= system(trimCMD);
          if (msg)
            fprintf('%d from trimCMD\n',msg)
            trimCMDPrintError = sprintf('trimvol -rx -y %d,%d %s.TMPPAD %s', ...
                                         1,floor(round(slab_list{iT}(iSection,5))),reconName,reconName);
            system(trimCMDPrintError);
            error('error during trimvol');
          end
          system(sprintf('rm %s.TMPPAD', reconName)); 
        end
     
      end % end loop over tomos for this section
      
      if isfile(outputStack)
        system(sprintf('rm %s',outputStack)); 
      end
    end % end loop over sectionsF()

    deltaZ = [];
    evalMask = [];
    maskedStack = [];
    
    for iT = 1:nTomos
      thisTomo = tomoIdx(iT);
      % Note that bh_global_turn_on_phase_plate could be true for any of the recon_for_stage bools, so it must
      % be checked first.
      if (bh_global_turn_on_phase_plate(1))
        reconNameFull = sprintf('cache/%s_%d_bin%d_filtered.rec', ...
          tiltList{iTilt},thisTomo,samplingRate);
      elseif recon_for_tomoCPR
        reconNameFull = sprintf('%scache/%s_%d_bin%d_backgroundEst.rec', ...
          CWD,tiltList{iTilt},thisTomo,samplingRate);
      else
        reconNameFull = sprintf('cache/%s_%d_bin%d.rec', ...
          tiltList{iTilt},thisTomo,samplingRate);
      end
      fprintf('in ctf3d reconNameFull is %s\n\n',reconNameFull);
          
          
      % Get the total number of sections for this tomo
      n_total_sections = 0;
      for iSection = 1:n_slabs_to_reconstruct
        if(slab_list{iT}(iSection,1))
          n_total_sections = n_total_sections + 1;
        end
      end

      if (n_total_sections == 0)
        fprintf('no sections for tomo %d\n',thisTomo);
        continue
      end

      file_of_outputs = sprintf('%s.filelist',reconNameFull);
      recombineCMD = fopen(file_of_outputs,'w');
      fprintf(recombineCMD,'%d\n', n_total_sections);

      cleanup3 = sprintf('rm %s',file_of_outputs);
      % if (use_inverted_newstack)
        % slab_order = n_slabs_to_reconstruct:-1:1;
      % else
        slab_order = 1:n_slabs_to_reconstruct;
      % end
      % for iSection = 1:n_slabs_to_reconstruct
      for iSection = slab_order
        if (slab_list{iT}(iSection,1))
          this_slab = sprintf('%s/%s_ali%d_%d_%d.rec', tmpCache, tiltList{iTilt}, mapBackIter+1, thisTomo, iSection);
          cleanup3 = sprintf('%s %s',cleanup3,this_slab);
          fprintf(recombineCMD, '%s\n', this_slab);
          fprintf(recombineCMD, '1-%d\n',floor(round(slab_list{iT}(iSection,5))));
        end
      end
      fclose(recombineCMD);
      pause(1);
      recCMD = sprintf('newstack -fromone -FileOfInputs %s -output %s\n', file_of_outputs, reconNameFull);
    
      [err_msg, ~] = system(sprintf('%s > /dev/null ',recCMD)); %/dev/null
      if (err_msg)
        fprintf('error during recombination %s\n',reconNameFull);
        system(recCMD);
        error('error during recombination');
      end
      
      system(cleanup3);  
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


function  [ slab_list ] = calc_slab_boundaries(iCoords, tomoIdx, pixel_size_angstroms, n_slabs_to_reconstruct, tiltName, ctf_3d_depth_si, samplingRate, use_inverted_newstack)

  %%% This function is to produce a list of z indices, starting from 1, to pass to imod for tilt based reconstruction
nTomos = length(tomoIdx);
slab_list = cell(nTomos,1);
for iTomo = 1:nTomos
  % min and max in absolute pixels min and max from 1:nZrecon
  slab_list{iTomo} = zeros(n_slabs_to_reconstruct,6);
end

% With rounding this could end up a bit short except the top and bottom are both
% half a section larger than minimally needed.
slab_size_pixels = floor(ctf_3d_depth_si * 10^10 / pixel_size_angstroms);
slab_size_pixels = slab_size_pixels + ~mod(slab_size_pixels,2);
oS = emc_get_origin_index(slab_size_pixels);

for iT = 1:length(tomoIdx)
  iTomo = tomoIdx(iT);
  % Origin + originshift
  tomo_origin_wrt_tilt_origin = iCoords{iT}.dZ_specimen_to_tomo ./ samplingRate;              
  tomo_origin_in_tomo_frame = emc_get_origin_index(iCoords{iT}.NZ ./ samplingRate); 
                                                    
  fraction_origin_shift = tomo_origin_wrt_tilt_origin - round(tomo_origin_wrt_tilt_origin);
  
  tomogram_lower_bound = floor((tomo_origin_wrt_tilt_origin - tomo_origin_in_tomo_frame));
  recon_range_z_in_specimen_frame = tomogram_lower_bound : tomogram_lower_bound + ceil(iCoords{iT}.NZ./samplingRate) - 1;
  % For each slab see if this tomogram has any sections in it
  for iSlab = 1:n_slabs_to_reconstruct
  
    slab_origin_in_specimen_frame = ((n_slabs_to_reconstruct-1)/-2+(iSlab-1)) * slab_size_pixels;
    
    slab_lower_bound = slab_origin_in_specimen_frame - oS + 1;
    slab_upper_bound = slab_origin_in_specimen_frame + oS - 1;
    slab_range = slab_lower_bound:slab_upper_bound;

    is_in_range = ismember(recon_range_z_in_specimen_frame, slab_range);
    valid_indices = recon_range_z_in_specimen_frame(is_in_range);

    slab_list{iT}(iSlab,5) = length(valid_indices);
    if (slab_list{iT}(iSlab,5) > 0)
      slab_list{iT}(iSlab,1) = 1;
    else
      continue;
    end
    valid_region_origin = emc_get_origin_index(slab_list{iT}(iSlab,5));
    % This is a vector from the origin of the sample to the origin of the slab.
    % The shift passed to imod-tilt moves the reconstructed area in the opposite sense.
    % All slabs need to be shifted to the specimen origin (0) from tilts perspective, and then the are assembled into the final volume.
    % This means a slab at Z > 0 needs to be shifted in the negative direction, which means supplying 
    % a shift that is also > 0, moving the volume "up" in the rotated coordinate system (imod -Z)
    % I know ... this is a shit show.
    dZ_for_reconstructed_slab = (valid_indices(valid_region_origin) + fraction_origin_shift);

    slab_list{iT}(iSlab,6) = dZ_for_reconstructed_slab; %dZ
  end

  % Check to ensure we don't have any tiny slabs leftover, if so, merge them into a neighboring slab
  biggest_slab = max(slab_list{iT}(:,5));
  for iSlab = 1:n_slabs_to_reconstruct
    if (slab_list{iT}(iSlab,1) && slab_list{iT}(iSlab, 5) / biggest_slab < 0.1)
      if (iSlab > 1 && slab_list{iT}(iSlab-1,1))
        delta = slab_list{iT}(iSlab,5);
        slab_list{iT}(iSlab-1,5) = slab_list{iT}(iSlab-1,5) + delta;
        slab_list{iT}(iSlab,1) = 0;
        % we are adding slices from above the specimen in Z so the z shift is positive
        slab_list{iT}(iSlab-1,6) = (slab_list{iT}(iSlab-1,6) + delta);
      elseif (iSlab < n_slabs_to_reconstruct && slab_list{iT}(iSlab+1,1))
        delta = slab_list{iT}(iSlab,5);
        slab_list{iT}(iSlab+1,5) = slab_list{iT}(iSlab+1,5) + slab_list{iT}(iSlab,5);
        slab_list{iT}(iSlab,1) = 0;
        % we are adding slices from below the specimen in Z so the z shift is negative
        slab_list{iT}(iSlab+1,6) = (slab_list{iT}(iSlab+1,6) - delta);
      end
    end
  end

  % TroubleShoot
  tSHT = fopen(sprintf('.tblSht_%s_i%d.txt',tiltName,iTomo),'w');
  fprintf(tSHT,'%2.2f %2.2f %2.2f %2.2f %2.2f %2.2f\n',  iCoords{iT}.NX, iCoords{iT}.NY, iCoords{iT}.NZ, iCoords{iT}.dX_specimen_to_tomo, iCoords{iT}.dY_specimen_to_tomo, iCoords{iT}.dZ_specimen_to_tomo);
  fprintf(tSHT,'%2.2f %2.2f %2.2f %2.2f %2.2f %2.2f\n', slab_list{iT}');
  fclose(tSHT);
end % end loop over tomos



end




function [correctedStack] = ctfMultiply_tilt(n_slabs_to_reconstruct,iSection,ctf3dDepth, ...
                                            avgZ,TLT,pixel_size_angstroms,maskedStack,...
                                            specimen_NX_nm,flgDampenAliasedFrequencies,...
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
    % fprintf('Warning, surfaceFit is not an sfit object\n');
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
  defocusOffset = (((n_slabs_to_reconstruct-1)/-2+(iSection-1))*ctf3dDepth);
  fprintf('Not using surface fit, so using offset %3.3e nm for section %d with COM offset %3.3e nm with ctf3dDepth %3.3e\n', defocusOffset*10^9, iSection, avgZ*10^9, ctf3dDepth*10^9);
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

if (filterProjectionsForTomoCPRBackground ~= 0)
  bpFilter = BH_bandpass3d(fastFTSize,0, 0, filterProjectionsForTomoCPRBackground, 'GPU',pixel_size_angstroms);
  fprintf('Filtering input projections to %f angstroms with %f pixel size\n',bpFilter,pixel_size_angstroms);

else
  bpFilter = 1;
end


for iPrj = 1:nPrjs
  
  maxEval = cosd(TLT(iPrj,4)).*(d1/2) + specimen_NX_nm./2*abs(sind(TLT(iPrj,4)));
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
  rA =  BH_defineMatrix(TLT(iPrj,4),'TILT','fwdVector') ;

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

function [avgZ, specimen_NX_nm, tomoIdx,surfaceFit] = calcAvgZ(subTomoMeta, ...
                                                               iCoords, ...
                                                               tiltName, ...
                                                               tomoList,...
                                                               nTomos,  ...
                                                               pixel_size_angstroms,...
                                                               samplingRate, ...
                                                               cycleNumber,...
                                                               slab_list, ...
                                                               calcMaxZ)

% Calculate the maximum extensions in Z and then how many separate sections
% need to be corrected.

surfaceFit = '';
avgZ = 0;


if isa(subTomoMeta,'struct')
  val_to_pass =  subTomoMeta.mapBackGeometry.tomoName;
else
  val_to_pass = 'dummy';
end

[ specimen_NX_nm, tomoIdx ] = emc_get_max_specimen_NZ(val_to_pass, ...
                                                      iCoords, ...
                                                      tomoList, ...
                                                      nTomos, ...
                                                      samplingRate);

specimen_NX_nm = specimen_NX_nm .* pixel_size_angstroms ./ 10;
fprintf('combining the thickness and shift on tilt %s, found a specimen_NX_nm  %3.3f nm\n',tiltName,specimen_NX_nm);

if (calcMaxZ)
  return;
end

% For now use cycle000, if adding a refinment focused on a specific set of
% particles, then consider that later.
try
  initGeom = subTomoMeta.(cycleNumber).RawAlign;
  fprintf('Loaded the geometry for RawAlign %s\n',cycleNumber);
catch
  fprintf('Failed to load the geometry for RawAlign %s\nTrying cycle000\n',cycleNumber);
  try
    initGeom = subTomoMeta.cycle000.geometry;
  catch
    error(['Could not load the initial geometry subTomoMeta.%s.geometry\n or--',...
      'subTomoMeta.cycle000.geometry\n'],cycleNumber);
  end
end

nSubTomos = 0;
totalZ = 0;
% for each tomogram get the size and origin in Z then find mean subTomo
% position.

n_slabs_to_reconstruct = size(slab_list{1},1);
xFull = cell(n_slabs_to_reconstruct,1);
yFull = cell(n_slabs_to_reconstruct,1);
zFull = cell(n_slabs_to_reconstruct,1);
surfaceFit = cell(n_slabs_to_reconstruct,1);

% Initialize with empty arrays
for iSection = 1:n_slabs_to_reconstruct
  xFull{iSection} = [];
  yFull{iSection} = [];
  zFull{iSection} = [];
end


for iT = 1:nTomos

  iTomo = tomoIdx(iT);

  % X in the Y frame means a vector from the Y lower left to the X origin
  % X origin wrt Y origin is a vector from the origin of Y to the X origin
  reconGeometry = subTomoMeta.mapBackGeometry.tomoCoords.(tomoList{iTomo});
  tomo_origin_wrt_tilt_origin = [reconGeometry.dX_specimen_to_tomo, ...
                                reconGeometry.dY_specimen_to_tomo, ...
                                reconGeometry.dZ_specimen_to_tomo];              
  tomo_origin_in_tomo_frame = emc_get_origin_index([reconGeometry.NX, ...
                                                    reconGeometry.NY, ...
                                                    reconGeometry.NZ]); 
 
  iTomoName = sprintf('%s_%d',tiltName,iTomo);
  
  % Get the z-coordinates of the origin for all included subtomograms relative to the lower left of the tomogram
  % shouldn't be any removed particles at this stage but later there would be.
  subtomo_origin_z_in_tomo_frame = initGeom.(iTomoName)(initGeom.(iTomoName)(:,26)~=-9999,13);

  % We should not get to this point if all subtomograms have been removed.
  if isempty(subtomo_origin_z_in_tomo_frame)
    error('No subtomograms found for %s',iTomoName);
  end

  % shift from lower left to centered and include the tomos offset from the
  subtomo_origin_wrt_specimen_origin = subtomo_origin_z_in_tomo_frame - tomo_origin_in_tomo_frame(3) + tomo_origin_wrt_tilt_origin(3);
  subtomo_origin_wrt_specimen_origin = subtomo_origin_wrt_specimen_origin ./ samplingRate;
  totalZ = totalZ + sum(subtomo_origin_wrt_specimen_origin);
  fprintf('%s tomo has %d subTomos with mean Z %3.3f nm\n', ...
    iTomoName, length(subtomo_origin_wrt_specimen_origin), ...
    mean(subtomo_origin_wrt_specimen_origin) * pixel_size_angstroms ./ 10);

  nSubTomos = nSubTomos + length(subtomo_origin_wrt_specimen_origin);

  for iSection = 1:n_slabs_to_reconstruct
    
    iSecOrigin = slab_list{iT}(iSection,6);
    iSecRadius = slab_list{iT}(iSection,5)/2;
    inSectionIDX = subtomo_origin_wrt_specimen_origin >  iSecOrigin - iSecRadius & subtomo_origin_wrt_specimen_origin <= iSecOrigin + iSecRadius;
    
    
    x = initGeom.(iTomoName)(initGeom.(iTomoName)(:,26)~=-9999,11);
    x = (x - tomo_origin_in_tomo_frame(1) + tomo_origin_wrt_tilt_origin(1))./samplingRate;
    y = initGeom.(iTomoName)(initGeom.(iTomoName)(:,26)~=-9999,12);
    y = (y - tomo_origin_in_tomo_frame(2) + tomo_origin_wrt_tilt_origin(2))./samplingRate;
    
    
    xFull{iSection} = [xFull{iSection} ; x(inSectionIDX)];
    yFull{iSection} = [yFull{iSection} ; y(inSectionIDX)];
    zFull{iSection} = [zFull{iSection} ; subtomo_origin_wrt_specimen_origin(inSectionIDX)];
    
  end % loop over sections
  clear subtomo_origin_wrt_specimen_origin
end % loop over tomos


avgZ = totalZ / nSubTomos*pixel_size_angstroms / 10*10^-9;


for iSection = 1:n_slabs_to_reconstruct
  
  if length(xFull{iSection}) >= 6
    surfaceFit{iSection} = fit([xFull{iSection}, yFull{iSection}],zFull{iSection},'poly22','Robust','on');
  else
    surfaceFit{iSection} = 0;
  end
end

fprintf('%s tilt-series has %d subTomos with mean Z %3.3f nm\n', tiltName, nSubTomos,avgZ*10^9);

end
