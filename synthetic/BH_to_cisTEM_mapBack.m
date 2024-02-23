function [ ] = BH_to_cisTEM_mapBack(PARAMETER_FILE, CYCLE, output_prefix, symmetry, MAX_EXPOSURE, varargin)

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


emc = BH_parseParameterFile(PARAMETER_FILE);
MAX_EXPOSURE = EMC_str2double(MAX_EXPOSURE)
if isnan(MAX_EXPOSURE)
  error('MAX_EXPOSURE is nan - if running from an interactive matlab session, did you enter as a string?');
end

% Ideally, we would transform fully and go back to the non-rotated stack. I think with the apoferritin test set,
% The resolution will be high-enough to sort this out.
useFixedNotAliStack = false;

CYCLE = EMC_str2double(CYCLE);
cycle_numerator = '';
cycle_denominator ='';

% Copied over from tomoCPR. I think I can use this to write out all the partial stacks,
% Then rather than "runAlignments" have a command that combines all the stacks.
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

% Always working at full binning, not emc.('Ali_samplingRate');
reconScaling = 1;

% TODO: use this and add a block to calculate the FSC of the output reconstruction prior to refinement
MOL_MASS = emc.('particleMass');



% Used to calc defocus values using tilt instead of manually. Convention
% diff.
flgInvertTiltAngles = 0;


if (skip_to_the_end_and_run)
  % The fractional runs have already copied everything to cache/mapback%d,
  % so override the tmpCache.
  tmpCache = '';
else
  tmpCache= emc.('fastScratchDisk');
end

[tmpCache, flgCleanCache, CWD] = EMC_setup_tmp_cache(tmpCache, '', 'cisTEM', true);

nGPUs = emc.('nGPUs');
pInfo = parcluster();
gpuScale=3;
nWorkers = min(nGPUs*gpuScale,emc.('nCpuCores')); % 18
fprintf('Using %d workers as max of %d %d*nGPUs and %d nWorkers visible\n', ...
  nWorkers,gpuScale,nGPUs*gpuScale,pInfo.NumWorkers);


load(sprintf('%s.mat', emc.('subTomoMeta')), 'subTomoMeta');
mapBackIter = subTomoMeta.currentTomoCPR;

% TODO: use these to add an optional defocus fitting step
% So translational, optional defocus, angles
ctfRange = emc.('tomo_cpr_defocus_range')*10^10;
ctfInc = emc.('tomo_cpr_defocus_step')*10^10;
calcCTF = emc.('tomo_cpr_defocus_refine');


[tiltNameList, nTiltSeries] = BH_returnIncludedTilts( subTomoMeta.mapBackGeometry );


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


load(sprintf('%s.mat', emc.('subTomoMeta')), 'subTomoMeta');
resForFitting = 1.3*mean(subTomoMeta.currentResForDefocusError);
tiltGeometry = subTomoMeta.tiltGeometry;

% TODO: this is a bit of an aritfact, can be removed.
outCTF = '_ctf';

is_first_run = true;

mbOUT = {[tmpCache],'dummy'};
tiltStart=1; 
firstTilt = true;



iCell = 0;
output_cell = {};
newstack_file = sprintf('%s/temp_particle_stack.newstack',mbOUT{1});
newstack_file_handle = fopen(newstack_file,'w');


for iTiltSeries = tiltStart:nTiltSeries
  n_particles_added_to_stack = 0;
  if (skip_to_the_end_and_run)
    continue;
  end
  
  if (useFixedNotAliStack)
    tilt_filename = sprintf('%sfixedStacks/%s.fixed',CWD,tiltNameList{iTiltSeries});
  else
    tilt_filename = sprintf('%saliStacks/%s_ali%d.fixed',CWD,tiltNameList{iTiltSeries},mapBackIter+1);
  end
  tilt_filename = sprintf('%saliStacks/%s_ali%d.fixed', CWD, tiltNameList{iTiltSeries}, mapBackIter + 1);

  mapBackRePrjSize = subTomoMeta.mapBackGeometry.(tiltNameList{iTiltSeries}).('tomoCprRePrjSize');
  % % %   iViewGroup = subTomoMeta.mapBackGeometry.viewGroups.(tiltNameList{iTiltSeries});
  nTomograms = subTomoMeta.mapBackGeometry.(tiltNameList{iTiltSeries}).nTomos
  if nTomograms == 0
    % No points were saved after template matching so skip this tilt series
    % altogether.
    continue
  end
  
  skip_this_tilt_series_because_it_is_empty = false(nTomograms,1);
  
  % tomoList = fieldnames(subTomoMeta.mapBackGeometry.tomoName);
  tomoList = {};
  tomoIDX = 1;
  fn = fieldnames(subTomoMeta.mapBackGeometry.tomoName);
  for iTomo = 1:numel(fn)
    if (subTomoMeta.mapBackGeometry.tomoName.(fn{iTomo}).tiltName == tiltNameList{iTiltSeries})
      % This is dumb, fix it to be explicit.
      if (subTomoMeta.mapBackGeometry.tomoCoords.(fn{iTomo}).is_active)
        tomoList{tomoIDX} = fn{iTomo};
        % Only increment if values found.
        tomoIDX = tomoIDX + 1;
      end
    end
  end
    
  [~,tiltBaseName,~] = fileparts(tilt_filename);
  mbOUT{2} = tiltBaseName;
  
  if (mapBackIter)
    localFile = sprintf('%s/%s_ali%d_ctf.local', CWD,mapBackIter,tiltNameList{iTiltSeries},mapBackIter);
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
  pixel_size = emc.pixel_size_angstroms;
  
  try
    eraseMaskType = emc.('Peak_mType');
    eraseMaskRadius = emc.('Peak_mRadius') ./ pixel_size;
    fprintf('Further restricting peak search to radius of [%f %f %f] pixels\n', eraseMaskRadius);
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
  try
    lowPassCutoff = emc.('tomoCprLowPass');
    fprintf('Using a user supplied lowpass cutoff of %3.3f Ang\n', lowPassCutoff);
  catch
    % TODO are these range limits okay?
    lowPassCutoff = 1.5.*mean(subTomoMeta.currentResForDefocusError);
    if (lowPassCutoff < 10)
      lowPassCutoff = 10;
    elseif (lowPassCutoff > 24)
      lowPassCutoff = 24;
    end
    fprintf('Using an internatlly determined lowpass cutoff of %3.3f Ang\n',...
      lowPassCutoff);
  end

  % FIXME: this can also be in parseParameterFile 
  if lowPassCutoff < 2* pixel_size
    fprintf('Psych, the cutoff is being set to Nyquist');
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
      fprintf('Warning the current resolution is too low to refine the defocus. Turning off this feature');
      calcCTF = false;
    end
  end
  

   
  % Get the thickest for recon
  maxZ = 0;

  % The 
  tiltHeader = getHeader(MRCImage(tilt_filename, 0));
  tiltName = subTomoMeta.mapBackGeometry.tomoName.(tomoList{1}).tiltName;
  [ maxZ ] = emc_get_max_specimen_NZ(subTomoMeta.mapBackGeometry.tomoName, ...
                                     subTomoMeta.mapBackGeometry.tomoCoords,  ...
                                     tomoList, ...
                                     nTomograms, ...
                                     1);

  fprintf('combining thickness and shift, found a maxZ of %d\n',maxZ);
  
  % xyzproj assumes centered in Z, so add extra height for z offsets to create
  % the true "in microsope" dimension
  
  reconstruction_size = [tiltHeader.nX, tiltHeader.nY, maxZ];
  originRec = emc_get_origin_index(reconstruction_size);

  TLT = tiltGeometry.(tomoList{1});
  
  
  iRawTltName = sprintf('%s/%s_align.rawtlt',mbOUT{1:2})
  iTiltFile = fopen(iRawTltName, 'w');
  rawTLT = sortrows(TLT(:,[1,4]),1);
  fprintf(iTiltFile,'%f\n',rawTLT(:,2)');
  fclose(iTiltFile);
  
  coordOUT = fopen(sprintf('%s/%s.coord',mbOUT{1:2}),'w');
  coordSTART = fopen(sprintf('%s/%s.coord_start',mbOUT{1:2}),'w');
  
  defOUT   = fopen(sprintf('%s/%s.defAng',mbOUT{1:2}),'w');
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

    % Extract a "defocus file" for tilt to calculate the defocus for each
    % fiducial also considering the local alignment. If this works, I can
    % get rid of defAng
    iDefocusFileName = sprintf('%s/%s_align.defocus',mbOUT{1:2});
    iDefocusFile = fopen(iDefocusFileName,'w');
    defTLT = sortrows(TLT(:,[1,15]),1);
    fprintf(iDefocusFile,'%f\n',abs(defTLT(:,2)').*10^9);
    fclose(iDefocusFile);
    
    % We also need the transform from the microscope frame in order to
    % get an accurate defocus value. Not sure if I should be binning?
    % Additionally, we do NOT want the model for alignment in the
    % microscope frame,
    iXFName = sprintf('%s/%s_align.XF',mbOUT{1:2});
    iXF = fopen(iXFName,'w');
    
    if (useFixedNotAliStack)
      
      % 20190509 - I think this is ry_startally screwing things up FIXME
      % Commenting this out invalidates the defocus vals
      % positionn in stack, imod rotation matrix (2x2), x,y shift (unbinned)
      xfTLT = sortrows(TLT(:,[1,7:10,2,3],1));
      fprintf(iXF,'%f %f %f %f %f %f\n',xfTLT(:,2:7)');
      fclose(iXF);
      % Odd size stacks are enforced which creates a shift prior to the
      % xform.
      if (useFixedNotAliStack)
        isEven = 1;
        
        iXFBase = sprintf('%s/%s_align_base.XF',mbOUT{1:2});
        iXFB = fopen(iXFBase,'w');
        for ix = 1:size(xfTLT,1)
          fprintf(iXFB,'%f %f %f %f %f %f\n',[1,0,0,1,-isEven,-isEven]);
        end
        
        fclose(iXFB);
        system(sprintf('xfproduct %s %s %s',iXFBase, iXFName,iXFName));
        
        % We need to invert this transform to map from the aligned stack to the
        % fixed stack
        iXFName_inv = sprintf('%s/%s_align_inv.XF',mbOUT{1:2});
        system(sprintf('xfinverse %s %s', iXFName, iXFName_inv));
      end
    else
      % Create an identity transform for the model
      % 20190509 - I think this is ry_startally screwing things up FIXME
      % Commenting this out invalidates the defocus vals
      xfTLT = zeros(size(TLT,1),6);
      xfTLT(:,[1,4]) = 1.0;
      fprintf(iXF,'%f %f %f %f %f %f\n',xfTLT');
      fclose(iXF);
      
      %       xfTLT = sortrows(TLT(:,[1,7:10,2,3],1));
      %       fprintf(iXF,'%f %f %f %f %f %f\n',xfTLT(:,2:7)');
      %       fclose(iXF);
    end
     
    positionList = geometry.(tomoList{iTomo});
    tomoIdx = subTomoMeta.mapBackGeometry.tomoName.(tomoList{iTomo}).tomoIdx;
    tiltName   = subTomoMeta.mapBackGeometry.tomoName.(tomoList{iTomo}).tiltName;
    
    positionList = positionList(positionList(:,26) ~= -9999,:);
    nFidsTotal = nFidsTotal + size(positionList,1);

    tiltHeader = getHeader(MRCImage(tilt_filename,0));
    
    fullTiltSizeXandY = [tiltHeader.nX,tiltHeader.nY];
    
    sTX = floor(tiltHeader.nX);
    sTY = floor(tiltHeader.nY);    
        
    tomoIdx = subTomoMeta.mapBackGeometry.tomoName.(tomoList{iTomo}).tomoIdx;
    tiltName = subTomoMeta.mapBackGeometry.tomoName.(tomoList{iTomo}).tiltName;
    reconGeometry = masterTM.mapBackGeometry.tomoCoords.(tomoList{iTomo});    
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

    modelRot = BH_defineMatrix([0,90,0],'Bah','fwdVector');
    
    for iSubTomo = 1:nSubTomos
      
      subtomo_rot_matrix = reshape(positionList(iSubTomo,17:25),3,3);
      subtomo_origin_in_tomo_frame = (positionList(iSubTomo,11:13));
      subtomo_origin_wrt_tilt_origin = subtomo_origin_in_tomo_frame - tomo_origin_in_tomo_frame + tomo_origin_wrt_tilt_origin;
        
      % This extra shift came from experiments with real data but is both anny_starting and not understood.
      subtomo_origin_wrt_tilt_origin = subtomo_origin_wrt_tilt_origin - emc.flgPreShift;

      % subTomo origin relative to reconLowerLeft
      subtomo_origin_in_sample = originRec + subtomo_origin_wrt_tilt_origin; 
      % Reproject using tilt, so just save the 3d coords.
      fprintf(coordOUT,'%0.4f %0.4f %0.4f %d\n', modelRot * subtomo_origin_wrt_tilt_origin' + [originRec(1),originRec(3),originRec(2)]'- emc.prjVectorShift([1,3,2])', fidIDX);
      
      nPrjsIncluded = 0;
      for iPrj = 1:nPrjs
        
        iPrj_nat = find(TLT(:,1) == iPrj);
        if (abs(TLT(iPrj_nat,11)) <= MAX_EXPOSURE)
          nPrjsIncluded = nPrjsIncluded + 1;
          
          % imod is indexing from zero
          zCoord = iPrj_nat;
          % For a positive angle, this will rotate the positive X axis farther from the focal plane (more underfocus)
          rTilt = BH_defineMatrix(TLT(iPrj_nat,4),'TILT','fwdVector') ;
          
          prjCoords = rTilt*subtomo_origin_wrt_tilt_origin';
          
          % I think this is for comparison with the values obtained from projecting using IMOD: FIXME
          fprintf(defOUT,'%d %d %6.6e\n', fidIDX, zCoord, abs(TLT(iPrj_nat,15)) - prjCoords(3).*pixel_size.*10^-10);
          
          % Defocus value adjusted for Z coordinate in the tomogram. nm
          d1 = (abs(TLT(iPrj_nat,15)) - subtomo_origin_wrt_tilt_origin(3).*pixel_size.*10^-10) * 10^9;
          d2 = TLT(iPrj_nat,12)*10^9; % half astigmatism value
          
          fprintf(coordSTART,'%d %d %d %3.3f %3.3f %3.3f %3.3f %3.3f %3.3f %3.3f %3.3f %3.3f %3.3f %3.3f %3.3f %3.3f %3.3f %d\n', ...
                                fidIDX, tomoIdx, positionList(iSubTomo,4), d1, d2, 180./pi.*TLT(iPrj_nat,13), reshape(subtomo_rot_matrix,1,9), preExposure(iPrj_nat), postExposure(iPrj_nat), positionList(iSubTomo,7));
        else
          fprintf(coordSTART,'%d %d %d %3.3f %3.3f %3.3f %3.3f %3.3f %3.3f %3.3f %3.3f %3.3f %3.3f %3.3f %3.3f %3.3f %3.3f %d\n',-9999, -9999,-9999,1.0,1.0,1.0,1,1,1,1,1,1,1,1,1,0,0,1);
        end
        % nFidsTotalDataSet = nFidsTotalDataSet + 1;
      end % loop over tilt projections
      
      fidIDX = fidIDX + 1;
      
    end % loop over subtomos  
  end % end of loop over tomograms on this tilt-series
  
  % No subtomos remain
  if all( skip_this_tilt_series_because_it_is_empty )
    continue;
  end
  
  fclose(coordOUT);
  fclose(coordSTART);
  
  p2m = sprintf(['point2model -zero -circle 3 -color 0,0,255 -values -1 ',...
                '%s/%s.coord %s/%s.3dfid'], ...
                mbOUT{1:2},mbOUT{1:2});
  system(p2m);
  
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
  
  
  if (lastLine1)
    lastLine2 = 'UseGPU 0';
    lastLine3 = 'ActionIfGPUFails 2,2';
  else
    lastLine1 = 'UseGPU 0';
    lastLine2 = 'ActionIfGPUFails 2,2';
    lastLine3 = '';
  end
  
  
  
  % Break this up into chunks since things hang even with the
  % ActionIfGPUFails option. Try 3 times 512,256,128
  %         refPrj = zeros(sTX,sTY,iTLT, 'single');
  
  iSave = 1;
  reModFileName = sprintf('%s/%s_%d_reMod.sh',mbOUT{1:2},iSave);
  reModFile = fopen(reModFileName,'w');
  invertTiltAngles = 0;
  fprintf(reModFile,['#!/bin/bash\n\n',...
    'tilt -StandardInput << EOF\n',...
    'input %s\n', ...
    'output %s/%s.fid\n', ...
    'COSINTERP 0\n', ...
    'THICKNESS %d\n', ...
    'TILTFILE %s/%s_align.rawtlt\n', ...
    'DefocusFile %s/%s_align.defocus\n', ...
    'PixelForDefocus %f,%f\n', ...
    'AngleOutputFile %s/%s.defAngTilt\n', ...
    'AlignTransformFile %s/%s_align.XF\n', ...
    'ProjectModel %s/%s.3dfid\n', ...
    '%s\n',...
    '%s\n',...
    '%s\n',...
    'EOF'],tilt_filename, mbOUT{1:2}, maxZ, ...
    mbOUT{1:2},...
    mbOUT{1:2},...
    pixel_size./10, ... % Ang --> nm
    0, ... % do not invert the tilt angles
    mbOUT{1:2},...
    mbOUT{1:2},...
    mbOUT{1:2},...
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
  
  if (useFixedNotAliStack)
    % transform the projected model back to the fixed stack frame, and
    % then convert to text.
    
    system(sprintf('imodtrans -2 %s %s/%s.fid %s/%s.invfid  > /dev/null', iXFName_inv, mbOUT{1:2},mbOUT{1:2}));
    
    system(sprintf(['model2point -contour -zero ',...
          '%s/%s.invfid %s/%s.coordPrj > /dev/null'],...
          mbOUT{1:2}, mbOUT{1:2}))
  else
    system(sprintf(['model2point  -contour -zero ',...
                '%s/%s.fid %s/%s.coordPrj > /dev/null'],...
                mbOUT{1:2}, mbOUT{1:2}))
  end
  
  
  
  try
    fidList = load(sprintf('%s/%s.coordPrj',mbOUT{1:2}));
  catch
    error('\nWarning, did not load the projected coords\nSkipping along');
  end
  % unused
  parList = load(sprintf('%s/%s.coord_start',mbOUT{1:2}));
  defList = load(sprintf('%s/%s.defAngTilt',mbOUT{1:2}));

  %   Need to shift again from the model coordinate system
  % Columns are 
  % particle idx, from 0
  % x
  % y
  % projection idx, from 0

  fidList(:,[2,3]) = fidList(:,[2,3]) + repmat(emc.prjVectorShift(1:2), size(fidList,1),1);
  foundNans = sum(isnan(fidList(:,3)));
  if (foundNans)
    fprintf('\n\t\tThere are %d NaNs in the projected fiducial list %3.3f\n\n',foundNans, foundNans/size(fidList,1)*100);
    fprintf('The only confirmed case that produced this were NaNs in the fixedStacks/tiltN.local file.\n');
    error("Exiting");
  end
  
  % Give every instance of each fiducial a unique identifier.
  fidList = [1:size(fidList,1);fidList']';
  
  particlePad = 2.0;
  tileRadius = floor(particlePad.*particle_radius);
  tileSize = BH_multi_iterator((2.*tileRadius).*[1,1],'fourier2d');
  tileOrigin = emc_get_origin_index(tileSize);
  
  nFidsTotal = numel(unique(fidList(parList(:,1)~=-9999,2)));

  
  % This will need to be changed to aggregate
  output_particle_stack = zeros([tileSize,nFidsTotal*nPrjsIncluded],'single');
  
  iGpuDataCounter = 1;
  

  if (firstTilt)
    iDataCounter = 1;
    starFile = fopen(sprintf('%s.star',output_prefix),'w');
    fprintf(starFile, [ ...
      '# Written by emClarity Version 2.0.0-alpha on %s\n\n' ...
      'data_\n\n' ...
      'loop_\n\n' ...
      '_cisTEMPositionInStack #1\n' ...
      '_cisTEMAnglePsi #2\n' ...
      '_cisTEMAngleTheta #3\n' ...
      '_cisTEMAnglePhi #4\n' ...
      '_cisTEMXShift #5\n' ...
      '_cisTEMYShift #6\n' ...
      '_cisTEMDefocus1 #7\n' ...
      '_cisTEMDefocus2 #8\n' ...
      '_cisTEMDefocusAngle #9\n' ...
      '_cisTEMPhaseShift #10\n' ...
      '_cisTEMOccupancy #11\n' ...
      '_cisTEMLogP #12\n' ...
      '_cisTEMSigma #13\n' ...
      '_cisTEMScore #14\n' ...
      '_cisTEMScoreChange #15\n' ...
      '_cisTEMPixelSize #16\n' ...
      '_cisTEMMicroscopeVoltagekV #17\n' ...
      '_cisTEMMicroscopeCsMM #18\n' ...
      '_cisTEMAmplitudeContrast #19\n' ...
      '_cisTEMBeamTiltX #20\n' ...
      '_cisTEMBeamTiltY #21\n' ...
      '_cisTEMImageShiftX #22\n' ...
      '_cisTEMImageShiftY #23\n' ...
      '_cisTEMBest2DClass #24\n' ...
      '_cisTEMBeamTiltGroup #25\n' ...
      '_cisTEMParticleGroup #26\n' ...
      '_cisTEMPreExposure #27\n' ...
      '_cisTEMTotalExposure #28\n' ...
      '#    POS     PSI   THETA     PHI       SHX       SHY      DF1      DF2  ANGAST  PSHIFT     OCC      LogP      SIGMA   SCORE  CHANGE    PSIZE    VOLT      Cs    AmpC  BTILTX  BTILTY  ISHFTX  ISHFTY 2DCLS  TGRP    PARGRP  PREEXP  TOTEXP\n' ...
      ], datetime);

    firstTilt = false;
  end
  
  if (useFixedNotAliStack)
    fullXform = load(iXFName_inv);
  end
  
  STACK = single(getVolume(MRCImage(tilt_filename)));
  
  for iPrj = 1:nPrjs
    
    if (abs(TLT(iPrj,11)) > MAX_EXPOSURE)
      continue;
    end
    

    % Both are ordered by fiducial (imod contour number) but are not
    % explicitly checked to correspond. Should this be done?
    
    % fid list is produced by projection of the 3dmodel using tilt with
    wrkPrjIDX = ( fidList(:,5) == TLT(iPrj,1) - 1 );
    wrkFid = fidList(wrkPrjIDX,:);
    wrkPar = parList(wrkPrjIDX,:);
    wrkDefAngTilt = defList(wrkPrjIDX,[7,6,5]); % Confirming with David but this should include the local adjustments to tilt/in-plane angle
    
    
    for iFid = 1:size(wrkFid,1)
      
      if (wrkPar(iFid,1) == -9999)
        continue;
      end

      pixelX = wrkFid(iFid,3) - emc.pixelShift + emc.flgPostShift(1);
      pixelY = wrkFid(iFid,4) - emc.pixelShift + emc.flgPostShift(2);
      
      x_start = floor(pixelX) - tileOrigin(1) + 1;
      y_start = floor(pixelY) - tileOrigin(2) + 1;
      
      sx = pixelX - floor(pixelX);
      sy = pixelY - floor(pixelY);
      
      particle_was_skipped = false;
      if  ( x_start > 0 && y_start > 0 && x_start + tileSize(1) - 1 < sTX && y_start + tileSize(2) - 1 < sTY )
        output_particle_stack(:,:,iGpuDataCounter) = STACK(x_start:x_start+tileSize(1)-1,y_start:y_start+tileSize(2)-1,TLT(iPrj,1));
      
        % The trasformation of the particle is e1,e2,e3,esym  into it's postion in the tomogram frame, then
        % the tomogram is tilted about the original Y axis and then the original Z
        % The angles stored are those used for interpolation, i.e. produced from BH_define_matrix([e1, e2, e3], 'Bah', 'inv' (or 'forwardVector'))
        % The angles are flipped in order so that the rotation matrix is R3*R2*R1 (really they should be flipped and negated, so the abvoe should be -e1, -e2, -e3, -esym)
        if (useFixedNotAliStack)
          rTilt = BH_defineMatrix([0,wrkDefAngTilt(iFid,3),wrkDefAngTilt(iFid,2)],'SPIDER','fwdVector');
          RF = fullXform(TLT(iPrj,1),1:4);
          rotFull = rTilt*[RF(1), RF(2), 0; RF(3), RF(4), 0; 0, 0, 1]*reshape(wrkPar(iFid,7:15),3,3);
        else
          % This gives us Rz*Ry
            rTilt = BH_defineMatrix([0,wrkDefAngTilt(iFid,3),wrkDefAngTilt(iFid,2)],'SPIDER','fwdVector');
          
          % this fives Rz*Ry*e3*e2*e1 * interpolant would rotate the particle by Rz*Ry*e1*e2*e3
          rotFull = rTilt*reshape(wrkPar(iFid,7:15),3,3);
        end
        
        eul = rotm2eul(rotFull,'ZYZ');
        e1 = 180./pi.*eul(1);
        
        e2 = 180./pi.*eul(2);
        e3 = 180./pi.*eul(3);
        
        
        phaseShift = 0.0;
        occupancy = 100.0; % TODO test replacement with CCC score?
        logp = -1000; % Tru -40000
        sigma = 10.0; % try 20
        score = 10.0; % TODO test with scaled CCC score?
        scoreChange = 0.0;
        pixelSize = emc.pixel_size_angstroms;
        micVoltage = emc.('VOLTAGE') * 10^-3;
        micCS = emc.('Cs') * 10^3;
        ampContrast = emc.('AMPCONT') * 10^0;
        beamTiltX = 0.0;
        beamTiltY = 0.0;
        beamTiltShiftX = 0.0;
        beamTiltShiftY = 0.0;
        best2dClass = 0.0;
        if (particle_was_skipped)
          beamTiltGroup = 0; % FSC half set, coopting this param for now.
        else
          beamTiltGroup = wrkPar(iFid,18); % FSC half set, coopting this param for now.
        end
        particleGroup = wrkPar(iFid,3);
        preExposure = wrkPar(iFid,16);
        totalExposure = wrkPar(iFid,17);
        pixelMultiplier = 0;
        xShift = pixelMultiplier*sx*pixelSize;
        yShift = pixelMultiplier*sy*pixelSize;

        % df1 = ( wrkPar(iFid,4) + wrkPar(iFid,5)) * 10;
        % df2 = ( wrkPar(iFid,4) - wrkPar(iFid,5)) * 10;
        % dfA = wrkPar(iFid,6)
        % fidIDX, tomoIdx, positionList(iSubTomo,4), d1, d2, 180./pi.*TLT(iPrj_nat,13), reshape(subtomo_rot_matrix,1,9), preExposure(iPrj_nat), postExposure(iPrj_nat), positionList(iSubTomo,7));

        df1 = (wrkDefAngTilt(iFid,1) + wrkPar(iFid,5)) * 10;
        df2 = (wrkDefAngTilt(iFid,1) - wrkPar(iFid,5)) * 10;
        dfA = wrkPar(iFid,6);
        
        fprintf(starFile, '%8u %7.2f %7.2f %7.2f %9.2f %9.2f %8.1f %8.1f %7.2f %7.2f %5i %7.2f %9i %10.4f %7.2f %8.5f %7.2f %7.2f %7.4f %7.3f %7.3f %7.3f %7.3f %5i %5i %8u %7.2f %7.2f\n', ...
          iDataCounter,-e1,-e2,-e3,xShift,yShift, ...
          df1,df2,dfA, ...
          phaseShift, occupancy, logp, sigma, score, scoreChange, ...
          pixelSize, micVoltage, micCS, ampContrast, ...
          beamTiltX, beamTiltY, beamTiltShiftX, beamTiltShiftY, ...
          best2dClass, beamTiltGroup, particleGroup, preExposure, totalExposure);
        
        
        iDataCounter = iDataCounter + 1;
        iGpuDataCounter = iGpuDataCounter + 1;
      end % if on windowing
    end % end of fiducial loop
    
  end % end of prj loop
  
  % Trim the stack to account for windowing skips
  output_particle_stack = gather(output_particle_stack(:,:,1:iGpuDataCounter - 1));
  tmp_stack_filename = sprintf('%s/%s_%d.mrc',mbOUT{1:2},iCell);
  SAVE_IMG(output_particle_stack, tmp_stack_filename, pixelSize);

  fprintf(newstack_file_handle, '%s\n',tmp_stack_filename);
  fprintf(newstack_file_handle, '0-%d\n',iGpuDataCounter-2);

  % output_cell{iCell}= gather(output_particle_stack);
  iCell = iCell + 1;

end % end of the loop over tilt series

fclose(newstack_file_handle);
fclose(starFile);

newstack_file_with_n_stacks = sprintf('%s/%s.newstack_full',mbOUT{1:2});
fh = fopen(newstack_file_with_n_stacks,'w');
fprintf(fh,'%d\n', iCell);
fclose(fh);

system(sprintf('cat %s >> %s',newstack_file,newstack_file_with_n_stacks));
system(sprintf('newstack -FileOfInputs %s %s.mrc > /dev/null',newstack_file_with_n_stacks,output_prefix));

% SAVE_IMG(cat(3,output_cell{:}),sprintf('%s.mrc',output_prefix),pixelSize);

maxThreads = emc.('nCpuCores');

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial reconstruction
%%%%%%%%%%%%%%%%%%%%%%%%%%%
system(sprintf('rm -f %s_rec.sh',output_prefix));
recScript = fopen(sprintf('%s_rec.sh',output_prefix), 'w');
fprintf(recScript,[ ...
  '#!/bin/bash\n\n', ...
  '%s << eof\n', ...
  '%s.mrc\n', ... sprintf('%s.mrc',output_prefix)
  '%s.star\n', ... sprintf('%s.star',output_prefix)
  'none.mrc\n', ...
  '%s_rec1.mrc\n',...
  '%s_rec2.mrc\n',...
  '%s_recFilt.mrc\n',...
  '%s_stats.txt\n',...
  '%s\n', ...
  '1\n', ...
  '0\n', ...
  '%3.3f\n', ... pixel size
  '%4.4f\n', ... molecularMass'
  '%3.3f\n', ... inermask ang
  '%3.3f\n', ... outermas ang
  '0.0\n', ... rec res limit
  '0.0\n', ... ref res limit
  '5.0\n', ... Particle weighting factor (A^2) [5.0]
  '1.0\n', ... Score threshold (<= 1 = percentage) [1.0]
  '1.0\n', ...Tuning parameter: smoothing factor [1.0]           :
  '1.0\n', ...Tuning parameters: padding factor [1.0]            :
  'Yes\n', ...Normalize particles [Yes]                          :
  'No\n', ...Adjust scores for defocus dependence [no]          :
  'No\n', ...Invert particle contrast [No]                      :
  'Yes\n', ...Exclude images with blank edges [yes]              :
  'No\n', ...Crop particle images [no]                          :
  'Yes\n', ...FSC calculation with even/odd particles [Yes]      :
  'No\n', ...Center mass [No]                                   :
  'No\n', ...Apply likelihood blurring [No]                     :
  'No\n', ...Threshold input reconstruction [No]                :
  'No\n', ...Dump intermediate arrays (merge later) [No]        :
  'dum_1.dat\n', ...Output dump filename for odd particle [dump_file_1.dat]                                  :
  'dum_2.dat\n', ...Output dump filename for even particle [dump_file_2.dat]                                  :
  '%2.2d\n', ...Max. threads to use for calculation [36]           :
  ], getenv('EMC_RECONSTRUCT3D'),output_prefix, output_prefix, output_prefix, output_prefix, output_prefix, output_prefix, ...
  symmetry,emc.pixel_size_angstroms, ...
  emc.('particleMass')*10^3, 0.0, mean(emc.('Ali_mRadius')), maxThreads);

fprintf(recScript, '\neof\n');

fclose(recScript);
system(sprintf('chmod a=wrx %s_rec.sh',output_prefix));
system(sprintf('./%s_rec.sh',output_prefix));


%%%%%%%%%%%%%%%%%%%%%%%%%
% Refine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
system(sprintf('rm -f %s_ref.sh',output_prefix));
refineScript = fopen(sprintf('%s_ref.sh',output_prefix), 'w');
fprintf(refineScript,[ ...
  '#!/bin/bash\n\n', ...
  '%s << eof\n', ...
  '%s.mrc\n', ... sprintf('%s.mrc',output_prefix)
  '%s.star\n', ... sprintf('%s.star',output_prefix)
  '%s_recFilt.mrc\n',...
  '%s_stats.txt\n',...
  'yes\n',... Use statistics [Yes]                               :
  'my_projection_stack.mrc\n',... not going to be used                       :
  '%s_refined.star\n', ...
  '%s_changes.star\n', ...Output parameter changes
  '%s\n',... Particle symmetry [C1]                             :
  '1\n', ...First particle to refine (0 = first in stack) [1]  :
  '0\n', ...Last particle to refine (0 = last in stack) [0]    :
  '1.0\n',...Percent of particles to use (1 = all) [1.0]        :
  '%3.3f\n', ... pixel size
  '%4.4f\n', ... molecularMass'
  '%3.3f\n', ... inermask ang
  '%3.3f\n', ... outermas ang
  '300.0\n',...Low resolution limit (A) [300.0]                   :
  '%3.3f\n',...High resolution limit (A) [8.0]                    :
  '0.0\n',...Resolution limit for signed CC (A) (0.0 = max [0.0]                                              :
  '0.0\n',...Res limit for classification (A) (0.0 = max) [0.0] :
  '0.0\n',...Mask radius for global search (A) (0.0 = max)[100.0]                                            :
  '%3.3f\n',...Approx. resolution limit for search (A) [8]        :
  '0.0\n',...Angular step (0.0 = set automatically) [0.0]       :
  '20\n',...Number of top hits to refine [20]                  :
  '10\n',...Search range in X (A) (0.0 = 0.5 * mask radius)[12]                                               :
  '10\n',...[12]                                               :
  '100.0\n',...2D mask X coordinate (A) [100.0]                   :
  '100.0\n',...2D mask Y coordinate (A) [100.0]                   :
  '100.0\n',...2D mask Z coordinate (A) [100.0]                   :
  '100.0\n',...2D mask radius (A) [100.0]                         :
  '500.0\n',...Defocus search range (A) [500.0]                   :
  '50.0\n',...Defocus step (A) [50.0]                            :
  '1.0\n',...Tuning parameters: padding factor [1.0]            :
  'no\n',...Global search [No]                                 :
  'yes\n',...  Local refinement [Yes]                             :
  'no\n',...Refine Psi [no]                                    :
  'no\n',...Refine Theta [no]                                  :
  'no\n',...Refine Phi [no]                                    :
  'yes\n',...Refine ShiftX [Yes]                                :
  'yes\n',...Refine ShiftY [Yes]                                :
  'no\n',...Calculate matching projections [No]                :
  'no\n',...Apply 2D masking [No]                              :
  'no\n',...Refine defocus [No]                                :
  'yes\n',...Normalize particles [Yes]                          :
  'no\n',...Invert particle contrast [No]                      :
  'yes\n',...Exclude images with blank edges [Yes]              :
  'yes\n',...Normalize input reconstruction [Yes]               :
  'no\n',...Threshold input reconstruction [No]                :
  '%2.2d\n', ...Max. threads to use for calculation [36]           :
  ],  getenv('EMC_REFINE3D'),output_prefix, output_prefix, output_prefix, output_prefix, output_prefix, output_prefix, ...
  symmetry,emc.pixel_size_angstroms, ...
  emc.('particleMass')*10^3, 0.0, mean(emc.('Ali_mRadius')), ...
  resForFitting,resForFitting,maxThreads);

fprintf(refineScript, '\neof\n');
fclose(refineScript);
pause(3);
system(sprintf('chmod a=wrx %s_ref.sh',output_prefix));
system(sprintf('./%s_ref.sh',output_prefix));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reconstruct refined
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

system(sprintf('rm -f %s_rec2.sh',output_prefix));
recScript = fopen(sprintf('%s_rec2.sh',output_prefix), 'w');
fprintf(recScript,[ ...
  '#!/bin/bash\n\n', ...
  '%s << eof\n', ...
  '%s.mrc\n', ... sprintf('%s.mrc',output_prefix)
  '%s_refined.star\n', ... sprintf('%s.star',output_prefix)
  'none.mrc\n', ...
  '%s_rec1.mrc\n',...
  '%s_rec2.mrc\n',...
  '%s_recFilt_refined.mrc\n',...
  '%s_stats_refined.txt\n',...
  '%s\n', ...
  '1\n', ...
  '0\n', ...
  '%3.3f\n', ... pixel size
  '%4.4f\n', ... molecularMass'
  '%3.3f\n', ... inermask ang
  '%3.3f\n', ... outermas ang
  '0.0\n', ... rec res limit
  '0.0\n', ... ref res limit
  '5.0\n', ... Particle weighting factor (A^2) [5.0]
  '1.0\n', ... Score threshold (<= 1 = percentage) [1.0]
  '1.0\n', ...Tuning parameter: smoothing factor [1.0]           :
  '1.0\n', ...Tuning parameters: padding factor [1.0]            :
  'Yes\n', ...Normalize particles [Yes]                          :
  'No\n', ...Adjust scores for defocus dependence [no]          :
  'No\n', ...Invert particle contrast [No]                      :
  'Yes\n', ...Exclude images with blank edges [yes]              :
  'No\n', ...Crop particle images [no]                          :
  'Yes\n', ...FSC calculation with even/odd particles [Yes]      :
  'No\n', ...Center mass [No]                                   :
  'No\n', ...Apply likelihood blurring [No]                     :
  'No\n', ...Threshold input reconstruction [No]                :
  'No\n', ...Dump intermediate arrays (merge later) [No]        :
  'dum_1.dat\n', ...Output dump filename for odd particle [dump_file_1.dat]                                  :
  'dum_2.dat\n', ...Output dump filename for even particle [dump_file_2.dat]                                  :
  '%2.2d\n', ...Max. threads to use for calculation [36]           :
  ], getenv('EMC_RECONSTRUCT3D'), output_prefix, output_prefix, output_prefix, output_prefix, output_prefix, output_prefix, ...
  symmetry,emc.pixel_size_angstroms, ...
  emc.('particleMass')*10^3, 0.0, mean(emc.('Ali_mRadius')), maxThreads);

fprintf(recScript, '\neof\n');

fclose(recScript);
pause(2)
system(sprintf('chmod a=wrx %s_rec2.sh',output_prefix));
pause(2)
system(sprintf('./%s_rec2.sh',output_prefix));


end

