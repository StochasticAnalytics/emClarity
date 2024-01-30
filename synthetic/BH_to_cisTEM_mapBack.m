function [ ] = BH_to_cisTEM_mapBack(PARAMETER_FILE, CYCLE,outputName, symmetry, MAX_EXPOSURE, varargin)

% Map back and align using the subtomograms as fiducial markers.

% If multiple classes are being mapped back, pass in the className (refName
% really b/c these will be the only re-weighted volumes) and then each class
% will be given a unique density value in a seperate volume used to visualize
% color in Chimera.
%
% Otherwise, pass a string that points at a single volume to use.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some flags that are worth keeping as options, but not accessible directly


global bh_global_imodProjectionShifts;
if isempty(bh_global_imodProjectionShifts)
  %  bh_global_imodProjectionShifts = [ -0.5, -0.5, 0.5 ; -0.5, -0.5, 0; 0.5,0.5,1.0 ];
  % From a more thorough sweep
  bh_global_imodProjectionShifts = [ 0.5, -0.5, 0.5 ; 0.0, -0.5, 0; 0.5,0.5,1.0 ];
  
end


tiltStart=1;
MAX_EXPOSURE = EMC_str2double(MAX_EXPOSURE)
if isnan(MAX_EXPOSURE)
  error('MAX_EXPOSURE is nan - if running from an interactive matlab session, did you enter as a string?');
end
CYCLE = EMC_str2double(CYCLE);

if CYCLE < 0
  CYCLE = abs(CYCLE);
  doFullXform = true;
else
  doFullXform = false;
end

pixelShift = -1;
pixelMultiplier = 0;

cacheAdd = '';
if nargin > 5
  preShift = varargin{1};
  postShift = varargin{2};
  prjVectorShift = varargin{3}';
  if length(varargin) > 3
    cacheAdd = varargin{4};
  else
    cacheAdd = '';
  end
else
  
  % in tilt the coords are transformed from the model coordinate system to
  % the tomogram by [0.5,0.5,1.0]. After transformation, these are then
  % subtracted back off.
  preShift = bh_global_imodProjectionShifts(1,:);
  postShift = bh_global_imodProjectionShifts(2,1:2);
  prjVectorShift = bh_global_imodProjectionShifts(3,:)';
  %   prjVectorShift = [0,0,-1];
  
  
end
% baseFile = sprintf('%s_%d_%2.2f_preShift_%2.2f_%2.2f_%2.2f_postShift_%2.2f_%2.2f_prjVect_%2.2f_%2.2f_%2.2f','microShiftsFollowup',MAX_EXPOSURE, preShift, postShift, prjVectorShift);
baseFile = outputName; %sprintf('%s_%d_%2.2f','withZeroedXF_IPFirst',MAX_EXPOSURE);

useFixedNotAliStack = false;

cycleNumber = sprintf('cycle%0.3u', CYCLE);

emc = BH_parseParameterFile(PARAMETER_FILE);
reconScaling = 1;
samplingRate = 1; % Always working at full binning. emc.('Ali_samplingRate');

load(sprintf('%s.mat', emc.('subTomoMeta')), 'subTomoMeta');
resForFitting = 1.3*mean(subTomoMeta.currentResForDefocusError);

% % Add error check onrange for reasonable values.
% ctfRange = emc.('tomoCprDefocusRange')*10^10;
% ctfInc = emc.('tomoCprDefocusStep')*10^10;

% calcCTF = emc.('tomoCprDefocusRefine');

nGPUs = emc.('nGPUs');
pInfo = parcluster();
gpuScale=3*samplingRate
nWorkers = min(nGPUs*gpuScale,emc.('nCpuCores')); % 18
fprintf('Using %d workers as max of %d %d*nGPUs and %d nWorkers visible\n', ...
  nWorkers,gpuScale,nGPUs*gpuScale,pInfo.NumWorkers);

tmpCache = sprintf('cache/to_cisTEM%s/',cacheAdd);
CWD = '';

system(sprintf('mkdir -p %s',tmpCache));


load(sprintf('%s.mat', emc.('subTomoMeta')), 'subTomoMeta');
mapBackIter = subTomoMeta.currentTomoCPR;


system(sprintf('mkdir -p %smapBack%d',tmpCache, mapBackIter+1));


tiltNameList = fieldnames(subTomoMeta.mapBackGeometry);

tiltNameList = tiltNameList(~ismember(tiltNameList,{'tomoName','viewGroups'}));
nTiltSeries = length(tiltNameList);

% Cycle 0 is named differently - I'll be deleting this in an overhaul of the way
% the subTomoMeta is written.
if (CYCLE)
  try
    geometry = subTomoMeta.(cycleNumber).RawAlign;
    fprintf('Using Alignment geometry %s\n',cycleNumber);
  catch
    geometry = subTomoMeta.(cycleNumber).Avg_geometry;
    fprintf('Using Average geometry %s\n',cycleNumber);
  end
else
  try
    geometry = subTomoMeta.(cycleNumber).RawAlign;
    fprintf('Using Alignment geometry %s\n',cycleNumber);
    
  catch
    geometry = subTomoMeta.(cycleNumber).geometry;
    fprintf('Using Average geometry %s\n',cycleNumber);
  end
end

tiltGeometry = subTomoMeta.tiltGeometry;

% Assume No 2d CTF until proven otherwise

firstTilt = 1;
nFidsTotalDataSet= 0;
iCell = 1;
output_cell = {};
% TODO split this up into chunks
for iTiltSeries = tiltStart:nTiltSeries
  
  nTomograms = subTomoMeta.mapBackGeometry.(tiltNameList{iTiltSeries}).nTomos;
  if nTomograms == 0
    % No points were saved after template matching so skip this tilt seoarries
    % altogether.
    continue
  end
  
  tomoList = {};
  tomoIDX = 1;
  for iTomo = 1:size(subTomoMeta.mapBackGeometry.(tiltNameList{iTiltSeries}).coords,1)
    % This is dumb, fix it to be explicit.
    if any(subTomoMeta.mapBackGeometry.(tiltNameList{iTiltSeries}).coords(iTomo,:))
      tomoList{tomoIDX} = sprintf('%s_%d',tiltNameList{iTiltSeries},iTomo);
      
      tiltList{tomoIDX} = sprintf('%saliStacks/%s_ali%d.fixed',...
        CWD,tiltNameList{iTiltSeries},mapBackIter+1);
      % Only increment if values found.
      tomoIDX = tomoIDX + 1;
    end
    
  end
  
  if (mapBackIter)
    localFile = sprintf('%smapBack%d/%s_ali%d_ctf.local', ...
      CWD,mapBackIter,tiltNameList{iTiltSeries},mapBackIter);
  else
    localFile = sprintf('%sfixedStacks/%s.local',CWD,tiltNameList{iTiltSeries});
  end
  
  if exist(localFile,'file')
    fprintf('Found local file %s\n.', localFile);
  else
    fprintf('No local transforms found.\n');
    localFile = 0;
  end
  
  
  % The model is scaled to full sampling prior to passing to tiltalign,
  % make sure the header in the synthetic stack is set appropriately.
  fullPixelSize = emc.pixel_size_angstroms;
  pixelSize = fullPixelSize.*samplingRate;
  
  PARTICLE_RADIUS = floor(max(emc.('particleRadius')./pixelSize));
  
  [~,tiltBaseName,~] = fileparts(tiltList{1});
  mbOUT = {[tmpCache],[mapBackIter+1],[tiltBaseName]};
  fprintf('\nmBOUT name is %smapBack%d/%s\n',mbOUT{1:3});
  
  
  tiltHeader = getHeader(MRCImage(tiltList{1},0));
  % This is only needed in the re-projection of the model. I don't think it
  % should affect anything, but double check. FIXME
  maxZ = 100;
  
  
  % The tilt angles are the same for each tomo, so it is okay to just use
  % number 1 here.
  TLT = tiltGeometry.(tomoList{1});
  
  
  iRawTltName = sprintf('%smapBack%d/%s_align.rawtlt',mbOUT{1:3})
  iTiltFile = fopen(iRawTltName, 'w');
  rawTLT = sortrows(TLT(:,[1,4]),1);
  fprintf(iTiltFile,'%f\n',rawTLT(:,2)');
  fclose(iTiltFile);
  
  coordOUT = fopen(sprintf('%smapBack%d/%s.coord',mbOUT{1:3}),'w');
  coordSTART = fopen(sprintf('%smapBack%d/%s.coord_start',mbOUT{1:3}),'w');
  
  defOUT   = fopen(sprintf('%smapBack%d/%s.defAng',mbOUT{1:3}),'w');
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
    
    
    
    if (useFixedNotAliStack || doFullXform)
      
      % 20190509 - I think this is royally screwing things up FIXME
      % Commenting this out invalidates the defocus vals
      xfTLT = sortrows(TLT(:,[1,7:10,2,3],1));
      fprintf(iXF,'%f %f %f %f %f %f\n',xfTLT(:,2:7)');
      fclose(iXF);
      % Odd size stacks are enforced which creates a shift prior to the
      % xform.
      if (useFixedNotAliStack)
        isEven = 1;
        
        iXFBase = sprintf('%smapBack%d/%s_align_base.XF',mbOUT{1:3});
        iXFB = fopen(iXFBase,'w');
        for ix = 1:size(xfTLT,1)
          fprintf(iXFB,'%f %f %f %f %f %f\n',[1,0,0,1,-isEven,-isEven]);
        end
        
        fclose(iXFB);
        system(sprintf('xfproduct %s %s %s',iXFBase, iXFName,iXFName));
        
        % We need to invert this transform to map from the aligned stack to the
        % fixed stack
        iXFName_inv = sprintf('%smapBack%d/%s_align_inv.XF',mbOUT{1:3});
        system(sprintf('xfinverse %s %s', iXFName, iXFName_inv));
      end
    else
      % 20190509 - I think this is royally screwing things up FIXME
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
    
    positionList = positionList(positionList(:,26) ~= -9999,:);
    nFidsTotal = nFidsTotal + size(positionList,1);
    
    sTX = floor(tiltHeader.nX );
    sTY = floor(tiltHeader.nY );
    
    originPrj = floor([sTX,sTY,1]./2) + 1;
    
    tomoNumber = subTomoMeta.mapBackGeometry.tomoName.(tomoList{iTomo}).tomoNumber;
    tiltName = subTomoMeta.mapBackGeometry.tomoName.(tomoList{iTomo}).tiltName;
    reconCoords = subTomoMeta.mapBackGeometry.(tiltName).coords(tomoNumber,:);
    
    flgLoad = 0;
    [~,tomoReconCoords] = BH_multi_loadOrBuild(tomoList{iTomo}, ...
      reconCoords, mapBackIter, ...
      samplingRate, 1,reconScaling,flgLoad, 'tomoCPR');
    
    originVol = floor(tomoReconCoords(1,1:3)./2) + 1;
    reconShift = tomoReconCoords(2,1:3);
    
    maxZ = 1000; % Does not seem to affect anything.
    reconstructionSize = [tiltHeader.nX,tiltHeader.nY,maxZ];
    originRec = floor(reconstructionSize./2) + 1;
    
    
    nPrjs = size(TLT,1);
    nSubTomos = size(positionList,1);
    
    
    if (iTomo == 1)
      fidIDX = 0;
    end
    
    modelRot = BH_defineMatrix([0,90,0],'Bah','forwardVector');
    
    for iSubTomo = 1:nSubTomos
      
      
      rSubTomo = reshape(positionList(iSubTomo,17:25),3,3);
      prjVector = (positionList(iSubTomo,11:13)./samplingRate) - originVol + reconShift;
      
      
      % % %       nRefs = 1;
      % % %       if (nRefs > 1)
      % % %         iClassIDX = positionList(iSubTomo,26);
      % % %       else
      % % %         iClassIDX = 1;
      % % %       end
      
      %       prjVector = prjVector - [0.5,0.5,1.0]; %prjVectorShift;
      prjVector = prjVector - preShift;
      % Reproject using tilt, so just save the 3d coords.
      fprintf(coordOUT,'%0.4f %0.4f %0.4f %d\n', modelRot*prjVector' + [originRec(1),originRec(3),originRec(2)]' - prjVectorShift([1,3,2]), fidIDX);
      
      nPrjsIncluded = 0;
      for iPrj = 1:nPrjs
        
        iPrj_nat = find(TLT(:,1) == iPrj);
        if (abs(TLT(iPrj_nat,11)) <= MAX_EXPOSURE)
          nPrjsIncluded = nPrjsIncluded + 1;
          
          
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
          
          nFidsTotalDataSet = nFidsTotalDataSet + 1;
          
        else
          fprintf(coordSTART,'%d %d %d %3.3f %3.3f %3.3f %3.3f %3.3f %3.3f %3.3f %3.3f %3.3f %3.3f %3.3f %3.3f %3.3f %3.3f %d\n',-9999, -9999,-9999,1.0,1.0,1.0,1,1,1,1,1,1,1,1,1,0,0,1);
          
          nFidsTotalDataSet = nFidsTotalDataSet + 1;
        end
      end % loop over tilt projections
      
      
      fidIDX = fidIDX + 1;
      
    end % loop over subtomos
    
    
  end
  
  
  fclose(coordOUT);
  fclose(coordSTART);
  
  
  p2m = sprintf(['point2model -zero -circle 3 -color 0,0,255 -values -1 ',...
    '%smapBack%d/%s.coord %smapBack%d/%s.3dfid > /dev/null'], ...
    mbOUT{1:3},mbOUT{1:3});
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
  reModFileName = sprintf('%smapBack%d/%s_%d_reMod.sh',mbOUT{1:3},iSave);
  reModFile = fopen(reModFileName,'w');
  invertTiltAngles = 0;
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
    pixelSize/10, invertTiltAngles,... % Ang --> nm
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
  
  if (useFixedNotAliStack)
    % transform the projected model back to the fixed stack frame, and
    % then convert to text.
    
    system(sprintf('imodtrans -2 %s %smapBack%d/%s.fid %smapBack%d/%s.invfid', iXFName_inv, mbOUT{1:3},mbOUT{1:3}));
    
    system(sprintf(['model2point -contour -zero ',...
      '%smapBack%d/%s.invfid %smapBack%d/%s.coordPrj'],...
      mbOUT{1:3}, mbOUT{1:3}))
  else
    system(sprintf(['model2point  -contour -zero ',...
      '%smapBack%d/%s.fid %smapBack%d/%s.coordPrj'],...
      mbOUT{1:3}, mbOUT{1:3}))
  end
  
  
  
  try
    fidList = load(sprintf('%smapBack%d/%s.coordPrj',mbOUT{1:3}));
  catch
    fprintf('\nWarning, did not load the projected coords\nSkipping along');
    continue;
  end
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
  
  fidList = [1:size(fidList,1);fidList']';
  
  
  
  particlePad = 2.0;
  tileRadius = floor(particlePad.*PARTICLE_RADIUS);
  tileSize = (2.*tileRadius).*[1,1];
  
  tileOrigin = floor(tileSize./2) + 1;
  
  nFidsTotal = numel(unique(fidList(parList(:,1)~=-9999,2)));
  
  
  if (useFixedNotAliStack)
    tiltSeries = sprintf('%sfixedStacks/%s.fixed',CWD,tiltName);
  else
    tiltSeries = sprintf('%saliStacks/%s_ali%d.fixed',CWD,tiltName,mapBackIter+1);
  end
  
  % This will need to be changed to aggregate
  output_particle_stack = zeros([tileSize,nFidsTotal*nPrjsIncluded],'single');
  
  iGpuDataCounter = 1;
  
  if (firstTilt)
    
    iDataCounter = 1;
    starFile = fopen(sprintf('%s.star',baseFile),'w');
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
    
    
    firstTilt = 0;
  end
  
  if (useFixedNotAliStack)
    fullXform = load(iXFName_inv);
  end
  
  STACK = single(getVolume(MRCImage(tiltSeries)));
  
  for iPrj = 1:nPrjs
    
    if (abs(TLT(iPrj,11)) > MAX_EXPOSURE)
      continue;
    end
    
    dataPrj = STACK(:,:,TLT(iPrj,1));
    
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
      
      
      pixelX = wrkFid(iFid,3) - pixelShift + postShift(1);
      pixelY = wrkFid(iFid,4) - pixelShift + postShift(2);
      
      ox = floor(pixelX) - tileRadius;
      oy = floor(pixelY) - tileRadius;
      
      
      sx = pixelX - floor(pixelX);
      sy = pixelY - floor(pixelY);
      
      particle_was_skipped = false;
      if  ( ox > 0 && oy > 0 && ox + 2*tileRadius < sTX && oy +2*tileRadius < sTY )
        output_particle_stack(:,:,iGpuDataCounter) = dataPrj(ox:ox+2.*tileRadius-1,oy:oy+2.*tileRadius-1);
      else
        particle_was_skipped = true;
        output_particle_stack(:,:,iGpuDataCounter) = randn(tileSize,'single').*0.1; % Why am I not just skipping these?
      end
      
      if (useFixedNotAliStack)
        rTilt = BH_defineMatrix([0,wrkDefAngTilt(iFid,3),wrkDefAngTilt(iFid,2)],'SPIDER','forwardVector');
        RF = fullXform(TLT(iPrj,1),1:4);
        rotFull = rTilt*[RF(1), RF(2), 0; RF(3), RF(4), 0; 0, 0, 1]*reshape(wrkPar(iFid,7:15),3,3);
      else
        
        
        %           rTilt = BH_defineMatrix([0,wrkDefAngTilt(iFid,3),wrkDefAngTilt(iFid,2)],'SPIDER','forwardVector');
        rTilt = BH_defineMatrix([wrkDefAngTilt(iFid,2),wrkDefAngTilt(iFid,3),0],'SPIDER','forwardVector');
        
        rotFull = rTilt*reshape(wrkPar(iFid,7:15),3,3);
      end
      
      eul = rotm2eul(rotFull,'ZYZ');
      e1 = 180./pi.*eul(1);
      
      e2 = 180./pi.*eul(2);
      e3 = 180./pi.*eul(3);
      
      
      phaseShift = 0.0;
      occupancy = 100.0; % TODO test replacement with CCC score?
      logp = -1000;
      sigma = 10.0;
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
      
      xShift = pixelMultiplier*sx*pixelSize;
      yShift = pixelMultiplier*sy*pixelSize;
      
      %         df1 = wrkPar(iFid,4);
      %         df2 = wrkPar(iFid,5);
      %         dfA = wrkPar(iFid,6);
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
      
    end
    
  end % end of prj loop
  
  
  output_cell{iCell}= gather(output_particle_stack);
  iCell = iCell + 1;
  
end

fclose(starFile);


SAVE_IMG(cat(3,output_cell{:}),sprintf('%s.mrc',baseFile),pixelSize);

maxThreads = emc.('nCpuCores');

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial reconstruction
%%%%%%%%%%%%%%%%%%%%%%%%%%%
system(sprintf('rm -f %s_rec.sh',baseFile));
recScript = fopen(sprintf('%s_rec.sh',baseFile), 'w');
fprintf(recScript,[ ...
  '#!/bin/bash\n\n', ...
  '%s << eof\n', ...
  '%s.mrc\n', ... sprintf('%s.mrc',baseFile)
  '%s.star\n', ... sprintf('%s.star',baseFile)
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
  ], getenv('EMC_RECONSTRUCT3D'),baseFile, baseFile, baseFile, baseFile, baseFile, baseFile, ...
  symmetry,emc.pixel_size_angstroms, ...
  emc.('particleMass')*10^3, 0.0, mean(emc.('Ali_mRadius')), maxThreads);

fprintf(recScript, '\neof\n');

fclose(recScript);
system(sprintf('chmod a=wrx %s_rec.sh',baseFile));
pause(3)
system(sprintf('./%s_rec.sh',baseFile));


%%%%%%%%%%%%%%%%%%%%%%%%%
% Refine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
system(sprintf('rm -f %s_ref.sh',baseFile));
refineScript = fopen(sprintf('%s_ref.sh',baseFile), 'w');
fprintf(refineScript,[ ...
  '#!/bin/bash\n\n', ...
  '%s << eof\n', ...
  '%s.mrc\n', ... sprintf('%s.mrc',baseFile)
  '%s.star\n', ... sprintf('%s.star',baseFile)
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
  ],  getenv('EMC_REFINE3D'),baseFile, baseFile, baseFile, baseFile, baseFile, baseFile, ...
  symmetry,emc.pixel_size_angstroms, ...
  emc.('particleMass')*10^3, 0.0, mean(emc.('Ali_mRadius')), ...
  resForFitting,resForFitting,maxThreads);

fprintf(refineScript, '\neof\n');
fclose(refineScript);
pause(3);
system(sprintf('chmod a=wrx %s_ref.sh',baseFile));
system(sprintf('./%s_ref.sh',baseFile));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reconstruct refined
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

system(sprintf('rm -f %s_rec2.sh',baseFile));
recScript = fopen(sprintf('%s_rec2.sh',baseFile), 'w');
fprintf(recScript,[ ...
  '#!/bin/bash\n\n', ...
  '%s << eof\n', ...
  '%s.mrc\n', ... sprintf('%s.mrc',baseFile)
  '%s_refined.star\n', ... sprintf('%s.star',baseFile)
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
  ], getenv('EMC_RECONSTRUCT3D'), baseFile, baseFile, baseFile, baseFile, baseFile, baseFile, ...
  symmetry,emc.pixel_size_angstroms, ...
  emc.('particleMass')*10^3, 0.0, mean(emc.('Ali_mRadius')), maxThreads);

fprintf(recScript, '\neof\n');

fclose(recScript);
pause(2)
system(sprintf('chmod a=wrx %s_rec2.sh',baseFile));
pause(2)
system(sprintf('./%s_rec2.sh',baseFile));


end

