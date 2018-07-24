function [ ] = BH_synthetic_mapBack(PARAMETER_FILE, CYCLE, tiltStart)

% Map back and align using the subtomograms as fiducial markers.

% If multiple classes are being mapped back, pass in the className (refName
% really b/c these will be the only re-weighted volumes) and then each class
% will be given a unique density value in a seperate volume used to visualize
% color in Chimera.
%
% Otherwise, pass a string that points at a single volume to use.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some flags that are worth keeping as options, but not accessible directly
% by the users (private methods-ish)
buildTomo=1;% % % % % % %
METHOD = 'GPU';
 flgRunAlignments = true;
COLOR_MAP= '0';
flgClassAvg = 0;


% Color map to reproject for something like the ribosome.
if ~isempty(str2num(COLOR_MAP))
  flgColorMap = 0;
else
  flgColorMap = 1;
end
CYCLE = str2num(CYCLE);

if CYCLE < 0
  % Additional node, for now you'll have to manually copy the results.
  % Writes the runAlignments.sh to runAlignments_alt.sh which manually cat
  % on end of main run.
  % Does not save subTomoMeta
  flgAltRun = 1
  CYCLE = abs(CYCLE)
else
  flgAltRun = 0
end
cycleNumber = sprintf('cycle%0.3u', CYCLE);



pBH = BH_parseParameterFile(PARAMETER_FILE);
reconScaling = 1;
samplingRate = pBH.('Ali_samplingRate');
% used to determine the number of fiducials/patch for local area. 
MOL_MASS = pBH.('particleMass');
molMass = MOL_MASS.*(25/samplingRate); 
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
  flgInvertTiltAngles = 1;


try
  flgLocalMag = pBH.('flgLocalMag');
  fprintf('localMag from param %d\n',flgLocalMag);
catch
  flgLocalMag = 1;
end

try 
  tiltAliOption = pBH.('tiltAliOption');
catch
  % Global option, Global Grouping, Local Opt, Local Grouping.
  % All early tests were Opt 5 (linear mapping) Default Grouping 5
  tiltAliOption = [5,5,5,5];

  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tmpCache= pBH.('fastScratchDisk');

% % % nWorkers = str2num(nWORKERS)
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
ctfRange = pBH.('tomoCprDefocusRange'); 
ctfInc = pBH.('tomoCprDefocusStep');

calcCTF = pBH.('tomoCprDefocusRefine');

tiltNameList = fieldnames(subTomoMeta.mapBackGeometry);

tiltNameList = tiltNameList(~ismember(tiltNameList,{'tomoName','viewGroups'}));
nTiltSeries = length(tiltNameList)

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

% Load in the reference images.
try
  refNameODD = sprintf('%s_%s_class0_REF_ODD.mrc', ...
                                             cycleNumber,pBH.('subTomoMeta'));
  refNameEVE = sprintf('%s_%s_class0_REF_EVE.mrc', ...
                                             cycleNumber,pBH.('subTomoMeta'));     
  refODD = getVolume(MRCImage(refNameODD));
  refEVE = getVolume(MRCImage(refNameEVE));                                           
catch
  fprintf('\nDid not find either %s or %s, trying Raw prefix\n',refNameODD,refNameEVE);
  try
    refNameODD = sprintf('%s_%s_class0_Raw_ODD.mrc', ...
                                             cycleNumber,pBH.('subTomoMeta'));
    refNameEVE = sprintf('%s_%s_class0_Raw_EVE.mrc', ...
                                             cycleNumber,pBH.('subTomoMeta'));                                             

    refODD = getVolume(MRCImage(refNameODD));
    refEVE = getVolume(MRCImage(refNameEVE));
  catch
    error('\nDid not find either %s or %s\n',refNameODD,refNameEVE)
  end
end

try
  conserveDiskSpace = pBH.('conserveDiskSpace');
catch
  conserveDiskSpace = 0;
end
tiltGeometry = subTomoMeta.tiltGeometry;

% Assume No 2d CTF until proven otherwise
flg2dCTF = 0;
outCTF = '';
for iTiltSeries = tiltStart:nTiltSeries
    
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
  tomoIDX = 1;
  for iTomo = 1:size(subTomoMeta.mapBackGeometry.(tiltNameList{iTiltSeries}).coords,1)
    % This is dumb, fix it to be explicit.
    if any(subTomoMeta.mapBackGeometry.(tiltNameList{iTiltSeries}).coords(iTomo,:))
      tomoList{tomoIDX} = sprintf('%s_%d',tiltNameList{iTiltSeries},iTomo);
     
    
    
     if ( flg2dCTF )
        tiltList{tomoIDX} = sprintf('%sctfStacks/%s_ali%d_ctf.fixed', ...
                                CWD,tiltNameList{iTiltSeries},mapBackIter+1);

        [noHeader,~] = system(sprintf('header %s > /dev/null',tiltList{tomoIDX}));
        if ( noHeader ) 
          tiltList{tomoIDX} = sprintf('%saliStacks/%s_ali%d.fixed', ...
                                   CWD,tiltNameList{iTiltSeries},mapBackIter+1);
          flg2dCTF = 0;
          outCTF='_ctf';
        else
          flg2dCTF = 1;
        end
      else
        tiltList{tomoIDX} = sprintf('%saliStacks/%s_ali%d.fixed',...
                                 CWD,tiltNameList{iTiltSeries},mapBackIter+1);
        outCTF='_ctf';                       
     end
    % Only increment if values found.
    tomoIDX = tomoIDX + 1;
    end
    
  end

  fprintf('\nComparing against non-CTF corrected stack\n');
  
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
  eraseMaskType = pBH.('peak_mType');
	eraseMaskRadius = pBH.('peak_mRadius')./pixelSize;
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
  peakSearchRad = floor(0.2*PARTICLE_RADIUS.*[1,1]);
  try
    lowPassCutoff = pBH.('tomoCprLowPass');
    fprintf('Using a user supplied lowpass cutoff of %3.3f Ang\n.',...
            lowPassCutoff);
  catch
    lowPassCutoff = 1.5.*mean(subTomoMeta.currentResForDefocusError);
    if (lowPassCutoff < 10)
      lowPassCutoff = 10;
    elseif (lowPassCutoff > 24)
      lowPassCutoff = 24;
    end
    fprintf('Using an internatlly determined lowpass cutoff of %3.3f Ang\n.',...
            lowPassCutoff);
  end
  nFiducialsPerPatch = ceil(100./sqrt(molMass))
  targetPatchSize = max(500, ceil(2.*(PARTICLE_RADIUS).*sqrt(nFiducialsPerPatch)))

  [~,tiltBaseName,~] = fileparts(tiltList{1});
  mbOUT = {[tmpCache],[mapBackIter+1],[tiltBaseName]};
  fprintf('\nmBOUT name is %smapBack%d/%s\n',mbOUT{1:3});

  % Check to see if this tilt has already been worked on, if so skip
  aliCmdFileCheck = sprintf('%smapBack%d/%s.align',mbOUT{1:3});
  if exist(aliCmdFileCheck,'file')
    fprintf('\n\nFound aliCmdFileCheck, skipping rather than overwrite.\n');
    continue
  end
% % %   tiltList_orig = tiltList;
% % %   tomoList_bin = tomoList;
  % Everything is there, now make sure all binned data is there
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



  fprintf('\n\nflg2Dctf check %d\n\n',flg2dCTF);
% move this into subTomo meta and track failure to optmize per tilt as
% alignments proceed (bigger local shifts require more memory)
% % %   mapBackRePrjSize = pBH.('tomoCprRePrjSize')
  
  % re-initialize the parpool for each tilt series to free up mem.
  if ~isempty(gcp('nocreate'))
    delete(gcp)
    parpool(nWorkers);
  else
    parpool(nWorkers);
  end
  fprintf('init with %d workers\n',nWorkers);

  outputPrefix = sprintf('%s_%s', cycleNumber, pBH.('subTomoMeta'));





  

  % Get the thickest for recon
  maxZ = 0;
  tiltList{1}
  tiltHeader = getHeader(MRCImage(tiltList{1},0));

  for iTomo = 1:nTomograms

% 
%     reconGeometry.(tomoList{iTomo}) = floor(reconGeometry.(tomoList{iTomo}) ./ ...
%                                                               samplingRate)
    tomoNumber = subTomoMeta.mapBackGeometry.tomoName.(tomoList{iTomo}).tomoNumber
    tiltName = subTomoMeta.mapBackGeometry.tomoName.(tomoList{iTomo}).tiltName
    nZdZ = subTomoMeta.mapBackGeometry.(tiltName).coords(tomoNumber,[4,6])./samplingRate
    
    % half the size in z plus the shift back to the microscope coords.
    sZneeded = 2.*ceil(nZdZ(1)/2+abs(nZdZ(2)));
    if sZneeded > maxZ
      maxZ = sZneeded;
    end
    clear tomoNumber tiltName nZdZ
  end
  maxZ = maxZ + (samplingRate*2);
  fprintf('combining thickness and shift, found a maxZ of %d\n',maxZ);

    % xyzproj assumes centered in Z, so add extra height for z offsets to create
    % the true "in microsope" dimension

    reconstructionSize = [tiltHeader.nX,tiltHeader.nY,maxZ]
    originRec = ceil((reconstructionSize+1)./2)
    avgTomo = zeros(reconstructionSize,'single');
    avgSampling = zeros(reconstructionSize,'uint8');
    % These two are mutually exclusive for now, but not enforced.
    if (flgClassAvg)
      avgColor = zeros(reconstructionSize, 'int16');
    end

    if (flgColorMap)
      avgColor = zeros(reconstructionSize, 'int16');
    end

    if (buildTomo)
      coordOUT = fopen(sprintf('%smapBack%d/%s.coord',mbOUT{1:3}),'w');
      defOUT   = fopen(sprintf('%smapBack%d/%s.defAng',mbOUT{1:3}),'w');
    end
  % Track the number of fiducials in order to scale the K-factor to more or less
  % aggressivley downweight outliers in the alignment
  nFidsTotal = 0;
  for iTomo = 1:nTomograms



    TLT = tiltGeometry.(tomoList{iTomo});



    % Extract a "raw tilt" for alignment, make sure it is ordered properly
    % as the projection of the 3dModel with tilt will use this file and it
    % must match the zCoords in the defAng file.
    iRawTltName = sprintf('%smapBack%d/%s_align.rawtlt',mbOUT{1:3});
    iTiltFile = fopen(iRawTltName, 'w');
    rawTLT = sortrows(TLT(:,[1,4]),1);
    fprintf(iTiltFile,'%f\n',rawTLT(:,2)');
    fclose(iTiltFile); 
    
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
    
    % We also need the transform from the microscope frame in order to get
    % an accurate defocus value. Not sure if I should be binning?
    iXFName = sprintf('%smapBack%d/%s_align.XF',mbOUT{1:3});
    iXF = fopen(iXFName,'w');
    xfTLT = sortrows(TLT(:,[1,7:10,2,3],1));
    fprintf(iXF,'%f %f %f %f %f %f\n',xfTLT(:,2:7)');
    fclose(iXF);
    
    positionList = geometry.(tomoList{iTomo});
    positionList = positionList(positionList(:,26) ~= -9999,:);
    nFidsTotal = nFidsTotal + size(positionList,1);

    % Need to store tilt name/path explicity in meta deta
    tiltName    = tiltList{iTomo};
    

    tiltHeader = getHeader(MRCImage(tiltName,0));

    fullTiltSizeXandY = [tiltHeader.nX,tiltHeader.nY].*samplingRate;
 

    sTX = floor(tiltHeader.nX );
    sTY = floor(tiltHeader.nY );
    iTLT = floor(tiltHeader.nZ);

    originPrj = ceil(([sTX,sTY,0]+1)./2)

    tomoNumber = subTomoMeta.mapBackGeometry.tomoName.(tomoList{iTomo}).tomoNumber
    tiltName = subTomoMeta.mapBackGeometry.tomoName.(tomoList{iTomo}).tiltName
    reconCoords = subTomoMeta.mapBackGeometry.(tiltName).coords(tomoNumber,:)
            
    iGPU=1;
    tomoList{iTomo}
    if (buildTomo)
       [tomo,tomoReconCoords] = BH_multi_loadOrBuild(tomoList{iTomo}, ...
                                              reconCoords, mapBackIter, ...
                                              samplingRate, iGPU,reconScaling,1);
    
                                            
        

      tomoTrim = 4;
      tomo = tomo(1+tomoTrim:end-tomoTrim, ...
                  1+tomoTrim:end-tomoTrim, ...
                  1+tomoTrim:end-tomoTrim);
                
      tomo = tomo ./ (rmsScale*rms(tomo(:)));

      size(tomo)    
      size(avgTomo)
      sX = size(tomo,1);%tomoReconCoords(1,1)
      sY = size(tomo,2);%tomoReconCoords(1,2)
      sZ = size(tomo,3);%tomoReconCoords(1,3)
      originVol = ceil((size(tomo)+1)./2);%ceil((tomoReconCoords(1,1:3)+1)./2)

      reconShift = tomoReconCoords(2,1:3)
     
      % vector from first pixel in tilt series to lower left corner of tomogram from
      % which the subTomo origin is described XYZ in the geometry                                            
      lowerLeftVol = originRec + reconShift - originVol + 1
      
      % There is sometimes overlap, particularly with things like viral
      % capsids, so add to rather than just insert. Creates a little higher
      % density from the tomo in background but this is better than
      % replacing high res model density and an easy solution.
      
      try
        
        avgSampling(lowerLeftVol(1):lowerLeftVol(1)+sX -1, ...
                    lowerLeftVol(2):lowerLeftVol(2)+sY -1, ...
                    lowerLeftVol(3):lowerLeftVol(3)+sZ -1) = ...
         avgSampling(lowerLeftVol(1):lowerLeftVol(1)+sX -1, ...
                    lowerLeftVol(2):lowerLeftVol(2)+sY -1, ...
                    lowerLeftVol(3):lowerLeftVol(3)+sZ -1) +  ones(size(tomo),'uint8');
        overlapMask = (avgSampling(lowerLeftVol(1):lowerLeftVol(1)+sX -1, ...
                                   lowerLeftVol(2):lowerLeftVol(2)+sY -1, ...
                                   lowerLeftVol(3):lowerLeftVol(3)+sZ -1) <= 1);
        tomo = tomo .* overlapMask;
        avgSampling(avgSampling > 1) = 1;
        overlapMask = [];
        
        avgTomo(lowerLeftVol(1):lowerLeftVol(1)+sX -1, ...
                lowerLeftVol(2):lowerLeftVol(2)+sY -1, ...
                lowerLeftVol(3):lowerLeftVol(3)+sZ -1) = ...
        avgTomo(lowerLeftVol(1):lowerLeftVol(1)+sX -1, ...
                lowerLeftVol(2):lowerLeftVol(2)+sY -1, ...
                lowerLeftVol(3):lowerLeftVol(3)+sZ -1)      +        tomo;
              

      catch
        overShoot = size(avgTomo) - (lowerLeftVol + [sX,sY,sZ]);
        CLIPSIZE = max(abs(overShoot(:)));
        fprintf('WARNING, your tomo is too close to an edge. Overshoot x,y,z, %d %d %d.\nClipping tomo by %d\n',overshoot,CLIPSIZE);
        avgSampling(lowerLeftVol(1)+CLIPSIZE:lowerLeftVol(1)+sX -CLIPSIZE+1, ...
                lowerLeftVol(2)+CLIPSIZE:lowerLeftVol(2)+sY -CLIPSIZE+1, ...
                lowerLeftVol(3)+CLIPSIZE:lowerLeftVol(3)+sZ -CLIPSIZE+1) = ...
        avgSampling(lowerLeftVol(1)+CLIPSIZE:lowerLeftVol(1)+sX -CLIPSIZE+1, ...
                lowerLeftVol(2)+CLIPSIZE:lowerLeftVol(2)+sY -CLIPSIZE+1, ...
                lowerLeftVol(3)+CLIPSIZE:lowerLeftVol(3)+sZ -CLIPSIZE+1) + 1;
        overlapMask = (avgSampling(lowerLeftVol(1)+CLIPSIZE:lowerLeftVol(1)+sX -CLIPSIZE-1, ...
                                   lowerLeftVol(2)+CLIPSIZE:lowerLeftVol(2)+sY -CLIPSIZE-1, ...
                                   lowerLeftVol(3)+CLIPSIZE:lowerLeftVol(3)+sZ -CLIPSIZE-1) <= 1);

        tomo = tomo(CLIPSIZE+1:end-CLIPSIZE,CLIPSIZE+1:end-CLIPSIZE,CLIPSIZE+1:end-CLIPSIZE) .* overlapMask;
        avgSampling(avgSampling > 1) = 1;
       
        
        
    
        overlapMask = []; 
        % Should match, but allow a couple pixels of wiggle room
        avgTomo(lowerLeftVol(1)+CLIPSIZE:lowerLeftVol(1)+sX -CLIPSIZE-1, ...
                lowerLeftVol(2)+CLIPSIZE:lowerLeftVol(2)+sY -CLIPSIZE-1, ...
                lowerLeftVol(3)+CLIPSIZE:lowerLeftVol(3)+sZ -CLIPSIZE-1) = ...
        avgTomo(lowerLeftVol(1)+CLIPSIZE:lowerLeftVol(1)+sX -CLIPSIZE-1, ...
                lowerLeftVol(2)+CLIPSIZE:lowerLeftVol(2)+sY -CLIPSIZE-1, ...
                lowerLeftVol(3)+CLIPSIZE:lowerLeftVol(3)+sZ -CLIPSIZE-1)      +  tomo;%(CLIPSIZE+1:end-CLIPSIZE,CLIPSIZE+1:end-CLIPSIZE,CLIPSIZE+1:end-CLIPSIZE); 
           
      end      
              
      clear tomo   
    end

    

    nPrjs = size(TLT,1)
    nSubTomos = size(positionList,1);



    
    % need to update this.
    if (flgColorMap)
      colorMap = single(getVolume(MRCImage(COLOR_MAP)));
      % should be the same size as the average

      if any(size(refODD)-size(colorMap))
        error('Color map and average vol must be the same size.\n');
      end
        colorMap = colorMap(avgOrigin(1)-maxRad:avgOrigin(1)+maxRad,...
                    avgOrigin(2)-maxRad:avgOrigin(2)+maxRad,...
                    avgOrigin(3)-maxRad:avgOrigin(3)+maxRad);
    end              



    % Switch from maskRadius to particleRadius 20180129
    sizeAvgVol = size(refODD);
    if strcmpi(METHOD,'GPU')

      refODD = gpuArray(refODD);
      refEVE = gpuArray(refEVE);
      particleMask = BH_mask3d('sphere',sizeAvgVol,PARTICLE_RADIUS.*[1,1,1],[0,0,0]).* ...
                     BH_mask3d(refODD+refEVE,pixelSize,'','');
             
    else
      particleMask = BH_mask3d_cpu('sphere',sizeAvgVol,PARTICLE_RADIUS.*[1,1,1],[0,0,0]);
    end

%     [ particleMask ] = BH_multi_randomizeTaper(particleMask);



    % Apply high (and low) pass in 3d - when applying the particles CTF, soften
    % the high-pass portion in 2d
    bandPass = BH_bandpass3d(sizeAvgVol,0,0,lowPassCutoff,METHOD,pixelSize);

    %Should read in FSC value, but just use generic for now.
    padVal = [0,0,0;0,0,0]


    if (iTomo == 1)
      fidIDX = 0;
    end

    %%%%%%%%%%
    if (buildTomo)

      modelRot = BH_defineMatrix([0,90,0],'Bah','forwardVector');

    for iSubTomo = 1:nSubTomos


      rSubTomo = reshape(positionList(iSubTomo,17:25),3,3);



      prjVector = (positionList(iSubTomo,11:13)./samplingRate) - ...
                                         tomoTrim - originVol + reconShift;
      
      
     
    


% % %       prjVector = prjVector + [0.5,0.0,-0.5];
      prjVector = prjVector + [0.0,0.0,1.0];
      recVector = (originPrj + [0,0,ceil((reconstructionSize(3)+1)/2)] + prjVector); % subTomo origin relative to reconLowerLeft

      %Resample a copy of the average to match the position in the tomogram
      % The third entry is a dummy, normally used to make sure at least the
      % particle was being extracted even if the surrounding density (where
      % some delocalized values may be located) are not.
      [ indVAL, padVAL, shiftVAL ] = ...
                          BH_isWindowValid(reconstructionSize, sizeAvgVol, sizeAvgVol./5, recVector);
       

      if ~ischar(indVAL)                  

       
        % Reproject using tilt, so just save the 3d coords. 
       fprintf(coordOUT,'%0.4f %0.4f %0.4f %d\n', modelRot*prjVector' + [originRec(1),originRec(3),originRec(2)]', fidIDX);
       if positionList(iSubTomo,7) == 1
        iAvgResamp = BH_resample3d(refODD,rSubTomo',shiftVAL,'Bah',METHOD,'forward');  
       elseif positionList(iSubTomo,7) ==2
         iAvgResamp = BH_resample3d(refEVE,rSubTomo',shiftVAL,'Bah',METHOD,'forward');
       else
         error('positionList iSubtomo %d col 7 is %d',iSubTomo,positionList(iSubTomo,7));
       end
       iMaskResamp = BH_resample3d(particleMask,rSubTomo',shiftVAL,'Bah',METHOD,'forward');
         

        iAvgResamp = gather(iMaskResamp.*iAvgResamp);

        if iSubTomo == 1
          SAVE_IMG(MRCImage(gather(iAvgResamp)),'testResampleMasked.mrc');
        end
        if (flgColorMap)
          iColorMap = gather(int16(particleMask.* BH_resample3d(colorMap, ...
                                   rSubTomo',shiftVAL,'Bah',METHOD,'forward')));  


          avgColor(indVAL(1,1):indVAL(2,1), ...
                   indVAL(1,2):indVAL(2,2), ...
                   indVAL(1,3):indVAL(2,3)) = avgColor(indVAL(1,1):indVAL(2,1), ...
                                                       indVAL(1,2):indVAL(2,2), ...
                                                       indVAL(1,3):indVAL(2,3)) + ...
                iColorMap(1+padVAL(1,1):end-padVAL(2,1),...
                                                      1+padVAL(1,2):end-padVAL(2,2),...
                                                      1+padVAL(1,3):end-padVAL(2,3)); 
        end
        avgTomo(indVAL(1,1):indVAL(2,1), ...
                indVAL(1,2):indVAL(2,2), ...
                indVAL(1,3):indVAL(2,3)) =  ...
                                    avgTomo(indVAL(1,1):indVAL(2,1), ...
                                            indVAL(1,2):indVAL(2,2), ...
                                            indVAL(1,3):indVAL(2,3)) .* ...
                        gather((1 -  iMaskResamp(1+padVAL(1,1):end-padVAL(2,1),... % zeros out region being replaced
                                          1+padVAL(1,2):end-padVAL(2,2),...
                                          1+padVAL(1,3):end-padVAL(2,3)))) + ...                                           
                                iAvgResamp(1+padVAL(1,1):end-padVAL(2,1),...
                                           1+padVAL(1,2):end-padVAL(2,2),...
                                           1+padVAL(1,3):end-padVAL(2,3)); 



        for iPrj = 1:nPrjs

          % imod is indexing from zero
          zCoord = TLT(iPrj,1) - 1; 
          % These shifts are a record of transformation from the raw data, but here
          % we are comparing with [CTF] corrected data, from which the
          % reconstructino was made directly


          rTilt = BH_defineMatrix([1.*TLT(iPrj,6),1.*TLT(iPrj,4),-1*TLT(iPrj,6)],'Bah','forwardVector');

           R = rTilt;
  

          % The prjCoords are relative to the origin of the tilt series, but the
          % model needs to be relative to lower left, but the shift vector to do
          % this must be transformed for each tilt.

          shiftCoords = originPrj';

         prjCoords = R*prjVector';

          fprintf(defOUT,'%d %d %6.6e\n', fidIDX, zCoord, samplingRate.*prjCoords(3).*pixelSize.*10^-10+TLT(iPrj,15));

        end % loop over tilt projections      


       fidIDX = fidIDX + 1;
      else
        fprintf('ignoring subTomo %d for out of bounds conditions.\n', iSubTomo);
      end
    end % loop over subtomos
    end %%%% temp condition to skip building full tomo
    

  end


  if (buildTomo)
   fclose(coordOUT);
   

    p2m = sprintf(['point2model -zero -circle 3 -color 0,0,255 -values -1 ',...
                   '%smapBack%d/%s.coord %smapBack%d/%s.3dfid'], ...
                   mbOUT{1:3},mbOUT{1:3})
    system(p2m);

    avgTomo = MRCImage(avgTomo);
    SAVE_IMG(avgTomo,sprintf('%smapBack%d/%s.tmpTomo', mbOUT{1:3}),4.0);
    clear avgTomo
    % If not planning on visualization, save only a binned copy of the synthetic
    % tomo.
    if ~(flgColorMap) && ~(conserveDiskSpace)
      tmpTomoBin = floor(1/samplingRate*6);
      system(sprintf(['binvol -bin %d %smapBack%d/%s.tmpTomo ',...
                       '%smapBack%d/%s.bin%dTomo'], ...
                       tmpTomoBin,mbOUT{1:3},mbOUT{1:3},tmpTomoBin));
  % % % % %     system(sprintf('rm mapBack/%s.tmpTomo', tiltBaseName));
    end

    % -90 is assumed for trim vol, so if rotate vol is used add 90
  
  
% % %       extraRot = 90;
% % %       reconRotation(iTomo,3)+90.0+extraRot
      rotSize = [tiltHeader.nX,maxZ,tiltHeader.nY]

      rotCMD = sprintf(['rotatevol -angles 0,0,90 -size %d,%d,%d ',...
               '%smapBack%d/%s.tmpTomo %smapBack%d/%s.tmpRot'], ...
               rotSize, mbOUT{1:3},mbOUT{1:3});
      system(rotCMD);
      system(sprintf('mv %smapBack%d/%s.tmpRot %smapBack%d/%s.tmpTomo', ...
                     mbOUT{1:3},mbOUT{1:3}));

 

    % fix me as above
  % % % % %   if (flgColorMap)
  % % % % %     SAVE_IMG(MRCImage(avgColor),sprintf('mapBack/%s_colorMap.mrc',tiltBaseName));
  % % % % %     % -90 is assumed for trim vol, so if rotate vol is used add 90
  % % % % %     if (rotateVol)
  % % % % %       system(sprintf('rotatevol -angles 0,0,%d mapBack/%s_colorMap.mrc mapBack/%s_colorMap.rot',reconRotation(iTomo,3)+90.0,tiltBaseName,tiltBaseName));
  % % % % %       system(sprintf('mv mapBack/%s_colorMap.rot mapBack/%s_colorMap.mrc',tiltBaseName,tiltBaseName));
  % % % % %       system(sprintf('rm mapBack/%s_colorMap.rot',tiltBaseName));
  % % % % %     end 
  % % % % %   end


    clear avgTomo  wgt
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
        
 
        keepItRunning = 1;
        outputStackName = sprintf('%smapBack%d/%s_mapBack.st',mbOUT{1:3});
        
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
              rePrjFileName = sprintf('%smapBack%d/%s_rePrj.sh',mbOUT{1:3});
              reModFileName = sprintf('%smapBack%d/%s_reMod.sh',mbOUT{1:3});
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
                                  'RecFileToReproject %smapBack%d/%s.tmpTomo\n',...
                                  'TOTALSLICES %d,%d\n',...
                                  'ZMinAndMaxReproj %d,%d\n',...
                                  '%s\n', ...
                                  '%s\n', ...                          
                                  '%s\n',...
                                  'EOF'],tiltList{1} ,outputStackName, maxZ, ...
                                         mbOUT{1:3},...
                                         taStr, mbOUT{1:3},...
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
              rePrjFileName = sprintf('%smapBack%d/%s_rePrj.sh',mbOUT{1:3});
              reModFileName = sprintf('%smapBack%d/%s_reMod.sh',mbOUT{1:3});
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
                                  'RecFileToReproject %smapBack%d/%s.tmpTomo\n',...
                                  'TOTALSLICES %d,%d\n',...
                                  'ZMinAndMaxReproj %d,%d\n',...
                                  '%s\n', ...
                                  '%s\n', ...                          
                                  '%s\n',...
                                  'EOF'],tiltList{1} ,outputStackName, maxZ, ...
                                         mbOUT{1:3},...
                                         taStr, mbOUT{1:3},...
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

    % Remove the full size tomo
    system(sprintf('rm %smapBack%d/%s.tmpTomo',mbOUT{1:3}));  

    fidList = load(sprintf('%smapBack%d/%s.coordPrj',mbOUT{1:3}));
    defList = load(sprintf('%smapBack%d/%s.defAngTilt',mbOUT{1:3}));

    % results.

    % I am pretty sure that the existing first column is the fiducial id
    % incremented from zero. This step could be removed and indexed over i
    % + 1 in later steps. (or saved as indexed from one, but I don't
    % remember what that would effect.)
    % Now from zero: fidIDX, X, Y, iPrj(0)
%     fidList = [1:size(fidList,1);fidList']';
%     defList = [1:size(defList,1);defList']';
    defList = defList(:,[1,7,3]);
    defList(:,[1,3]) = defList(:,[1,3]) - 1;
    defList(:,2) = defList(:,2).*(-1*10^-9);

     % Give every instance of each fiducial a unique identifier.
     fidList = [1:size(fidList,1);fidList']';
%     defList = [1:size(defList,1);defList']';


    % for center of mass
    COM = 3;
    [bx,by] = ndgrid(-COM:COM,-COM:COM)
    % add optional half radius for edge case and make the padding more
    % logical, twice the particle radius, and then CTF size using mulit_iter
    % with an optimization step
    particlePad = 1.5;
    tileRadius = floor(particlePad.*PARTICLE_RADIUS);
    tileSize = (2.*tileRadius + 1).*[1,1];
    padTile = floor(PARTICLE_RADIUS./1 .*[1,1]);
%     padCTF = 3.*padTile; 
    CTFSIZE = BH_multi_iterator([1 + 2* tileRadius + 6.*padTile,1], 'fourier')
    CTFSIZE = CTFSIZE(1:2)
    padCTF = BH_multi_padVal(tileSize,CTFSIZE);
    
   
    if (eraseMask)
      peakMask = BH_mask3d(eraseMaskType,(tileSize(1)+2.*padTile(1)).*[1,1],eraseMaskRadius,[0,0],'2d');
      peakMaskDefSearch = BH_mask3d(eraseMaskType,CTFSIZE,eraseMaskRadius,[0,0],'2d');
      peakMask(peakMask < 0.99) = 0;
      peakMaskDef(peakMaskDefSearch < 0.99) = 0;
    else
      peakMask = BH_mask3d('sphere',(tileSize(1)+2.*padTile(1)).*[1,1],peakSearchRad,[0,0],'2d');
      peakMaskDefSearch = BH_mask3d('sphere',CTFSIZE,peakSearchRad,[0,0],'2d');
    end

    fftMask = BH_fftShift(0,(tileSize(1)+2.*padTile(1)).*[1,1],1);
    fftMaskDefSearch = BH_fftShift(0,CTFSIZE,1);
    [dU, dV] = BH_multi_gridCoordinates(CTFSIZE,'Cartesian',METHOD, ...
                                                    {'none'},1,1,0);
    dU = dU .* (2i*pi);
    dV = dV .* (2i*pi);  


    
%     [ peakMask ] = BH_multi_randomizeTaper(peakMask);

    tileOrigin = ceil((size(peakMask)+1)./2);
    defSearchOrigin = ceil((size(peakMaskDefSearch)+1)./2);
%     bandPassFilter = BH_bandpass3d([tileSize+2.*padCTF,1], 0,0,lowPassCutoff,METHOD,pixelSize);
    bandPassFilter = BH_bandpass3d([CTFSIZE,1], 0,0,lowPassCutoff,METHOD,pixelSize);
    
    bandPassPrj = BH_bandpass3d([sTX,sTY,1],0,0,lowPassCutoff,'cpu',pixelSize);

    diagnosticCell = cell(nPrjs,1);
    evalMaskCell = cell(nPrjs,1);

    for iPrj = 1:nPrjs
      evalMaskCell{iPrj} = zeros(gather([sTX,sTY,1]),'uint8');
    end

    cccPrecisionTaper = 'singleTaper';
    cccPrecision = 'single';

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

    nToCheck = floor(ctfRange./ctfInc);
    defShiftVect = ctfInc.*[-nToCheck:nToCheck]';
    nDefTotal = length(defShiftVect);
    defocusCCC = cell(nPrjs,1);
    expectedDefocusPerFiducial=cell(nPrjs,1);

    nFidsTotal = numel(unique(fidList(:,1)));
    for iPrj = 1:nPrjs
      % I must specify the number of fiducials somehwere else, replace the
      % unique when there is time.
      defocusShifts{iPrj} = 0;
      defocusCCC{iPrj} = zeros(nDefTotal, nFidsTotal,'single');
      expectedDefocusPerFiducial{iPrj} = zeros(nDefTotal,nFidsTotal,'single');
    end
    
    

      if samplingRate > 1
        if ( flg2dCTF )
          tiltSeries = sprintf('%scache/%s_ali%d_ctf_bin%d.fixed',CWD,tiltName,mapBackIter+1,samplingRate)
        else
          tiltSeries = sprintf('%scache/%s_ali%d_bin%d.fixed',CWD,tiltName,mapBackIter+1,samplingRate)
        end
      else
        if ( flg2dCTF )
          
          tiltSeries = sprintf('%sctfStacks/%s_ali%d_ctf.fixed',CWD,tiltName,mapBackIter+1)
        else
          tiltSeries = sprintf('%saliStacks/%s_ali%d.fixed',CWD,tiltName,mapBackIter+1)
        end
      end

      
parfor iPrj = 1:nPrjs     
% for iPrj = 20;%1:nPrjs
	    % For some reason if these mrc objects are created before the parfor
	    % loop begins, they fail to load. It is fine as a regular for loop
	    % though - annoying, but very little overhead. It would be nice
	    % to know what is going on here.
      iMrcObj = MRCImage(tiltSeries);
      iMrcObjRef = MRCImage(sprintf('%smapBack%d/%s_mapBack.st',mbOUT{1:3}));
      
      % Matching the "natural" or sequential order 
      iTLT = find(TLT(:,1) == iPrj);
     
      tic
      while toc < 300
        try
          dataPrj = single(getVolume(iMrcObj,[1,sTX],[1,sTY],iPrj));
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
          refPrj  = single(getVolume(iMrcObjRef,[1,sTX],[1,sTY],iPrj));
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
      mRms = mean(dataRMS(:))
      sRms = rms(dataRMS(:)-mRms)
      
      
      
      
      figure, imshow3D(gather(dataPrj));
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
                                  
      meanDef = TLT(iTLT,15);          
      defAst = TLT(iTLT,12);
      angAst = TLT(iTLT,13);
      defVect = [meanDef - defAst, meanDef + defAst, angAst];
      
      [Hqz, HqzUnMod] = BH_ctfCalc(TLT(iTLT,16).*samplingRate,TLT(iTLT,17), ...
                                   TLT(iTLT,18),defVect,size(refPrj), ...
                                   1.*TLT(iTLT,19),-0.15);

      Hqz = gather(Hqz);
      HqzUnMod = gather(HqzUnMod);
      
        if ( flg2dCTF )
          % If not ctf corrected projections, just use the ctf directly, otherwise...
          Hqz = abs(Hqz.*HqzUnMod);
        end

      % Add abs(HqzUnMod) to make the data prj amplitudes match those of th
      % reference more accurately.
      cccPrj = fftshift(real(ifftn(bandPassPrj.*fftn(dataPrj).* abs(HqzUnMod).*...
                              conj(fftn(refPrj).*Hqz))));

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
        error('sdf')
      end
      cccPRJ = cccPRJ - min(cccPRJ(:));
      comPRJX = sum(sum(bx.*cccPRJ))./sum(cccPRJ(:));
      comPRJY = sum(sum(by.*cccPRJ))./sum(cccPRJ(:));




      estPeak = [mRx, mRy] - originPrj(1:2) + [comPRJX, comPRJY];
      glbList = fopen(sprintf('%smapBack%d/%s_%03d.global',mbOUT{1:3},iPrj),'w');
    % Add unique indicies to prevent ambiquity when comparing with paral
      fprintf(glbList,'%f  degree tilt at %f %f\n', TLT(iTLT,4),estPeak);
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
      % The defocus list is ordered like the TLT info which is not
      % sequential, so it needs it's own logical.
      wrkDefIDX = ( defList(:,3) == iPrj - 1 );
      wrkDef = defList(wrkDefIDX,:);
      
      coordOUT = fopen(sprintf('%smapBack%d/%s_%03d.coordFIT',mbOUT{1:3},iPrj),'w');
% % defInter = fopen(sprintf('%smapBack%d/%s%s_%03d.defInter',mbOUT{1:3},outCTF,iPrj),'w');
      [radialForCTF,phi,~,~,~,~] = BH_multi_gridCoordinates([CTFSIZE,1],'Cylindrical',METHOD,{'none'},1,0,0);


      radialForCTF = {radialForCTF./(pixelSize.*10^-10),0,phi};
      phi = [];

      calcPeakShifts = 1;
 
      
      if (calcPeakShifts && calcCTF)
        % should probably calculate the low pass closer to the first zero.
        % Here again just assuming < 8um. Also need a check in case the
        % lowPassCutoff is too low to produce meaningful ctf comparison.
        ctfMask = BH_bandpass3d([CTFSIZE,1],10^-1,40,max(5,2*pixelSize),METHOD,pixelSize);
        % Testing an additional translational step which requires full fft
        % band mask for just the hermitian symmetry
        ctfMask(:,ceil((CTFSIZE(1)+1)/2):end) = 0;
        ctfMask = (ctfMask > 10^-2 );
        
        % find range of defocus for this projection.

        defToCheck = min(wrkDef(:,2))-ctfRange-ctfInc:ctfInc:max(wrkDef(:,2))+ctfRange+ctfInc;
        nCTFs = length(defToCheck)
        if strcmpi(METHOD,'GPU')
          ctfStack = zeros(sum(ctfMask(:)),nCTFs,'single','gpuArray');
        else
          ctfStack = zeros(sum(ctfMask(:)),nCTFs,'single');
        end % % %         else
  
% % % % %         if strcmpi(METHOD,'GPU')
% % % % %           ctfStack = zeros([CTFSIZE,nCTFs],'single','gpuArray');
% % % % %         else
% % % % %           ctfStack = zeros([CTFSIZE,nCTFs],'single');
% % % % %         end
% % % % %         
        for iCTF = 1:nCTFs
          defAst = TLT(iTLT,12);
          angAst = TLT(iTLT,13);
          defVect = [defToCheck(iCTF) - defAst, defToCheck(iCTF) + defAst, angAst];
          [Hqz, ~] = BH_ctfCalc(radialForCTF,TLT(iTLT,17),TLT(iTLT,18), ...
                                    defVect, ...
                                    CTFSIZE, ...
                                    -1.*TLT(iTLT,19), ...
                                    -1.0);
          ctfStack(:,iCTF) = Hqz(ctfMask);
% % % % %           ctfStack(:,:,iCTF) = Hqz;

          
        end
        
        % Set peakShift calculation to false, this is set to be true after
        % the first loop.
        calcPeakShifts = 0;
      end
      
      for iFidLoop = 1:1+calcCTF
        
      for iFid = 1:size(wrkFid,1)

        ox = floor(wrkFid(iFid,3)) - tileRadius;
        oy = floor(wrkFid(iFid,4)) - tileRadius;
        oxEval = [floor(wrkFid(iFid,3) - PARTICLE_RADIUS),floor(wrkFid(iFid,3) + PARTICLE_RADIUS)];
        oyEval = [floor(wrkFid(iFid,4) - PARTICLE_RADIUS),floor(wrkFid(iFid,4) + PARTICLE_RADIUS)];
        % it would be good to try a smaller tile.
	% First check that the data are found in this given projection
       
        try 
          % If any zeros values within the particle radius, do not evaluate
          iSkipEval = any(any(evalMask(oxEval(1):oxEval(2),oyEval(1):oyEval(2)) == 0));
        catch
          % If the particle was outof bounds, do not evaluate
          iSkipEval = 1;
        end
        
        if ( iSkipEval )
          fprintf('\nThe current fiducial %d is not sampled in this projection %d\n',iFid,iPrj);
          if iFidLoop == 1+calcCTF
            % Only print out to file if doing the alignment
            fprintf(coordOUT,'%d %d %0.4f %0.4f %d\n', wrkFid(iFid,1:2), [2,-2], -9999);
          end
          continue
        end
  
        if  (ox < 1 || oy < 1 || ox +2*tileRadius > sTX || oy +2*tileRadius > sTY )
          fprintf('\nThe current fiducial is too close to the edge, ox %d oy %d 2x Rad %d\n',ox,oy,2*tileRadius)
          if iFidLoop == 1+calcCTF
            % Only print out to file if doing the alignment
            fprintf(coordOUT,'%d %d %0.4f %0.4f %d\n', wrkFid(iFid,1:2), [2,-2], -9999);
          end
          continue
        end

        dataTile = dataPrj(ox:ox+2.*tileRadius,oy:oy+2.*tileRadius);
        dataTile = dataTile - mean(dataTile(:));
       
        
        refTile = refPrj(ox:ox+2.*tileRadius,oy:oy+2.*tileRadius);         
        refTile = refTile - mean(refTile(:));
       
        
        dataTile = BH_padZeros3d(dataTile,padCTF(1,:),padCTF(2,:), ...
                                                 METHOD,cccPrecisionTaper); 

        refTile = BH_padZeros3d(refTile,padCTF(1,:),padCTF(2,:), ...
                                                 METHOD,cccPrecisionTaper);  
                            
        if (calcPeakShifts)               
          % The original idea was to adjust the mean defocus applied as
          % will be done in the ctf correction. The shifts in XY depend
          % strongly on the defocus, and there is no reason to anticipate
          % the sample to move uniformly (the entire point of this
          % program.) Try shifting by max scoring defocus per fiducial, and
          % maybe eventually correcting that way too.
% % % % %           meanDef = wrkDef(iFid,4)+defocusShifts{iPrj}; 
% % %           [~,iMaxDef] = max(defocusCCC{iPrj}(:,iFid));
% % %            
% % %           meanDef = wrkDef(iFid,4)+defShiftVect(iMaxDef); 
          meanDef = wrkDef(iFid,2) + expectedDefocusPerFiducial{iPrj}(1,iFid);
          defAst = TLT(iTLT,12);
          angAst = TLT(iTLT,13);
          defVect = [meanDef - defAst, meanDef + defAst, angAst];

          % If true, sets the ampCont to -ampCont which tells ctfCalc to use
          % an envelope function on the CTF.
          flgDampen = 0;
          [Hqz, HqzUnMod] = BH_ctfCalc(radialForCTF,TLT(iTLT,17),TLT(iTLT,18), ...
                                      defVect, ...
                                      CTFSIZE, ...
                                      TLT(iTLT,19).*(1-2*flgDampen), ...
                                      -0.15);

          if ( flg2dCTF )
            % If not ctf corrected projections, just use the ctf directly, otherwise...
            Hqz = abs(Hqz.*HqzUnMod);
          end
        

          refTile = real(ifftn(fftn(refTile).*Hqz.*bandPassFilter));
          refTile =  BH_padZeros3d(refTile,-1.*padCTF(1,:), ...
                                           -1.*padCTF(2,:),METHOD, ...
                                                              cccPrecision);            


          dataTile = real(ifftn(fftn(dataTile).*bandPassFilter));
          dataTile =  BH_padZeros3d(dataTile,-1.*padCTF(1,:), ...
                                           -1.*padCTF(2,:),METHOD, ...
                                                              cccPrecision);  

	  try
            dataTile = dataTile - mean(dataTile(:));
            rmsData = rms(dataTile(:))
            dataTile = dataTile ./ rmsData;

            refTile = refTile - mean(refTile(:));
            rmsRef = rms(refTile(:));
            refTile = refTile ./ rmsRef;


            % There shouldn't be out of bounds, but since the data is masked, there
            % might be occassional zero values
            if ( rmsData && rmsRef )

              phaseOnly = 0;
              if (phaseOnly)
                dataTile = exp(1i.*angle(fftn(BH_padZeros3d(dataTile, padTile, padTile, METHOD, cccPrecision))));
                refTile = exp(1i.*angle(conj((fftn(BH_padZeros3d(refTile, padTile, padTile, METHOD, cccPrecision))))));
              else
                dataTile = fftn(BH_padZeros3d(dataTile, padTile, padTile, METHOD, cccPrecision));
                refTile =  conj(fftn(BH_padZeros3d(refTile, padTile, padTile, METHOD, cccPrecision)));
              end

              cccMap = real(ifftn(dataTile.*refTile));
             
              cccMap = cccMap(fftMask) .* peakMask;

              [~,maxMap] = max(cccMap(:));

              [mMx, mMy] = ind2sub(size(cccMap), maxMap);


              cccMap = cccMap(mMx-COM:mMx+COM, mMy-COM:mMy+COM);

              cccMap = cccMap - min(cccMap(:));

              comMapX = sum(sum(bx.*cccMap))./sum(cccMap(:));
              comMapY = sum(sum(by.*cccMap))./sum(cccMap(:));

                 % peak in Map is where query is relative to ref, dXY then is the shift
                 % needed to move the predicted position to the measured.
     %           dXY = -1.*[mRx,mRy] + -1.*[comRefX,comRefY]+[mMx,mMy]+[comMapX,comMapY];
                % Data moved from a position of estPeak, so add this to dXY
                dXY = [mMx,mMy]+[comMapX,comMapY] - tileOrigin(1:2)+ estPeak;
     %           dXY = dXY + estPeak;

                fprintf(coordOUT,'%d %d %0.4f %0.4f %d\n', wrkFid(iFid,1:2), dXY, wrkFid(iFid,5));
            else
              fprintf('rms is out for this tile %d %d\n',rmsData,rmsRef);
              fprintf(coordOUT,'%d %d %0.4f %0.4f %d\n', wrkFid(iFid,1:2), [0,0], -9999);
            end % if condition on name

          catch
            fprintf('falling out of the try/catch over peak search\n');
            fprintf(coordOUT,'%d %d %0.4f %0.4f %d\n', wrkFid(iFid,1:2), [-1,-1], -9999);
          end
        
       
        else
 
          [~,meanIDX] = min(abs(defToCheck-wrkDef(iFid,2)));
          
          iData  = fftn(dataTile);
          iData(1) = 0;
          iData = iData(ctfMask);
          iDataNorm = sum(abs(iData(:)).^2);
% % % % %           dataTile = dataTile(ctfMask);
% % % % %           dataNorm = sum(abs(dataTile(:)).^2);
               
          iRef = conj(fftn(refTile));
          iRef(1) = 0;
          iRef = iRef(ctfMask);
% % % % %           refTile = refTile(ctfMask);
          
% % % % %           for iCTF = -nToCheck:nToCheck
% % % % %             iRef = conj(refTile.*ctfStack(:,meanIDX+iCTF));
% % % % %             iCCC = real(sum(iRef.*dataTile))./sqrt(sum(abs(iRef(:)).^2).* dataNorm);
% % % % %             defocusCCC{iPrj}(iCTF+nToCheck+1,iFid) = gather(iCCC);
% % % % %           end
          for iCTF = -nToCheck:nToCheck
            iRefCTF = iRef.*ctfStack(:,meanIDX+iCTF);
            iRefNorm = sum(abs(iRefCTF(:)).^2);
% % % % %             cccMap = real(ifftn(iRef.*iData));
% % % % %             cccMap = cccMap(fftMaskDefSearch).*peakMaskDefSearch;
% % % % %             
% % % % %             try
% % % % %               [~,maxMap] = max(cccMap(:));
% % % % %               [mMx, mMy] = ind2sub(size(cccMap), maxMap);
% % % % %               cccMap = cccMap(mMx-COM:mMx+COM, mMy-COM:mMy+COM);
% % % % %               cccMap = cccMap - min(cccMap(:));
% % % % %               comMapX = sum(sum(bx.*cccMap))./sum(cccMap(:));
% % % % %               comMapY = sum(sum(by.*cccMap))./sum(cccMap(:));
% % % % %               dXY = [mMx,mMy]+[comMapX,comMapY] - defSearchOrigin(1:2);
% % % % %             catch
% % % % %               fprintf('peak search in def search went awry');
% % % % %               dXY = [0,0];
% % % % %             end
% % % % %             % Shift prior to calculating the CCC
            
% % % %             iCCCun = real(sum(sum(iRef.*iData)))./sqrt(iUnShift.*iRefNorm);
% % % % %             iData = exp(dU.*dXY(1)+dV.*dXY(2)).*iData;  
% % % % %             iDataNorm = sum(abs(iData(:)).^2);
            iCCC = real(sum(iRefCTF.*iData))./sqrt(iDataNorm.*iRefNorm);
             
            defocusCCC{iPrj}(iCTF+nToCheck+1,iFid) = gather(iCCC);
%             fprintf(defInter,'%d %3.3e %3.3e %3.3f %3.3f %4.4f %4.4f\n',iFid, ...
%               defToCheck(meanIDX),defToCheck(meanIDX+iCTF),dXY,iCCC, iCCCun);
          end
        end
      end % end of loop over fiducials

      % First identify expected value for each fiducial. The distriution of
      % maximums has much greater variance than the distribution of
      % expected values, and shifting each fiducial by the prior degrades
      % the fitting.
      tmpCCC = defocusCCC{iPrj};
%       tmpCCC( tmpCCC < 0 ) = 0;
      % Truncating negative values of CCC to zero actually makes those
      % scores MORE likely
      tmpCCC = tmpCCC - min(tmpCCC(:));
      tmpCCC = (tmpCCC ./ max(tmpCCC,[],1)) .^ exp(probabilityPeakiness);
      tmpCCC = tmpCCC ./ sum(tmpCCC,1);
      expectedDefocusPerFiducial{iPrj} = sum(tmpCCC.* ...
                                       repmat(reshape(defShiftVect,nDefTotal,1), ...
                                       1,size(tmpCCC,2)),1);
      tmpCCC = [];
  
      fineDefocusVect = defShiftVect(1):ctfInc/10:defShiftVect(end);
      maxProbHist = hist(expectedDefocusPerFiducial{iPrj},fineDefocusVect);
      maxProbHist = maxProbHist ./ sum(maxProbHist);
      
      expectedDefocus = sum(maxProbHist.*fineDefocusVect);
      % Max for each fiducial, return index 
% % %       [~,mHist] = max(defocusCCC{iPrj});
% % %       maxProbHist = hist(mHist+1,1:length(shiftVect));  
% % %       maxProbHist = hist(mHist,1:nDefTotal);
% % % expectedShiftIDX = round(sum(maxProbHist .* [1:length(defShiftVect)]));      
      
      
% % %       if expectedShiftIDX == 0
% % %         expectedShiftIDX = 1;
% % %       elseif expectedShiftIDX > length(defShiftVect)
% % %         expectedShiftIDX = length(defShiftVect);
% % %       end
% % %       for iDefShift = 1:length(shiftVect)
% % %         iProb = defocusCCC{iPrj}(:,iDefShift);
% % %         iProb(iProb < 0) = 0;
% % %         iProb(isnan(iProb)) =  0;
% % %         if any(iProb)
% % %           iProb = iProb .* maxProbHist';
% % %           if (sum(iProb))
% % %             iProb = iProb ./ sum(iProb);
% % %             nReasonable = nReasonable + 1;
% % %             expectedDefocusShift = expectedDefocusShift + sum(shiftVect .* iProb);
% % %           end
% % %         end
% % %         end
% % %       if (nReasonable)
% % %         expectedDefocusShift = (expectedDefocusShift./nReasonable);
% % %       else
% % %         expectedDefocusShift = 0;
% % %       end
% % %       defocusShifts{iPrj} = expectedDefocusShift;
% % %       defocusShifts{iPrj} = defShiftVect(expectedShiftIDX);
      defocusShifts{iPrj} = expectedDefocus;

      calcPeakShifts = 1;
      fprintf('prj %d delDef %3.3e\n',expectedDefocus);
      end % if ctfs are calculated loop again over fiducials
    evalMaskCell{iPrj} = uint8(evalMask); evalMask = [];
    fclose(coordOUT);
 end % end of the parfor loop
    
    if ( calcCTF )
      defShifts = fopen(sprintf('%smapBack%d/%s%s.defShifts',mbOUT{1:3},outCTF),'w');
      defCCC = sprintf('%smapBack%d/%s%s_defCCC.mat',mbOUT{1:3},outCTF);
      save(defCCC,'defocusCCC','expectedDefocusPerFiducial');
      for iPrj = 1:nPrjs
        fprintf(defShifts,'%6.6e\n',defocusShifts{iPrj});
      end
    end
    evalMaskStack = zeros(sTX,sTY,nPrjs);
    diagnosticStack = zeros([(globalPeak.*2+1).*[1,1],nPrjs],'single');
    for iPrj = 1:nPrjs
      diagnosticStack(:,:,iPrj) = gather(diagnosticCell{iPrj});
      evalMaskStack(:,:,iPrj) = int16(gather(evalMaskCell{iPrj}));
    end
    clear diagnosticCell evalMaskCell
    if ~(conserveDiskSpace)
      SAVE_IMG(MRCImage(diagnosticStack),sprintf('%smapBack%d/%s_diagnostic.mrc',mbOUT{1:3}));
      SAVE_IMG(MRCImage(evalMaskStack),sprintf('%smapBack%d/%s_evalMask.mrc',mbOUT{1:3}));
    end
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
    fCombine = [fidShifts(:,1),fidList(:,2:3)+fidShifts(:,2:3),fidShifts(:,4)];
    fprintf('\n\n%d/%d pts ignored\n\n',sum(fCombine(:,4)==-9999),size(fCombine,1));
    
    
    fFull = fCombine;
    fDefFull = [fCombine,zeros(size(fCombine,1),1)];
% % % % %     fDefFull(:,2:3) = fDefFull(:,2:3).*pixelSize;
    fDefFull(:,2:3) = fDefFull(:,2:3).*samplingRate;
    for iPrj = 1:nPrjs
      % Create a file that has the X,Y,defocus positions for all fiducials
      % in each tilt to use in ctf correction.
      wrkDefIDX = ( defList(:,3) == iPrj - 1 );
      wrkFidIDX = ( fidList(:,4) == iPrj - 1 );
      fDefFull(wrkFidIDX,5) = defList(wrkDefIDX,2) + defocusShifts{iPrj};
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
    fprintf(aliCom,['#!/bin/bash\n\n',...
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
                    'SurfacesToAnalyze	2\n',...
                    'RotOption	1\n',... % def solve all rotations
                    'RotDefaultGrouping	3\n',... % if rot option --> 5 use def group size
                    'TiltOption	%d\n',... % Tilts are harder use automapping
                    'TiltDefaultGrouping	%d\n',...                 
                    'MagOption	1\n',... % def solve all mags
                    'MagDefaultGrouping	3\n',...
                    'XStretchOption	0\n',...
                    'SkewOption	0\n',...          
                    'BeamTiltOption	0\n',...  
                    'XTiltOption	0\n',...
                    'ResidualReportCriterion	0.001\n',...
                    'ShiftZFromOriginal\n',...
                    'AxisZShift 0.0\n',...
                    'RobustFitting\n',...
                    'KFactorScaling %3.3f\n',...
                    'LocalAlignments\n',...
                    'LocalRotOption 1\n',...
                    'LocalRotDefaultGrouping 3\n',...
                    'LocalTiltOption %d\n',...
                    'LocalTiltDefaultGrouping %d\n',...
                    'LocalMagOption %d\n',...
                    'LocalMagDefaultGrouping 5\n',...
                    'OutputLocalFile %smapBack%d/%s%s.local\n',...
                    'TargetPatchSizeXandY %d,%d\n', ...
                    'MinFidsTotalAndEachSurface %d,%d\n',...
                    'MinSizeOrOverlapXandY 0.5,0.5\n',...
                    'LocalOutputOptions 1,0,1\n', ...
                    'EOF'],mbOUT{1:3},fullTiltSizeXandY,...
                           fullPixelSize,fullPixelSize,...
                           mbOUT{1:3},outCTF,mbOUT{1:3},outCTF,mbOUT{1:3},outCTF,...
                           mbOUT{1:3},outCTF,mbOUT{1:3},outCTF,mbOUT{1:3},outCTF, ...
                           iRawTltName,tiltAliOption(1:2),...
                           10 / sqrt(nFidsTotal),tiltAliOption(3:4),flgLocalMag, ...
                           mbOUT{1:3},outCTF,targetPatchSize,targetPatchSize,...
                           nFiducialsPerPatch,floor(nFiducialsPerPatch/3));
                         
% % % Assume that any backlash was solved well enough that there are no major
% % % discontinuities in the coarse alignment. Mag and rot are solved/ tilt in
% % % the global solution anyhow, so this shouldn't be a bit deal.                         
% % % 'SeparateGroup	1-%d\n',... 
% % % iViewGroup,
%        fprintf(aliCom,'\n\ngrep -A %d  " At minimum tilt" ./mapBack%d/%s_ta.log >  ./mapBack%d/tmp.log',nPrjs+2,mbOUT{1:3},mbOUT{1:3});
%        fprintf(aliCom,'\nawk ''{if(NR >3) print $5}'' ./mapBack%d/tmp.log > mapBack%d/%s.mag',mbOUT{1:3},mbOUT{1:3});

        fclose(aliCom);
        system(sprintf('chmod a=wrx %smapBack%d/%s.align',mbOUT{1:3}));
        
        if (iTiltSeries == 1)
          if ( flgAltRun )
            fOUT = fopen(sprintf('%smapBack%d/runAlignments_alt.sh',mbOUT{1:2}),'w');
          else
            fOUT = fopen(sprintf('%smapBack%d/runAlignments.sh',mbOUT{1:2}),'w');
          end
          fprintf(fOUT,['#!/bin/bash\n\n%smapBack%d/%s.align > ',...
                        '%smapBack%d/%s.align_ta.log &\n'],mbOUT{1:3},mbOUT{1:3});
          % Since we send to the background in a shell, makes sure the
          % function waits on children.
          %if (iTiltSeries == nTiltSeries)
          %  fprintf(fOUT,'\nwait\n');
          %end
          fclose(fOUT);
        else
          if ( flgAltRun )
            fOUT = fopen(sprintf('%smapBack%d/runAlignments_alt.sh',mbOUT{1:2}),'a');
          else
            fOUT = fopen(sprintf('%smapBack%d/runAlignments.sh',mbOUT{1:2}), 'a');
          end
          fprintf(fOUT,['%smapBack%d/%s.align > ',...
                        '%smapBack%d/%s.align_ta.log &\n'], ...
                        mbOUT{1:3},mbOUT{1:3});  
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

end % loop over tomos

if ~( flgAltRun )
 % fOUT = fopen(sprintf('%smapBack%d/runAlignments.sh', ...
  %                               mbOUT{1:2}), 'a');
  altFile =  sprintf('%smapBack%d/runAlignments_alt.sh',mbOUT{1:2});
  mainFile = sprintf('%smapBack%d/runAlignments.sh',mbOUT{1:2});
  
  if exist(altFile,'file')
    fprintf('Combining Results from alt and main\n');
    system(sprintf('cat %s >> %s',altFile,mainFile));
  end

  fOUT = fopen(mainFile,'a');

  fprintf(fOUT,'\nwait\n');
  fclose(fOUT);
  system(sprintf('chmod a=wrx %smapBack%d/runAlignments.sh', mbOUT{1:2}));
  if (flgRunAlignments)
    system(sprintf('%smapBack%d/runAlignments.sh', mbOUT{1:2}));
  end

end 
clear refODD refEVE
if ( conserveDiskSpace )
  system(sprintf('rm %smapBack%d/%s_mapBack.st', mbOUT{1:3}));
end



if (tmpCache)
  if ( flgAltRun )
    system(sprintf('mv %smapBack%d mapBack%d_alt', mbOUT{1:2}, mbOUT{2}));
  else
    system(sprintf('mv %smapBack%d mapBack%d', mbOUT{1:2}, mbOUT{2}));
  end
end

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
% Since we've updated (potentially) mapBackRePrjSize, save the new metaData.
subTomoMeta.currentTomoCPR =  subTomoMeta.currentTomoCPR + 1;
if ~( flgAltRun )
  save(pBH.('subTomoMeta'), 'subTomoMeta');
end

end

