function [  ] = BH_ctf_Correct3d( PARAMETER_FILE, varargin )
%Apply a full 3d CTF ala Jensen and Kornberg approach
%   Detailed explanation goes here

% Read in 2dCtf stacks to trouble shoot
PosControl2d=0;
pBH = BH_parseParameterFile(PARAMETER_FILE);

try
  tiltWeight = pBH.('ctf_tiltWeight');
catch
  tiltWeight = [0.2,0];
end

try 
  % -1, whiten before ctf, 1 whiten after - test both.
  flgWhitenPS = pBH.('whitenPS');
catch
  flgWhitenPS = 0;
end

try
  % Shift the defocus origin from the sample origin in the microscope, to
  % the plane defined by the center of mass of subtomograms (in Z);
  shiftDefocusOrigin = pBH.('shiftDefocusToSubTomoCOM');
catch
  % This seems to be working well, so set default to true
  shiftDefocusOrigin = 1;
end
%default to cycle number zero for
%determining mean z height of particles
if nargin > 1
  % Start part way through on a seperate node
  tiltStart = str2double(varargin{1});
  fprintf('\n\nYou have asked to start with tilt-series # %d\n',tiltStart);
  fprintf('The cycle number is now automatically determined, no need to enter\n\n');
  if tiltStart == 0
    error('tiltStart of zero asked for!')
  end
else
  tiltStart = 1;
end

%cycleNumber = sprintf('cycle%0.3d',CYCLE);
%fprintf('cycle is %d\n',CYCLE);


fprintf('tiltweight is %f %f\n',tiltWeight);




tmpCache= pBH.('fastScratchDisk');

% Check to make sure it even exists
if isempty(dir(tmpCache))
    fprintf('\n\nIt appears your fastScratchDisk\n\t%s\ndoes not exist!\n\n',tmpCache);
    tmpCache = '';
end

reconScaling = 1;

if isempty(tmpCache)
  tmpCache='cache'; 
  flgCleanCache = 0;
  CWD='';
else
  flgCleanCache = 1;
  CWD = sprintf('%s/',pwd);
  % Check for a trailing slash
  slashCheck = strsplit(tmpCache,'/');
  if isempty(slashCheck{end})
    % This means the final character was a slash, strip it
    tmpCache = sprintf('%scache',tmpCache); %strjoin(slashCheck(1:end-1),'/');
  else
    tmpCache = sprintf('%s/cache',tmpCache); 
  end   
end
fprintf('tmpCache is %s\n',tmpCache);
system(sprintf('mkdir -p %s',tmpCache));

load(sprintf('%s.mat', pBH.('subTomoMeta')), 'subTomoMeta');
mapBackIter = subTomoMeta.currentTomoCPR;
masterTM = subTomoMeta; clear subTomoMeta

CYCLE = masterTM.currentCycle;
cycleNumber = sprintf('cycle%0.3d',CYCLE);
fprintf(' %s \n',cycleNumber);


try
  usePreCombDefocus = pBH.('usePreCombDefocus')
catch
  usePreCombDefocus = 0
end

try
  flgDampenAliasedFrequencies = pBH.('flgDampenAliasedFrequencies')
catch
  flgDampenAliasedFrequencies = 0;
end

try
  flg2dCTF = pBH.('flg2dCTF');
catch
  flg2dCTF = 0;
end   
  

% ctf3dDepth=pBH.('defocusErrorEst')
%mean in case cones.
resTarget = mean(masterTM.('currentResForDefocusError')*.5)                                  

%%%%% Take these from param file later.
% % % % % samplingRate = 2; nGPUs = 1; ctf3dDepth = 50e-9; usableArea = [3710,3710,1000]./samplingRate;
samplingRate = pBH.('Ali_samplingRate')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Need to change this to
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% either be precomputed, or
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% get x,y from TLT 20,21 and
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% max Z from
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% recon/name.coords

nGPUs = pBH.('nGPUs');
% Optionally specify gpu idxs
if numel(nGPUs) == 1
  gpuList = 1:nGPUs;
else
  gpuList = nGPUs;
  nGPUs = length(gpuList);
end

pixelSize = pBH.('PIXEL_SIZE').*10^10 .* samplingRate;
if pBH.('SuperResolution')
  pixelSize = pixelSize * 2;
end


tiltList_tmp = fieldnames(masterTM.mapBackGeometry);
tiltList_tmp = tiltList_tmp(~ismember(tiltList_tmp,{'viewGroups','tomoName'}));
nST = 1; tiltList = {};
for iStack = 1:length(tiltList_tmp)
  if masterTM.mapBackGeometry.(tiltList_tmp{iStack}).nTomos
    tiltList{nST} = tiltList_tmp{iStack};
    nST = nST +1;
  end
end
clear tiltList_tmp
  
tomoList = fieldnames(masterTM.mapBackGeometry.tomoName);

nTilts = length(tiltList);




% Divide the tilt series up over each gpu
iterList = cell(nGPUs,1);
for iGPU = 1:nGPUs
  iterList{iGPU} = iGPU+(tiltStart-1):nGPUs:nTilts;
  iterList{iGPU}
end

try
  parpool(nGPUs) 
catch
  delete(gcp)
  parpool(nGPUs)
end

% Instead of trying to do some imod_wait situation, just loop over the
% tilts in serial prior to entering the parallel reconstruction to make sure the
% prebinned stacks are there.

for iGPU = 1:nGPUs
  for iTilt = iterList{iGPU}
    
    nTomos = masterTM.mapBackGeometry.(tiltList{iTilt}).nTomos;
    iCoords = masterTM.mapBackGeometry.(tiltList{iTilt}).coords ./samplingRate;
    iCoords(:,1:4) = fix(iCoords(:,1:4));
    iTomoList = cell(nTomos,1);
    
    % For now, since the tilt geometry is not necessarily updated (it is manual)
    % in the masterTM, check that newer (possible perTilt refined) data is
    % not present.
    try
      % make sure there isn't a refined version first.
      TLTNAME = sprintf('fixedStacks/ctf/%s_ali%d_ctf_refine.tlt',tiltList{iTilt},mapBackIter+1);
      TLT = load(TLTNAME);
      fprintf('using refined TLT %s\n', TLTNAME);
    catch
      TLTNAME = sprintf('fixedStacks/ctf/%s_ali%d_ctf.tlt',tiltList{iTilt},mapBackIter+1);
      TLT = load(TLTNAME);
      fprintf('using TLT %s\n', TLTNAME);
    end
    
       
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
     checkRecon = sprintf('cache/%s_%d_bin%d.rec', ...
                          tiltList{iTilt},iTomo,samplingRate);
     if exist(checkRecon, 'file')
       fprintf('found %s to already exits\n',checkRecon);
       alreadyMade = alreadyMade +1;
     end
      
    end

    if alreadyMade == nTomos
      fprintf('All tomos 1-%d found to exist for tilt-series %s\n',nTomos,tiltList{iTilt});
      continue
    end
  
  
  
    preBinStacks(TLT, tiltList{iTilt}, mapBackIter,1,...
                                                      samplingRate,...
                                                      PosControl2d,...
                                                      tiltWeight);

  end                                                 
end

% All data is handled through disk i/o so everything unique created in the 
% parfor is also destroyed there as well.
parfor iGPU = 1:nGPUs
% for iGPU = 1:nGPUs 
  gpuDevice(gpuList(iGPU));
  % Loop over each tilt 
  for iTilt = iterList{iGPU}
% % % % %   for iTilt = 6 
    nTomos = masterTM.mapBackGeometry.(tiltList{iTilt}).nTomos;
    iCoords = masterTM.mapBackGeometry.(tiltList{iTilt}).coords ./samplingRate;
    iCoords(:,1:4) = fix(iCoords(:,1:4));
    iTomoList = cell(nTomos,1);
    
    % For now, since the tilt geometry is not necessarily updated (it is manual)
    % in the masterTM, check that newer (possible perTilt refined) data is
    % not present.
    try
      % make sure there isn't a refined version first.
      TLTNAME = sprintf('fixedStacks/ctf/%s_ali%d_ctf_refine.tlt',tiltList{iTilt},mapBackIter+1);
      TLT = load(TLTNAME);
      fprintf('using refined TLT %s\n', TLTNAME);
    catch
      TLTNAME = sprintf('fixedStacks/ctf/%s_ali%d_ctf.tlt',tiltList{iTilt},mapBackIter+1);
      TLT = load(TLTNAME);
      fprintf('using TLT %s\n', TLTNAME);
    end
    
    

    
   
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
     checkRecon = sprintf('cache/%s_%d_bin%d.rec', ...
                          tiltList{iTilt},iTomo,samplingRate);
     if exist(checkRecon, 'file')
       fprintf('found %s to already exits\n',checkRecon);
       alreadyMade = alreadyMade +1;
     end
      
    end

    if alreadyMade == nTomos
      fprintf('All tomos 1-%d found to exist for tilt-series %s\n',nTomos,tiltList{iTilt});
      continue
    end
    
    [ avgZ, maxZ, tomoNumber ] = calcAvgZ(masterTM,iCoords,tiltList{iTilt}, ...
                                          iTomoList,nTomos, pixelSize, ...
                                          samplingRate, cycleNumber);
    
    if ( shiftDefocusOrigin )
      fprintf('Using avgZ %3.3e nm as the defocus origin\n',avgZ*10^9);
    else
      avgZ = 0;
      fprintf('Using sample origin as the defocus origin\n');
      fprintf('If you want to use the COM of subTomos, set shiftDefocusToSubTomoCOM=1\n');
    end
    % Load in the aligned stack, check for outliers, and mask usable area with
    % a soft edge, also mask with high-pass filter.
    [ maskedStack, ~, ~ ] = loadAndMaskStack(TLT,tiltList{iTilt},...
                                                         mapBackIter,...
                                                         maxZ*10/pixelSize,...
                                                         samplingRate,...
                                                         PosControl2d,...
                                                         tiltWeight,flgWhitenPS,pixelSize);
      
    [d1, d2, nPrjs] = size(maskedStack);
  

                                        %0;%mean(iCoords(:,6))*pixelSize*10^-10;
    
    if ( flg2dCTF )
      nSections = 1;
      ctf3dDepth = maxZ * 10 ^ -9;
    else
      dampeningMax = 0.90;                                  
      [ ctf3dDepth ] = BH_ctfCalcError( samplingRate*mean(TLT(:,16)), ...
                TLT(1,17),TLT(1,18),TLT(1,15), ...
                                              2048, TLT(1,19), ...
                                              resTarget,maxZ*10, ...
                                              dampeningMax,CYCLE);
      fprintf('\n\nUsing a ctfDepth of %2.2f nm for %s\n\n',ctf3dDepth*10^9,tiltList{iTilt});
      % sections centered at 0, which for now is also supposed to coincide with
      % the mean defocus determination, although this could be corrected using
      % knowledge of particle positions given assurance that particles are the
      % primary source of signal (and not carbon for example).
      nSections = ceil(maxZ/(ctf3dDepth*10^9));
      % max odd number
      nSections = nSections + ~mod(nSections,2);
    end
    fprintf('with %3.3f nm sections, correcting %d tilt-series\n',ctf3dDepth*10^9,nSections);
    
    % For each tomo create a list of slices that are to be reconstructed 
    % for every section section.
    
    [ sectionList ] = calcTomoSections(iCoords, tomoNumber,pixelSize, ...
                                        nSections,tiltList{iTilt}, ctf3dDepth);
      

    
   
    % Correct a tilt series for earch section which requires writing each to
    % disk for use of IMOD.
    if (mapBackIter)
      tiltErrorFile = sprintf('mapBack%d/%s_ali%d_ctf.beamTiltError', ...
                                     mapBackIter,tiltList{iTilt}, mapBackIter);
      try      
        tiltError = load(tiltErrorFile); 
        if numel(tiltError) ~= 1
          error('tiltError should be a single number in degrees.\n');
        else
          fprintf('\nUsing %f degrees for beam tilt error.\n',tiltError)
        end        
        
      catch
        fprintf('\nTiltErrorFile %s not found.\n', tiltErrorFile);
        fprintf('\nUsing 0 degrees for beam tilt error.\n')
        tiltError = 0;
      end
      

    else
      fprintf('\nUsing 0 degrees for beam tilt error b/c no polishing yet.\n')
      tiltError = 0;
    end
    
    for iSection = 1:nSections
    

      defFitFull = '';
      preCombDefocus = 0;
      if (mapBackIter)
        defFitFull = sprintf('mapBack%d/%s_ali%d_ctf.defFidFull',mapBackIter, ...
                                                       tiltList{iTilt},mapBackIter);          
        if exist(defFitFull,'file')
          preCombDefocus = load(defFitFull);
          fprintf('3dCTF using pre calc combined per tilt defocus %s\n',defFitFull);
        else
          fprintf('Did not find %s\n!!',defFitFull);
        end
      end
      
      if ~(usePreCombDefocus)
        preCombDefocus = 0;
      end

      if (PosControl2d)
        correctedStack = maskedStack;
      else

      [ correctedStack ] = ctfMultiply_tilt(nSections,iSection,ctf3dDepth, ...
                                            avgZ,TLT,pixelSize,maskedStack,...
                                            maxZ*10/pixelSize,flgDampenAliasedFrequencies,...
                                            preCombDefocus,samplingRate);  
      end
      % Write out the stack to the cache directory as a tmp file
     
      outputStack = sprintf('%s/%s_ali%d_%d.fixed', ...
                            tmpCache,tiltList{iTilt},mapBackIter+1,iSection)
      SAVE_IMG(MRCImage(gather(correctedStack)),outputStack,pixelSize);
      correctedStack = [];
   
      % Loop over tomos reconstructing section and appending a file to 
      for iT = 1:nTomos
        iTomo = tomoNumber(iT);   
        if any(sectionList{iT}(iSection,:)+9999)
          


          reconName = sprintf('%s/%s_ali%d_%d_%d.rec', ...
                              tmpCache,tiltList{iTilt},mapBackIter+1,iTomo,iSection);

          if (mapBackIter)
            rawTLT = sprintf('%smapBack%d/%s_ali%d_ctf.tlt',CWD,mapBackIter,tiltList{iTilt},...
                                                       mapBackIter);
            LOCAL = sprintf('%smapBack%d/%s_ali%d_ctf.local',CWD,mapBackIter,tiltList{iTilt}, ...
                                                           mapBackIter);
                                          
                   
          else 
            rawTLT = sprintf('%sfixedStacks/%s.tlt',CWD,tiltList{iTilt});
            LOCAL = sprintf('%sfixedStacks/%s.local',CWD,tiltList{iTilt});
          end
          
          % Put a local copy if using a nondefault cache
%           if ( flgCleanCache )
%          sprintf('cp %s/%s %s/%s',CWD,rawTLT,tmpCache,rawTLT)
%          sprintf('cp %s/%s %s/%s',CWD,LOCAL,tmpCache,LOCAL)
%             system(sprintf('cp %s/%s %s/%s',CWD,rawTLT,tmpCache,rawTLT));
%             system(sprintf('cp %s/%s %s/%s',CWD,LOCAL,tmpCache,LOCAL));
%             rawTLT = sprintf('%s/%s',tmpCache,rawTLT)
%             LOCAL = sprintf('%s/%s',tmpCache,LOCAL)
%             
%           end
            

          fprintf('Local file %s\n',LOCAL);
          
          if exist(LOCAL,'file')
            flgLocal = 1;
          else
            fprintf('Did not find local alignment information at %s\n',LOCAL);
            flgLocal = 0;
          end

          % hangover from slab padding, remove later.
          padRec = 0;
          
          

          nTiltWorkers = 4;
          nTotalSlices = (iCoords(iTomo,3)-iCoords(iTomo,2)+1)
          tiltChunkSize = ceil(nTotalSlices/nTiltWorkers)
          tiltChunks = iCoords(iTomo,2):tiltChunkSize:iCoords(iTomo,3)
          tiltChunks(end) = iCoords(iTomo,3)
          totalSlices = [tiltChunks(1),tiltChunks(end)]

          rCMD = sprintf(['tilt -input %s -output %s.TMPPAD -TILTFILE %s -UseGPU %d ', ...
                       '-WIDTH %d -COSINTERP 0 -THICKNESS %d -SHIFT %f,%f '],...
                       outputStack, reconName, rawTLT, gpuList(iGPU), ...
                       iCoords(iTomo,1),floor(sectionList{iT}(iSection,5))+2*padRec,...
                       iCoords(iTomo,5),sectionList{iT}(iSection,6));

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
          [recError] = system(sprintf('%s.sh  > /dev/null ',reconName)); % /dev/null
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
      iTomo = tomoNumber(iT);
      reconNameFull = sprintf('cache/%s_%d_bin%d.rec', ...
                              tiltList{iTilt},iTomo,samplingRate)
                            
      recCMD = 'newstack -fromone';
      for iSection = 1:nSections
        reconName = sprintf('%s/%s_ali%d_%d_%d.rec', ...
                             tmpCache, tiltList{iTilt},mapBackIter+1,iTomo,iSection)
        
        if any(sectionList{iT}(iSection,:)+9999)  
          recCMD = [recCMD,sprintf(' -secs 1-%d %s', ...
                                   floor(sectionList{iT}(iSection,5)), ...
                                   reconName)];
        else
          
          fprintf('no info for section %d for tomo %d\n',iSection,iTomo); 
        end
      end
      
      system([recCMD, sprintf(' %s > /dev/null ',reconNameFull)]); %/dev/null
      
     
      for iSection = 1:nSections
        cleanUp3 = sprintf('rm %s/%s_ali%d_%d_%d.rec', ...
                         tmpCache,tiltList{iTilt},mapBackIter+1,iTomo,iSection);
        system(cleanUp3);
        cleanUp4 = sprintf('rm %s/%s_ali%d_%d_%d.rec.sh', ...
                         tmpCache,tiltList{iTilt},mapBackIter+1,iTomo,iSection);
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

function [STACK, evalMask, deltaZ] = loadAndMaskStack(TLT, STACK_PRFX, ...
                                                      mapBackIter,maxZpix,...
                                                      samplingRate,...
                                                      PosControl2d,...
                                                      tiltWeight,flgWhitenPS,pixelSize)

if (PosControl2d) 
  prefix = 'ctf';
  suffix = '_ctf'
else
  prefix = 'ali';
  suffix = '';
end

if samplingRate > 1
  fullStack = sprintf('%sStacks/%s_ali%d%s.fixed', ...
                       prefix,STACK_PRFX,mapBackIter+1,suffix);
  inputStack = sprintf('cache/%s_ali%d%s_bin%d.fixed',...
                       STACK_PRFX,mapBackIter+1,suffix,samplingRate);
  if ~exist(inputStack, 'file')
%     binCMD = sprintf('newstack -bin %d -antialias 6 %s %s > /dev/null',samplingRate,fullStack,inputStack);
% %     binCMD = sprintf('newstack -bin %d -antialias 6 %s %s ',samplingRate,fullStack,inputStack);
% 
%     system(binCMD);
    BH_multi_loadOrBin(fullStack,-1.*samplingRate,2);
   
  end
else
  inputStack = sprintf('%sStacks/%s_ali%d%s.fixed',...
                        prefix,STACK_PRFX,mapBackIter+1,suffix)
end

system(sprintf('header %s',inputStack));
% iHeader = MRCImage(inputStack,0);
% STACK = gpuArray(single(getVolume(iHeader)));

STACK = single(getVolume(MRCImage(inputStack)));

% iHeader = getHeader(iHeader);
% iPixelHeader = [iHeader.cellDimensionX/iHeader.nX, ...
%                 iHeader.cellDimensionY/iHeader.nY, ...
%                 iHeader.cellDimensionZ/iHeader.nZ];


[d1,d2,d3] = size(STACK);
nPrjs = d3;

useableArea = [d1-128,d2-128,maxZpix];
  

  % Keeping all on the GPU won't fit if a full 4k4k is reconstucted. run and then pull.
   % This should be double checked because it is annoying.
  %[exposureFilter] = gather((BH_exposureFilter([d1,d2], TLT,'GPU',samplingRate,0)));

  
  [evalMask, deltaZ ] = BH_multi_projectionMask( [d1,d2,d3;useableArea], TLT, 'cpu' );
  
 
  
 
  % Local normalization doesn't address any large scale gradients in the
  % images. Do a simple high pass over the lowest 7 frequencyBinns
  bandNyquist = BH_bandpass3d([d1,d2,1],0,0,1,'GPU','nyquistHigh');

  taperMask = gpuArray(fspecial('gaussian',[9,9],1.5));
  

  for iPrj = 1:nPrjs


    iEvalMask = gpuArray(evalMask(:,:,TLT(iPrj,1)));
    %iExpFilter = gpuArray(exposureFilter(:,:,TLT(iPrj,1)));
    iProjection = gpuArray(STACK(:,:,TLT(iPrj,1)));
    
%     iMask = convn(single(iEvalMask),taperMask,'same');


     iProjection = iProjection - mean(iProjection(iEvalMask));
     iProjection  = real(ifftn(fftn(iProjection).*bandNyquist));

%      iProjection  = real(ifftn(fftn(iProjection.*iMask).*bandNyquist));


    inFin = ~(isfinite(iProjection)); nInf = sum(inFin(:));
    if (nInf)
    % fprintf('Removing %d (%2.4f) inf from prj %d\n',nInf,100*nInf/numel(iProjection),TLT(iPrj,1));
     iProjection(inFin) = 0;
    end

    iRms = rms(iProjection(iEvalMask));
    outliers = (iProjection > 6 * iRms); nOutliers = sum(outliers(:));
    tiltScale = 1- ( abs(sind(TLT(iPrj,4))).* tiltWeight(1));
    if (nOutliers)

%       iProjection(outliers) = sign(iProjection(outliers)).*3.*iRms.*((rand(size(iProjection(outliers)))./2)+0.5);
      iProjection(outliers) = 6.*iRms.* (rand(size(iProjection(outliers)))-0.5);
      iProjection = iProjection ./ ( rms(iProjection(iEvalMask)) ./ tiltScale);
    else
      iProjection = iProjection ./ ( iRms ./ tiltScale);
    end

    if ( flgWhitenPS )
      %fprintf('confirm whitening PS\n.'); 
      [iProjection,~] = BH_whitenNoiseSpectrum(iProjection,'',pixelSize,1);
    end

    if tiltWeight(2)
      % I don't think this makes sense, but test keeping the power constant
      % after application of the exposure filter.
      iProjection = fftn(iProjection);
      iPower = sum(abs(iProjection(:)));
      iProjection = iProjection .* iExpFilter;
      STACK(:,:,TLT(iPrj,1)) = gather(single(real(ifftn(iProjection.* ...
                                    (iPower./sum(abs(iProjection(:))))))));
    else
      STACK(:,:,TLT(iPrj,1)) = gather(iProjection);%gather(single(real(ifftn(fftn(iProjection) .* iExpFilter))));
    end
%     STACK(:,:,TLT(iPrj,1)) = gather(single(iMask.*real(ifftn(fftn(iProjection) .* iExpFilter))));
  %   STACK(:,:,TLT(iPrj,1)) = gather(single(iProjection.*iMask));
    clear iProjection iMask iExpFilter iEvalMask 
  end


  clear bandNyquist iMask exposureFilter  iProjection lowRMSMAsk
% % Push to gpu when initializing workers
% deltaZ = gather(deltaZ);
% evalMask = gather(evalMask);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = preBinStacks(TLT, STACK_PRFX, mapBackIter,usableArea,...
                                                      samplingRate,...
                                                      PosControl2d,...
                                                      tiltWeight)
                                                    
                                                   

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
                     prefix,STACK_PRFX,mapBackIter+1,suffix,samplingRate);
if ~exist(inputStack, 'file')
  BH_multi_loadOrBin(fullStack,-1.*samplingRate,2);
end


end


function  [ sectionList ] = calcTomoSections(iCoords, tomoNumber, pixelSize,...
                                              nSections,tiltName, ctf3Depth)

nTomos = length(tomoNumber);
sectionList = cell(nTomos,1);
for iTomo = 1:nTomos
  % min and max in absolute pixels min and max from 1:nZrecon
  sectionList{iTomo} = zeros(nSections,6);
end

% With rounding this could end up a bit short except the top and bottom are both
% half a section larger than minimally needed.
nSec = floor(ctf3Depth*10^10/pixelSize)      ;
nSec = nSec + ~mod(nSec,2);
halfSec = (nSec-1)/2;

for iT = 1:length(tomoNumber)
  iTomo = tomoNumber(iT);
  % Origin + originshift
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
                                              avgZ,TLT,pixelSize,maskedStack,...
                                              maxZ,flgDampenAliasedFrequencies,...
                                              preCombDefocus,samplingRate)
% Correct in strips which is more expensive but (hopefully) more accurate.  



[d1,d2,nPrjs] = size(maskedStack)
PIXEL_SIZE = pixelSize*10^-10
% This is just going to be written out to disk so keep in main memory.
correctedStack = zeros(d1,d2,nPrjs,'single');

% This should probably come from the calcSections. If avgZ has not been
% forced to zero, then assume the defocus determined is the distance from
% the center of mass of subtomograms in Z to the focal plane, rather than
% the center of mass of the tomograms (specimen)

defocusOffset = (((nSections-1)/-2+(iSection-1))*ctf3dDepth);
fprintf('using offset %3.3e for section %d with COM offset %3.3e\n',defocusOffset,iSection,avgZ);
defocusOffset = defocusOffset - avgZ;

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
radialGrid = {radialGrid./PIXEL_SIZE,0,phi};
phi = [];

%[exposureFilter] = ((BH_exposureFilter([d1,d2], TLT,'GPU',samplingRate,0)));

for iPrj = 1:nPrjs

  maxEval = cosd(TLT(iPrj,4)).*(d1/2) + maxZ./2*abs(sind(TLT(iPrj,4)));
  oX = ceil((d1+1)./2);
  iEvalMask = floor(oX-maxEval):ceil(oX+maxEval);
  iExposureFilter = BH_exposureFilter(fastFTSize,TLT(iPrj,:),'GPU',samplingRate,0);
% % %    iEvalMask = gpuArray(evalMask(:,:,TLT(iPrj,1)) );
% % %    iDeltaZ = gpuArray(deltaZ(:,:,TLT(iPrj,1)));
  
%   STRIPWIDTH = 2*floor(((tileSize/pixelSize)*abs(cosd(TLT(iPrj,4))).^1.5/2));
  STRIPWIDTH = min(floor((0.5*ctf3dDepth/PIXEL_SIZE)/abs(tand(TLT(iPrj,4)))),512);
  STRIPWIDTH = STRIPWIDTH + mod(STRIPWIDTH,2);
  % take at least 1200 Ang & include the taper if equal to STRIPWIDTH
  tileSize   = floor(max(600./pixelSize, STRIPWIDTH + 28));
  tileSize = tileSize + mod(tileSize,2);
  %fprintf('stripwidth tilesize %d %d\n',STRIPWIDTH,tileSize);
  incLow = ceil(tileSize./2);
  incTop = tileSize - incLow;
  border = ceil(incLow+apoSize)+1;
  
  ddF = TLT(iPrj,12);
  dPhi = TLT(iPrj,13);
  D0 = TLT(iPrj,15);
  %TLT(iPrj,16); 
 

  padVal = BH_multi_padVal([d1,d2],fastFTSize);
  trimVal = BH_multi_padVal(fastFTSize,[d1,d2]);
  

  iProjection = BH_padZeros3d(maskedStack(:,:,TLT(iPrj,1)),padVal(1,:),padVal(2,:),'GPU','singleTaper');
  iProjectionFT = fftn(iProjection).*iExposureFilter; clear iExposureFilter
  correctedPrj = zeros([d1,d2],'single','gpuArray');
%   tiltedDefocusOffset = defocusOffset/cosd(TLT(iPrj,4));
  tiltedDefocusOffset = defocusOffset.*cosd(TLT(iPrj,4));

  envFilter = gpuArray(fspecial('gaussian',[12,12],9));
  envDiff = 60e-9.*[1,1,0];
  
  if any( preCombDefocus )
   wrkFids = ( preCombDefocus(:,4) == TLT(iPrj,1) - 1 );
   wrkDefocus = preCombDefocus(wrkFids,[2,3,5]);
   wrkDefocus(:,1:2) = wrkDefocus(:,1:2)./samplingRate;
  end
  stripDefocusOffset = floor(STRIPWIDTH/2);
  
  for i = 1: STRIPWIDTH : d1
    
    if (i+tileSize-1) < d1
        endIDX = (i+tileSize-1);
        endCUT = i + STRIPWIDTH - 1 + 7;  
        trimmedSIZE = STRIPWIDTH;
    elseif any(ismember(i:endIDX,iEvalMask))
        endIDX = d1;
        endCUT = d1;
        trimmedSIZE = endCUT-i+1 -7;
    end    
    
    % The eval mask condition can be replaced once the per tomo condition
    % is trusted.
% % %     if any( iEvalMask(i:endIDX,floor(d2/2)) )
    if any(ismember(i:endIDX,iEvalMask))
      %tile = iProjection(i:endIDX,:);

      iDeltaZ = (i + stripDefocusOffset - oX)*PIXEL_SIZE*-1.*tand(TLT(iPrj,4));
      avgDefX = [];
      if any(preCombDefocus) %%%%%
        stripIDX = wrkDefocus(:,1) > i & wrkDefocus(:,1) <= endIDX;
        if ( sum(stripIDX) )
          % If there are no particles in the strip, leave avgDefX as empty
          % and don't correct.
          avgDefX = mean(wrkDefocus(stripIDX,3));
          fprintf('Using %2.6e from comb rather than %2.6e\n',avgDefX, ...
            D0+ tiltedDefocusOffset + iDeltaZ);
% % %             D0+ tiltedDefocusOffset + (double(iDeltaZ(i+stripDefocusOffset,floor(d2/2))) .* PIXEL_SIZE));
              
        end
      end
      
        
      
      %if size(tile,1) < 14
        % At the end of the projection a very narrow strip can be cut out,
        % it is probably better just to skip this, but for now at least
        % don't try to taper, which would cause a crash since it is smaller
        % than the taper size.
     %   tile = BH_padZeros3d(tile, padVal(1,:), padVal(2,:),'GPU', 'single',mean(tile(:))); 
     % else
      %  tile = BH_padZeros3d(tile, padVal(1,:), padVal(2,:),'GPU', 'singleTaper',mean(tile(:)));
      %end


      if any(preCombDefocus) %%%%%
        DF = avgDefX + tiltedDefocusOffset;
      else
% % %         DF = D0 + tiltedDefocusOffset + (double(iDeltaZ(i,floor(d2/2))) .* PIXEL_SIZE);
        DF = D0 + tiltedDefocusOffset + iDeltaZ;

      end
      
      if ~( isempty(DF) )
        
        iDefocus = [DF - ddF, DF + ddF, dPhi];
        %fprintf('correcting strip centered on %d for prj %d with %3.3f %3.3f %3.3f',i,TLT(iPrj,1),iDefocus);
        % Assume reconstruction is padded by ~ 10 % on both ends.
  % % %       Hqz = BH_ctfCalc(PIXEL_SIZE,Cs,WAVELENGTH,iDefocus,[padTileSize,d2],AMPCONT,-1.0, 0.8*maxZ*10^-9);
  
         if PIXEL_SIZE < 2.0e-10
           % use double precision - this is not enabled, but needs to be -
           % requires changes to radial grid as well.
           
           Hqz = BH_ctfCalc(radialGrid,Cs,WAVELENGTH,iDefocus,fastFTSize,AMPCONT,-1,-1);
         else
           Hqz = BH_ctfCalc(radialGrid,Cs,WAVELENGTH,iDefocus,fastFTSize,AMPCONT,-1);
         end
         
        if ( flgDampenAlias )
          % Try this out with real space convolution, which should of
          % course be slow, but easy and direct.
          % The premise is that aliasing will be bad when the rings are too
          % close, and this can be dampened simply by averaging amplitudes
          % over a small window.
          dampSize = abs(floor(DF*10^6)); % Use the size of defocus to dermine severity.
          dampWindow = ones([1,1].*dampSize,'single','gpuArray')./dampSize.^2;
          Hqz = ifftshift(convn(Hqz,dampWindow,'same'));
        end
        
        tile  = iProjectionFT.*Hqz;

        
        tile = BH_padZeros3d(real(ifftn(tile)), ...
                             trimVal(1,:),trimVal(2,:),'GPU','single');

        % trim prior to pulling off gpu to minimize xfer
      else
        
        % No particles in this strip, so just replace with simple inversion
        % to keep the global image statistics ~ correct.

        tile = -1.*BH_padZeros3d(iProjection, trimVal(1,:),trimVal(2,:),'GPU','single');
        
      end
       
        %correctedStack(i + 7 : endCUT,:,TLT(iPrj,1)) = ...
         %             gather(tile(8:trimmedSIZE+7,:));
                    
      correctedPrj(i:endIDX,:) = tile(i:endIDX,:); clear tile
        
    else
    %fprintf('ignoring strip centered on %d for prj %d',i,TLT(iPrj,1));
    end
  end
 correctedStack(:,:,TLT(iPrj,1)) = gather(correctedPrj); clear correctedPrj
 clear iProjection iProjectionFT
end % end loop over projections
clear tile Hqz
end

function [avgZ, maxZ, tomoNumber] = calcAvgZ(masterTM,iCoords, ...
                                             tiltName,tomoList,...
                                             nTomos, pixelSize,...
                                             samplingRate,cycleNumber)

% Calculate the maximum extensions in Z and then how many separate sections
% need to be corrected.

maxZ = 0;
tomoNumber = zeros(nTomos,1);
for iTomo = 1:nTomos
  % The tomograms may not be listed monotonically so explicitly get their
  % id number
  tomoNumber(iTomo) = masterTM.mapBackGeometry.tomoName.(tomoList{iTomo}).tomoNumber;     
  nZdZ = iCoords(tomoNumber(iTomo),[4,6]);

  % half the size in z plus the shift back to the microscope coords.
  sZneeded = 2.*ceil(nZdZ(1)/2+abs(nZdZ(2))+1);
  if sZneeded > maxZ
    maxZ = sZneeded;
  end
end

maxZ = maxZ + (samplingRate*2);

maxZ = maxZ.*pixelSize./10;
fprintf('combining thickness and shift on tilt %s, found a maxZ  %3.3f nm\n',tiltName,maxZ);

% For now use cycle000, if adding a refinment focused on a specific set of
% particles, then consider that later.

try
  initGeom = masterTM.(cycleNumber).Avg_geometry;
catch
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

for iT = 1:nTomos
  iTomo = tomoNumber(iT);
  % Already scaled to sampled pixels
  tomoOrigin = ceil((iCoords(iTomo,4)+1)/2);
  micOrigin = iCoords(iTomo,6);
   
  iTomoName = sprintf('%s_%d',tiltName,iTomo);
  % shouldn't be any removed particles at this stage but later there would be.
  zList = initGeom.(iTomoName)(initGeom.(iTomoName)(:,26)~=-9999,13)./samplingRate;
  
  % shift from lower left to centered and include the tomos offset from the
  % microscope frame
  zList = zList - tomoOrigin + micOrigin;
  nSubTomos = nSubTomos + length(zList);
  totalZ = totalZ + sum(zList);
  fprintf('%s tomo has %d subTomos with mean Z %3.3f nm\n', ...
          iTomoName, length(zList), mean(zList)*pixelSize./10);
  clear zList
end

avgZ = totalZ/nSubTomos*pixelSize/10*10^-9;

fprintf('%s tilt-series has %d subTomos with mean Z %3.3f nm\n', ...
        tiltName, nSubTomos,avgZ*10^9);



end
