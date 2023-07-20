 function [  ] = BH_alignRaw3d(PARAMETER_FILE, CYCLE, varargin)
                                                               
%Extract and align class averages and references from 4D montages derived.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Goals & Limitations:
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   TODO
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global bh_global_print_shifts_in_particle_basis;
if isempty(bh_global_print_shifts_in_particle_basis)
  bh_global_print_shifts_in_particle_basis = true;
end

global bh_global_zero_lag_score;
if isempty(bh_global_zero_lag_score)
  bh_global_zero_lag_score = false
end

if (nargin ~= 2 && nargin ~= 3)
  error('args = PARAMETER_FILE, CYCLE, [1,abs(ccc),2,weighted,3,abs(weighted)]')
else
  parentFunc = mfilename;
  resumeVars = struct();
end

if nargin == 3
  
  flgWeightCCC = EMC_str2double(varargin{1});
else
  % default to linear ccc (which is actually weighted by the SNR though)
  flgWeightCCC = 0;
end

% Explicit reference to location of variables in main memory, or on the GPU. As
% in pcaPub, looking ahead to re-write in c++ for cuda, no cells allowed.
cpuVar = struct();
GPUVar = struct();

startTime =  clock;
CYCLE = str2num(CYCLE);
cycle_numerator = '';
cycle_denominator ='';
  flgStartThird = 0;
  flgReverseOrder = 0;
if numel(CYCLE) == 3
  cycle_numerator = CYCLE(2);
  cycle_denominator = CYCLE(3);
  CYCLE = CYCLE(1);
  flgStartThird = true;
elseif CYCLE < 0
  % Simple option to process in reverse order so that the load can be run on two
  % physically distinct systems at once.
  flgReverseOrder = 1;
  flgStartThird = 0;
  CYCLE = abs(CYCLE);


end



pBH = BH_parseParameterFile(PARAMETER_FILE);
cycleNumber = sprintf('cycle%0.3u', CYCLE);
load(sprintf('%s.mat', pBH.('subTomoMeta')), 'subTomoMeta');
mapBackIter = subTomoMeta.currentTomoCPR; 
reconScaling = 1;
try
  nPeaks = pBH.('nPeaks');
catch
  nPeaks = 1;
end

try
  flgCutOutVolumes=pBH.('flgCutOutVolumes')
catch
  flgCutOutVolumes=0
end

% TODO decide on a "reasonable" padding based on expected shifts.
try
  CUTPADDING = subTomoMeta.('CUTPADDING')
catch
  CUTPADDING=20
end

maxGoldStandard = subTomoMeta.('maxGoldStandard');


nGPUs = pBH.('nGPUs')


flgClassify= pBH.('flgClassify');
try
  flgMultiRefAlignment=pBH.('flgMultiRefAlignment');
catch
  flgMultiRefAlignment = 0;
end
try
  flgCenterRefCOM = pBH.('flgCenterRefCOM');
catch
  flgCenterRefCOM = 1;
end

try
  flgSymmetrizeSubTomos = pBH.('flgSymmetrizeSubTomos');
catch 
  flgSymmetrizeSubTomos = 0;
end
flgRaw_shapeMask =  0;%= pBH.('experimentalOpts')(3)
samplingRate = pBH.('Ali_samplingRate');

pixelSize      = pBH.('PIXEL_SIZE').*10^10.*samplingRate;
if pBH.('SuperResolution')
  pixelSize = pixelSize * 2;
end

flgPrecision = 'single'; %pBH.('flgPrecision');
angleSearch  = pBH.('Raw_angleSearch');
peakSearch   = (pBH.('particleRadius')./pixelSize);
peakCOM      = [1,1,1].*3;
className    = pBH.('Raw_className');  

try 
  loadTomo = pBH.('loadTomo')
catch
  loadTomo = 0;
end
try 
  eraseMaskType = pBH.('Peak_mType');
	eraseMaskRadius = pBH.('Peak_mRadius')./pixelSize;
  fprintf('Further restricting peak search to radius %f %f %f\n',...
          eraseMaskRadius);
  eraseMask = 1;
catch
  eraseMask = 0;
  fprintf('Using particle radius for peak search\n');
end

rotConvention = 'Bah';
% Check and override the rotational convention to get helical averaging.
% Replaces the former hack of adding a fifth dummy value to the angular search
try
  doHelical = pBH.('doHelical');
catch
  doHelical = 0;
end
if ( doHelical )
  rotConvention = 'Helical'
end

rotConvention

bFactor = pBH.('Fsc_bfactor');
if length(bFactor) > 1
  fprintf('multiple bFactors specified, using the first for alignment.\n');
  bFactor = bFactor(1);
end

try
  scaleCalcSize = pBH.('scaleCalcSize');
catch
  scaleCalcSize = 1.5;
end
% % % % if (flgClassify || flgMultiRefAlignment)
if (flgClassify)
  refName      = pBH.('Ref_className');
else
  refName = pBH.('Raw_className');
end
 
outputPrefix = sprintf('%s_%s', cycleNumber, pBH.('subTomoMeta'));



classVector{1}  = pBH.('Raw_classes_odd')(1,:);
classSymmetry{1}= pBH.('Raw_classes_odd')(2,:);


classVector{2}  = pBH.('Raw_classes_eve')(1,:);
classSymmetry{2}= pBH.('Raw_classes_eve')(2,:);


% % % % if (flgClassify || flgMultiRefAlignment)
if (flgClassify)
  geometry = subTomoMeta.(cycleNumber).ClassAlignment;
  refVectorFull{1}= [pBH.('Ref_references_odd');1]
  refVectorFull{2}= [pBH.('Ref_references_eve');1]
elseif (flgMultiRefAlignment)
  geometry = subTomoMeta.(cycleNumber).ClusterRefGeom;
  refVectorFull{1}= [pBH.('Raw_classes_odd');classVector{1} ]
  refVectorFull{2}= [pBH.('Raw_classes_eve');classVector{2} ]
else
  geometry = subTomoMeta.(cycleNumber).Avg_geometry;
  refVectorFull{1} = [pBH.('Raw_classes_odd');1];
  refVectorFull{2} = [pBH.('Raw_classes_eve');1];
end


% % % pathList= subTomoMeta.mapPath;
% % % extList = subTomoMeta.mapExt;
masterTM = subTomoMeta; clear subTomoMeta




refVector = cell(2,1);
refGroup = cell(2,1);
refSym = cell(2,1);

for iGold = 1:2
  % Sort low to high, because order is rearranged as such unstack
  refVectorFull{iGold} = sortrows(refVectorFull{iGold}', 1)';
  % class id corresponding to membership in ???_refName
  refVector{iGold} = refVectorFull{iGold}(1,:)
  % reference id, so multiple classes can be merged into one
  refGroup{iGold}  = refVectorFull{iGold}(3,:)
  % axial symmetry to apply, negative value indicates creating a mirrored ref
  % accros the corresponding axis
  refSym{iGold}    = refVectorFull{iGold}(2,:)
end

% make sure the number of references match the unique groups in the classVector
% and also that the class/group pairs match the class/ref pairs.
nReferences(1:2) = [length(unique(refGroup{1})),length(unique(refGroup{1}))];
nReferences = nReferences .* [~isempty(refGroup{1}),~isempty(refGroup{2})]


nRefOut(1:2) = [length(unique(refGroup{1})) + sum(( refSym{1} < 0 )),...
                length(unique(refGroup{2})) + sum(( refSym{2} < 0 ))];


%%%%%%%%%%%%%%%%%%%%%%%

% Get the number of tomograms to process.
tomoList = fieldnames(geometry);
nTomograms = length(tomoList);
tiltList = masterTM.tiltGeometry;
ctfGroupList = masterTM.('ctfGroupSize');



% mask defines area for angular search, peakRADIUS restricts translational


[ maskType, maskSize, maskRadius, maskCenter ] = ...
                                  BH_multi_maskCheck(pBH, 'Ali', pixelSize)

[ sizeWindow, sizeCalc, sizeMask, padWindow, padCalc ] = ...
                                       BH_multi_validArea( maskSize, maskRadius, scaleCalcSize  )


try 
  flgLimitToOneProcess = pBH.('flgLimitToOneProcess');
catch
  flgLimitToOneProcess = 0;
end

if ( loadTomo )
  limitToOne = loadTomo;
  if (flgLimitToOneProcess)
    limitToOne = min(limitToOne, flgLimitToOneProcess);
  end
elseif (flgLimitToOneProcess)
  limitToOne = flgLimitToOneProcess;
else
  limitToOne = pBH.('nCpuCores');
end

[ nParProcesses, iterList] = BH_multi_parallelJobs(nTomograms,nGPUs, sizeCalc(1),limitToOne);                                   
if ( flgReverseOrder )
  % Flip the order for reverse processing on a second machine. This will also disable saving of 
  % of the metadata so there aren't conflicts.
  for iParProc = 1:nParProcesses
    iterList{iParProc} = flip(iterList{iParProc});
  end
  
elseif ( flgStartThird )

  % Shift to start at one third through to process on a third machine. This will also disable saving of 
  % of the metadata so there aren't conflicts.
  for iParProc = 1:nParProcesses
    % Note the use of floor is more like ceiling here (rounds away from
    % zero)
    nParts = ceil(length(iterList{iParProc}) ./ cycle_denominator);
    fIDX = 1+(cycle_numerator - 1)*nParts;
    lIDX = min(cycle_numerator*nParts,length(iterList{iParProc}));
    iterList{iParProc} = iterList{iParProc}(fIDX:lIDX);
  end
  
end

if any(peakSearch > maskRadius)
  fprintf('\n\n\tpeakRADIUS should be <= maskRADIUS!!\n\n')
  peakSearch( (peakSearch > maskRadius) ) = ...
                                        maskRadius( (peakSearch > maskRadius) );
end

if ( any(classSymmetry{1}~=1) || any(classSymmetry{2}~=1) ) && flgSymmetrizeSubTomos
  flgSymmetry = true
else
  flgSymmetry = false
end
  
% Read in the references.
% Read in the references.
refIMG = cell(2,1);
refWGT = cell(2,1);
refWgtROT = cell(2,1);
imgCounts = cell(2,1);
for iGold = 1:2

  if iGold == 1
    halfSet = 'ODD';
  else
    halfSet = 'EVE';
  end


  imgNAME = sprintf('class_%d_Locations_REF_%s', refName, halfSet)  

  
  weightNAME = sprintf('class_%d_Locations_REF_%s_Wgt', refName, halfSet);
  imgCounts{iGold} = masterTM.(cycleNumber).(imgNAME){3};


  [ refTMP ] = BH_unStackMontage4d(1:nReferences(iGold), ...
                                   masterTM.(cycleNumber).(imgNAME){1}, ...
                                   masterTM.(cycleNumber).(imgNAME){2},...
                                   sizeWindow);
                                 
  [ wdgTMP ] = BH_unStackMontage4d(1:nReferences(iGold), ...
                                masterTM.(cycleNumber).(weightNAME){1},...
                                masterTM.(cycleNumber).(weightNAME){2},...
                                sizeCalc);

  sizeREF = masterTM.(cycleNumber).(imgNAME){2}{1}(2:2:6)';
                               
  if (flgCenterRefCOM)
% % % % % % %    [ comMask ] = BH_mask3d(maskType, sizeMask, maskRadius, maskCenter);
   [ comMask ]  = EMC_maskShape(maskType, sizeMask, maskRadius, 'gpu', {'shift', maskCenter});
  end
                            
  % get boxSize 
  n = 1 ; tIMG = cell(numel(refVector{iGold})); tWDG = cell(numel(refVector{iGold}));tWDG_r = tWDG;
  for iP = 1:numel(refTMP)
    if ~isempty(refTMP{iP})
      tIMG{n} = refTMP{iP}; refTMP{iP} = [];
      if (flgCenterRefCOM)
        % Not sure if this is always the best approach, but it may be
        % useful in some cases.
% % % % % % %         [~,iCOM] = BH_mask3d(gpuArray(tIMG{n}).*comMask,pixelSize,'','',1);
        
        [~, ~, ~,iCOM] = EMC_maskReference(gpuArray(tIMG{n}).*comMask, pixelSize, {'fsc',true; 'com', true});
        fprintf('centering ref %d on COM %3.3f %3.3f %3.3f \n',n,iCOM);

        tIMG{n} = BH_resample3d(tIMG{n},[0,0,0],gather(iCOM), ...
                                           {'Bah',1,'spline'},'cpu','inv');
      
      end
      tWDG{n} = wdgTMP{iP}; wdgTMP{iP} = [];
      tWDG{n} = tWDG{n} - min(tWDG{n}(:)) + 1e-6;
      tWDG{n} = tWDG{n} ./ max(tWDG{n}(:));
      n = n + 1;
    end
  end

  wdgPAD = BH_multi_padVal(size(tWDG{1}), sizeCalc)
  for iWdg = 1:n-1
    tWDG_r{iWdg} = BH_padZeros3d(tWDG{iWdg},wdgPAD(1,:),wdgPAD(2,:),...
                                                           'cpu',flgPrecision);
    tWDG{iWdg} = ifftshift(tWDG_r{iWdg});
  end
  
  refIMG{iGold} = tIMG ; clear tIMG refTMP
  refWGT{iGold} = tWDG; clear tWDG wdgTMP
  refWgtROT{iGold} = tWDG_r; clear tWDG_r

  clear comMask
  
end

[ refIMG ] = BH_multi_combineLowResInfo( refIMG, imgCounts, pixelSize, maxGoldStandard );



% optimize the fft for the given size. Padding to the next power of 2 is usually
% slower given the dimensionalityl of the volume data.
fftPlanner = rand(sizeCalc);
fftw('planner', 'exhaustive');
fftn(fftPlanner);
clear fftPlanner

% Make a mask, and apply to the average motif && save a masked,
% binned copy of the average for inspection. 
% % % 
% % %   [ volMask ] = gather(BH_mask3d(maskType, sizeMask, maskRadius, maskCenter));
% % %     volBinary = (volMask >= 0.01);


% In principle the window and mask could be different sizes, however, I
% think I am currently forcing them to be the same.

  
  % make rotationally invariant
% % % % % % %   [ peakMask] = gather(BH_mask3d('sphere', sizeWindow, [1,1,1].*max(peakSearch), maskCenter));
  [ peakMask ]  = gather(EMC_maskShape('sphere', sizeWindow, [1,1,1].*floor(max(peakSearch)), 'gpu', {'shift', maskCenter}));

  if (eraseMask)
    % Mask could be smaller than the normal taper would allow, so instead
    % of thresholding a normal mask, take this alt route.
    eraseMask = ones(ceil(2.*eraseMaskRadius),'single');
    padEraseMask = BH_multi_padVal(size(eraseMask),sizeCalc);
    eraseMask = BH_padZeros3d(eraseMask,padEraseMask(1,:),padEraseMask(2,:),'cpu','single');
    eraseMask = single(find(eraseMask < 1));
  else
    eraseMask = [];
  end
  
   
  if ( flgRaw_shapeMask )
% % % % % % %     [ volMask ] = BH_mask3d(maskType, sizeWindow, maskRadius, maskCenter);
    [ volMask ] = EMC_maskShape(maskType, sizeWindow, maskRadius, 'gpu', {'shift', maskCenter});

    % Currently not set up for mult-ref alignment
    iRef = 1;
    % Use the geometric mean so that excluded areas mask out
% % % % % % %     [ volMask ] = gather(sqrt(volMask .* ...
% % % % % % %               BH_mask3d(refIMG{1}{iRef}+refIMG{2}{iRef},pixelSize,'','')));
 
    [ volMask ] = gather(sqrt(volMask .* ...            
              EMC_maskReference(refIMG{1}{iRef}+refIMG{2}{iRef}, pixelSize, {'fsc', true})));

  else
% % % % % % %     [ volMask ] = gather(BH_mask3d(maskType, sizeWindow, maskRadius, maskCenter));
    [ volMask ] = gather(EMC_maskShape(maskType, sizeWindow, maskRadius, 'gpu', {'shift', maskCenter}));
    
  end     
                                                   
% % %   [ peakMask] = gather(BH_mask3d(maskType, sizeMask, peakSearch, maskCenter));
% % %     peakBinary = (peakMask >= 0.01);   
    
%   [ refInterp] = gather(BH_mask3d(maskType, sizeREF, peakSearch, maskCenter));
%     refInterp = (refInterp >= 0.01);   

  bandpassFilt = cell(nReferences(1),1);
  bandpassFiltREF = bandpassFilt;
  wCCC  = cell(nReferences(1),1);
  for iWccc = 1:length(nReferences(1));
    wCCC{iWccc} = 0;
  end
  if (flgClassify || flgMultiRefAlignment)
    for iRef = 1:nReferences(1)
      if (flgClassify)
        fscINFO = masterTM.(cycleNumber).('fitFSC').(sprintf('REF%d',iRef));
      else
        fscINFO = masterTM.(cycleNumber).('fitFSC').(sprintf('Raw%d',iRef)); % % % %
      end

      [radialGrid,~,~,~,~,~ ] = BH_multi_gridCoordinates(sizeCalc, 'Cartesian', ...
                                                         'GPU', {'none'}, 1, 0, 1 );
      radialGrid = single(radialGrid./pixelSize);
      % returns a cpu array
      if (flgWeightCCC)
        [ bandpassFilt{iRef}, ~,wCCC] =  BH_multi_cRef( fscINFO, radialGrid, bFactor, 1, 1); 
      else
        [ bandpassFilt{iRef}, ~] =  BH_multi_cRef( fscINFO, radialGrid, bFactor, 1); 
      end
      

      bandpassFiltREF{iRef} = 1;

    end
  else


    
    for iRef = 1
      fscINFO = masterTM.(cycleNumber).('fitFSC').('Raw1');
      [radialGrid,~,~,~,~,~ ] = BH_multi_gridCoordinates(sizeCalc, 'Cartesian', ...
                                                         'GPU', {'none'}, 1, 0, 1 );
      radialGrid = single(radialGrid./pixelSize);
      % returns a cpu array
      if (flgWeightCCC)
        [ bandpassFilt{iRef},~,wCCC{iRef} ] =  BH_multi_cRef( fscINFO, radialGrid, bFactor, 1, 1 ); 
      else
        [ bandpassFilt{iRef},~ ] =  BH_multi_cRef( fscINFO, radialGrid, bFactor, 1 ); 
      end

      bandpassFiltREF{iRef} = 1;

      
    end

  end
  
%   if (flgWeightCCC)
%     for i = 1:length(wCCC{1})
%       i
%       length(wCCC{1}{i})
%     end
%   end
    
  % This is just used to limit the interpolation search so use the most
  % permissive bandpass, while the appropriate bandpass (given a multi-ref
  % alignment) will still be applied.
  mostPermissive = zeros(1,nReferences(1));
  for iRef = 1:nReferences(1)
    mostPermissive(iRef) = sum(bandpassFilt{iRef}(:));
  end
  [~,mPidx] = max(mostPermissive);
  
 wdgBinary = single(find(fftshift(bandpassFilt{mPidx} > 10^-2)));
                                        

ref_FT1 = cell(2,1);
ref_FT2 = cell(2,1);

  
for iGold = 1:2

  if iGold == 1
    halfSet = 'ODD';
  else
    halfSet = 'EVE';
  end

  nOut = 1;
  refOUT = cell(2.*nReferences(iGold),2);

  for iRef = 1:nReferences(iGold)
    
    refTMP_2 = refIMG{iGold}{iRef}; refIMG{iGold}{iRef} = [];
    refTMP = refTMP_2(padWindow(1,1) + 1: end - padWindow(2,1), ...
                      padWindow(1,2) + 1: end - padWindow(2,2), ...
                      padWindow(1,3) + 1: end - padWindow(2,3));


    % if not using a weighted average (adapted SPW filter), apply an
    % approximation the cRef from Rosenthal/Henderson. This is currently always set to one
    % and is just doing the masking and normalization. It should be okay to just apply the mask
    % and rely on the normalization during the CCC calc. TODO
      ref_FT1{iGold}{iRef} = gather(conj(BH_bandLimitCenterNormalize(...
                             refTMP.*volMask, bandpassFiltREF{iRef}, (volMask>0.01), padCalc, flgPrecision)));

     
                           

    ref_FT2{iGold}{iRef} = gather(refTMP_2);
    % Trim for output reference
        refTMP_2 = refTMP_2(padWindow(1,1) + 1: end - padWindow(2,1), ...
                      padWindow(1,2) + 1: end - padWindow(2,2), ...
                      padWindow(1,3) + 1: end - padWindow(2,3));
    
    % Overwrite a copy of the filtered, bandpassed ref for output
    refOUT{nOut} = real(ifftn(conj(ref_FT1{iGold}{iRef})));
    refOUT{nOut} = gather(refOUT{nOut}(padCalc(1,1) + 1: end - padCalc(2,1), ...
                          padCalc(1,2) + 1: end - padCalc(2,2), ...
                          padCalc(1,3) + 1: end - padCalc(2,3)) .* volMask);
                       

    
    refOUT{nOut} = refOUT{nOut}.*volMask;
    
    refOUT{nOut+1} = real(ifftn(BH_bandLimitCenterNormalize(...
                          refTMP_2.*peakMask, '', (peakMask > 0.01), padCalc, 'single')));
    refOUT{nOut+1} = gather(refOUT{nOut+1}(padCalc(1,1) + 1: end - padCalc(2,1), ...
                            padCalc(1,2) + 1: end - padCalc(2,2), ...
                            padCalc(1,3) + 1: end - padCalc(2,3)) .* peakMask );  
    nOut = nOut + 2;
    
    refOUT{nOut} = refOUT{nOut} - mean(refOUT{nOut}(:));
    refOUT{nOut} = refOUT{nOut} ./ rms(refOUT{nOut}(:));
    
    refOUT{nOut+1} = refOUT{nOut+1} - mean(refOUT{nOut+1}(:));
    refOUT{nOut+1} = refOUT{nOut+1} ./ rms(refOUT{nOut+1}(:));
  end
 

  % Save a montage of the masked reference & shape masks if requested.

  maskedOUTFILE = sprintf('%s_maskedRef-mont_%s.mrc',outputPrefix,halfSet);
  [ maskedReferences, ~ ] = BH_montage4d(refOUT, '');
  SAVE_IMG(MRCImage(single(maskedReferences)), maskedOUTFILE);
  

end

clear refIMG refWDG refOUT iRef

%%%%%%%%%%%%%%%%%%%%% Determine the angular search, if any are zero, don't
%%%%%%%%%%%%%%%%%%%%% search at all in that dimension.

[  nInPlane, inPlaneSearch, angleStep, nAngles] ...
                                      = BH_multi_gridSearchAngles(angleSearch)

[masterTM] = BH_recordAngularSampling( masterTM, cycleNumber, angleStep, inPlaneSearch);                               
                                    
% set truth value for refinement during out of plane search

if any(angleStep(:,1))
  flgRefine = true;
  fprintf('flgRefine set to %s','True');
else
  flgRefine = false;
  fprintf('flgRefine set to %s','False');
end

angleStep(:,1)
any(angleStep(:,1))
nCount = 1;

firstLoop = true;
nIgnored = 0;

bestAnglesResults = cell(nParProcesses,1);
geometryResults   = cell(nParProcesses,1);



try
  EMC_parpool(nParProcesses+1)
catch
  delete(gcp('nocreate'))
  EMC_parpool(nParProcesses+1)
end

size(ref_FT2)

system('mkdir -p alignResume');

system(sprintf('mkdir -p alignResume/%s',outputPrefix));
softenWeight = 1/sqrt(samplingRate);
for iParProc = 1:nParProcesses

  % Caclulating weights takes up a lot of memory, so do all that are necessary
  % prior to the main loop -- CHANGE THE CHECK TO JUST READ THE HEADER NOT LOAD
  % THE WEIGHT INTO GPU MEMORY
  iParProc
  iterList{iParProc}
  for iTomo = iterList{iParProc}
    
    BH_multi_loadOrCalcWeight(masterTM,ctfGroupList,tomoList{iTomo},samplingRate ,...
                              sizeCalc,geometry,flgPrecision,1);
                            
     
  end
end
  
% Clear all of the GPUs prior to entering the main processing loop
for iGPU = 1:nGPUs
  g = gpuDevice(iGPU);
  fprintf('\n\nClear gpu %d mem prior to main loop, %3.3e available\n\n',iGPU,g.AvailableMemory);
  clear g
end 

parVect = 1:nParProcesses;
parfor iParProc = parVect
%for iParProc = 1:nParProcesses
%profile on
  bestAngles_tmp = struct();
  geometry_tmp = geometry;

% % %     % Get the gpuIDX assigned to this process
% % %     iGPUidx = gpuDevice();
% % %     iGPUidx = iGPUidx.Index;
    gpuIDXList = mod(parVect+nGPUs,nGPUs)+1;
    iGPUidx = gpuIDXList(iParProc);
    gpuDevice(iGPUidx);
    fprintf('parProc %d/%d assigned to GPU %d\n',iParProc,nParProcesses,iGPUidx);
  for iTomo = iterList{iParProc}


    
  nCtfGroups = ctfGroupList.(tomoList{iTomo})(1);
  % Check for interupted alignment.
  previousAlignment = sprintf('alignResume/%s/%s.txt',outputPrefix,tomoList{iTomo});
  if exist(previousAlignment,'file')
    % Sometimes when multiple nodes are used, an extra line is added.
    % TODO fix this workaround
    system(sprintf('awk ''{if($10 != "") print $0 }'' %s > %s_clean; mv %s_clean %s',...
          previousAlignment,previousAlignment,previousAlignment,previousAlignment));
    bestAngles_tmp.(tomoList{iTomo}) = load(previousAlignment);
    fprintf('Using existing alignment info for %s\n', tomoList{iTomo});
  else
    % There is some memory leak somewhere that I haven't been able to figure
    % out. I am clearing all vars but output in the children functions ... this
    % isn't ideal, but for now is an acceptable stop gap.
    %D = gpuDevice(gpuList(iGPU));
   
    % shake up the random number generator for phi and theta
    rng('shuffle');
    
    bandpassFilt_tmp = cell(nReferences(1),1);
    bandpassFiltREF_tmp = cell(nReferences(1),1);
    for iRef = 1:nReferences(1)
      if flgMultiRefAlignment <= 2
      bandpassFilt_tmp{iRef} = gpuArray(bandpassFilt{iRef});
      bandpassFiltREF_tmp{iRef} = gpuArray(bandpassFiltREF{iRef});
      else
      bandpassFilt_tmp{iRef} = (bandpassFilt{iRef});
      bandpassFiltREF_tmp{iRef} = (bandpassFiltREF{iRef});
      end
    end
    
    
    
    ref_FT1_tmp = cell(2,1);
    ref_FT2_tmp = cell(2,1);
    ref_WGT_tmp = cell(2,1);
    ref_WGT_rot = cell(2,1);
    
    
    volMask_tmp = gpuArray(volMask);
    volBinary_tmp = single(find( volMask_tmp > 0.01 ));
    peakMask_tmp = gpuArray(peakMask);
    peakBinary_tmp = single(find( peakMask_tmp > 0.01 )); 
    wdgBinary_tmp = gpuArray(wdgBinary);
    eraseMask_tmp = gpuArray(eraseMask);
    
    wCCC_tmp = cell(length(wCCC));
    

    
    for iRef = 1:nReferences(1)
      for iWccc = 1:length(wCCC{iRef})
        if (flgWeightCCC)
          wCCC_tmp{iRef}{iWccc} = gpuArray(wCCC{iRef}{iWccc});
        else
          % The check in xcf_rotational looks for a cell 
          wCCC_tmp{iRef} = 0;
        end
      end
    end
      
    
          
    for iGold = 1:2
      for iRef = 1:nReferences(iGold)
        if flgMultiRefAlignment <= 2
        ref_FT1_tmp{iGold}{iRef} = gpuArray(ref_FT1{iGold}{iRef});
        ref_FT2_tmp{iGold}{iRef} = gpuArray(ref_FT2{iGold}{iRef});
        ref_WGT_tmp{iGold}{iRef} = gpuArray(refWGT{iGold}{iRef});
        ref_WGT_rot{iGold}{iRef} = gpuArray(refWgtROT{iGold}{iRef});
        else
        % Temp workaround, six big ribo refs crashing
        ref_FT1_tmp{iGold}{iRef} = (ref_FT1{iGold}{iRef});
        ref_FT2_tmp{iGold}{iRef} = (ref_FT2{iGold}{iRef});
        ref_WGT_tmp{iGold}{iRef} = (refWGT{iGold}{iRef});
        ref_WGT_rot{iGold}{iRef} = (refWgtROT{iGold}{iRef});
        end
      end
    end   
   

    sprintf('\nWorking on %d/%d volumes',iTomo,nTomograms)
    tic;

    % Load the tomo into gpu
    tomoName = tomoList{iTomo};
      %fprintf('gpu %d working on tomoName %s\n', iGPU, tomoName);

      tiltGeometry = masterTM.tiltGeometry.(tomoList{iTomo});
    % Load in the geometry for the tomogram, and get number of subTomos.
    positionList = geometry_tmp.(tomoList{iTomo});
    
    tomoNumber = masterTM.mapBackGeometry.tomoName.(tomoList{iTomo}).tomoNumber;
    tiltName   = masterTM.mapBackGeometry.tomoName.(tomoList{iTomo}).tiltName;
    coords = masterTM.mapBackGeometry.(tiltName).coords(tomoNumber,1:4);
    
%     [ binShift, ~ ] = BH_multi_calcBinShift( coords, samplingRate); 
    binShift = [0,0,0];
    nSubTomos = size(positionList,1);



      iTiltName = masterTM.mapBackGeometry.tomoName.(tomoName).tiltName;
      wgtName = sprintf('cache/%s_bin%d.wgt',iTiltName,samplingRate);       
%         wgtName = sprintf('cache/%s_bin%d.wgt', tomoList{iTomo},...
%                                                 samplingRate);
        maxWedgeMask = BH_unStackMontage4d(1:nCtfGroups,wgtName,...
                                          ceil(sqrt(nCtfGroups)).*[1,1],'');
        maxWedgeIfft = maxWedgeMask;

        for iWdg = 1:length(maxWedgeMask)
          if ~isempty(maxWedgeMask{iWdg})           
            maxWedgeIfft{iWdg} = ifftshift(maxWedgeIfft{iWdg}.^softenWeight);
            maxWedgeMask{iWdg} = maxWedgeMask{iWdg}.^softenWeight;
          end
        end

      
      fprintf('loaded %s.\n',wgtName);

        
              % Can't clear inside the parfor, but make sure we don't have two tomograms
      % in memory at once.
     
     tomoNumber = masterTM.mapBackGeometry.tomoName.(tomoList{iTomo}).tomoNumber;
     tiltName = masterTM.mapBackGeometry.tomoName.(tomoList{iTomo}).tiltName;
     reconCoords = masterTM.mapBackGeometry.(tiltName).coords(tomoNumber,:);
     
     if (flgCutOutVolumes)
       volumeData = [];
     else
      [ volumeData, ~ ] = BH_multi_loadOrBuild( tomoList{iTomo}, ...
                                       reconCoords, mapBackIter, ...
                                       samplingRate,iGPUidx,reconScaling,loadTomo);  
       if ( loadTomo )
         volHeader = struct();
         volHeader.('nX') = size(volumeData,1);
         volHeader.('nY') = size(volumeData,2);
         volHeader.('nZ') = size(volumeData,3);
       else
         volHeader = getHeader(volumeData);  
       end
     end


    % For now, set up for full grid-search only, as I intend to just do
    % translational and in-plane searches for now anyhow.
  
    [~,iv1,iv2,iv3] = BH_resample3d(volMask_tmp,eye(3),[0,0,0],...
                                    {'Bah',1,'linear',1,volBinary_tmp}, ...
                                                         'GPU', 'inv');   
    inputVectors = {iv1,iv2,iv3};
    iv1 = []; iv2 = []; iv3 = [];    
    cccStorageBest = cell(nPeaks,1);
    cccStorageRefine = cell(nPeaks,1);
    for iPeak = 1:nPeaks
      cccStorageBest{iPeak} = zeros(nSubTomos,10); 
      cccStorageRefine{iPeak}= zeros(nSubTomos,10);
    end
    % reset for each tomogram
    wdgIDX = 0;
    
    for iSubTomo = 1:nSubTomos
      breakPeak = 0; % for try catch on cut out vols    
      if (wdgIDX ~= positionList(iSubTomo,9))
        % Geometry is sorted on this value so that tranfers are minimized,
        % as these can take up a lot of mem. For 9 ctf Groups on an 80s
        % ribo at 2 Ang/pix at full sampling ~ 2Gb eache.
        wdgIDX = positionList(iSubTomo,9);
        fprintf('pulling the wedge %d onto the GPU\n',wdgIDX);
        % Avoid temporar
        iMaxWedgeMask = []; iMaxWedgeIfft = [];
        iMaxWedgeMask = gpuArray(maxWedgeMask{wdgIDX});
        iMaxWedgeIfft = gpuArray(maxWedgeIfft{wdgIDX});       
      end
      

        [~,iw1,iw2,iw3] = BH_resample3d(iMaxWedgeMask, eye(3), [0,0,0], ...
                              {'Bah',1,'linear',1,wdgBinary_tmp}, ...
                                                         'GPU', 'inv');  
        inputWgtVectors = {iw1,iw2,iw3};
        iw1 = []; iw2 = []; iw3 = [];
                                                       
      for iPeak = 1:nPeaks
        if (breakPeak) 
          continue;
        end
        getInitialCCC = 1;
        cccInitial = zeros(nReferences(1),10,flgPrecision, 'gpuArray');
        cccStorage2= zeros(nAngles(1).*nReferences(1),10,'gpuArray');
        powerOut = zeros(nAngles(1).*nReferences(1),1,'gpuArray');
      
    
    
        % Used in refinment loop
        angCount = 1;
      
        % Check that the given subTomo is not to be ignored
        classIDX = positionList(iSubTomo, 26+26*(iPeak-1));
        particleIDX = positionList(iSubTomo, 4);
        iGold = positionList(iSubTomo, 7);
     

        if classVector{iGold}(1,:) == 0
          classPosition = 1;
          flgAllClasses = true;
        else
          classPosition = find(classVector{iGold}(1,:) == classIDX);
          flgAllClasses = false;
        end



        if (classIDX ~= -9999) && ... % All previously ignored particles
         ( flgAllClasses ||  ismember(classIDX, classVector{iGold}(1,:)) ) 

            
        center = positionList(iSubTomo,[11:13]+26*(iPeak-1))./samplingRate + binShift;
        angles = positionList(iSubTomo,[17:25]+26*(iPeak-1));
        
        % Find range to extract, and check for domain error.
       if (flgCutOutVolumes)
         % Need some check that the windowsize has not changed! TODO TODO
         
        [ indVAL, padVAL, shiftVAL ] = ...
                      BH_isWindowValid(2*CUTPADDING+sizeWindow, ...
                                       sizeWindow,maskRadius, center);
       else
        [ indVAL, padVAL, shiftVAL ] = ...
                      BH_isWindowValid([volHeader.nX,volHeader.nY,volHeader.nZ], ...
                                       sizeWindow,maskRadius, center);         
       end




      if ischar(indVAL)
        fprintf('\nnow ignoring particle %d from tomo %d', iSubTomo,iTomo)
        nIgnored = nIgnored + 1;
        geometry_tmp.(tomoList{iTomo})(iSubTomo, 26) = -9999;     
      else
        
        
        if (flgCutOutVolumes)
          % Test with some generic padding , only to be used on bin 1 at
          % first!!! TODO add a flag to check this.
          try
            particleOUT_name = sprintf('cache/subtomo_%0.7d_%d.mrc',positionList(iSubTomo,4),iPeak);
            iparticle = gpuArray(getVolume(MRCImage(particleOUT_name),[indVAL(1,1),indVAL(2,1)], ...
                                                                    [indVAL(1,2),indVAL(2,2)], ...
                                                                    [indVAL(1,3),indVAL(2,3)],'keep'));
          catch
            fprintf('\n\nDid not load cut out vol. on subTomo %d FixMEEEEEE\n\n',iSubTomo);
            geometry_tmp.(tomoList{iTomo})(iSubTomo, 26) = -9999;
            breakPeak = 1;
            continue;
          end
        else

          if ( loadTomo )
            iparticle = gpuArray(volumeData(indVAL(1,1):indVAL(2,1), ...
                                            indVAL(1,2):indVAL(2,2), ...
                                            indVAL(1,3):indVAL(2,3)));

          else
            iparticle = gpuArray(getVolume(volumeData,[indVAL(1,1),indVAL(2,1)], ...
                                                      [indVAL(1,2),indVAL(2,2)], ...
                                                      [indVAL(1,3),indVAL(2,3)],'keep'));
          end
          
        end
        [ iparticle ] = BH_padZeros3d(iparticle,  padVAL(1,1:3), ...
                                      padVAL(2,1:3), 'GPU', 'singleTaper');
         

        
        for iAngle = 1:size(angleStep,1)
     
          theta    = angleStep(iAngle,1);
          thetaInc = angleStep(iAngle,4);       
          % Calculate the increment in phi so that the azimuthal sampling is
          % consistent and equal to the out of plane increment.

          phiInc = angleStep(iAngle,3);

          % To prevent only searching the same increments each time in a limited
          % grid search, radomly offset the azimuthal angle by a random number
          % between 0 and 1/2 the azimuthal increment.
        
          azimuthalRandomizer = (rand(1)-0.5)*phiInc;

          for iAzimuth = 0:angleStep(iAngle,2)
            phi = rem((phiInc * iAzimuth)+azimuthalRandomizer,360);

           for iInPlane = inPlaneSearch
            psi    = iInPlane;
            psiInc = angleStep(iAngle,5);
            %[phi,theta,psi-phi];


            RotMat = BH_defineMatrix([phi, theta, psi - phi],rotConvention, 'inv');
            RotMat = reshape(angles,3,3) * RotMat;

            cccStorageTrans= zeros(1.*nReferences(1),10,'gpuArray');

            for alignLoop = 1:2


              switch alignLoop
                
                case 1
                  % This takes care of non-inter shift in the origin that is
                  % ignored during the windowing of the particle.
                  estPeakCoord = shiftVAL;
                  % Estimate the peakshift by rotating the ref not the particle.
                  iTrimParticle = ...
                         iparticle(padWindow(1,1) + 1:end - padWindow(2,1) , ...
                                   padWindow(1,2) + 1:end - padWindow(2,2) , ...
                                   padWindow(1,3) + 1:end - padWindow(2,3) );
                               

                case 2
  
                  
                  bestOfRefs = sortrows(gather(cccStorageTrans), -6);
                  %sortrows(gather(cccStorage1(angCount:angCount+nReferences(1)-1,:)),-6);

                  estPeakCoord = bestOfRefs(1,8:10);
                 

             

                % Assuming if class specific symmetry, then some not just 1
                if (flgSymmetry) 
                  symmetry = classSymmetry{iGold}(1, classPosition);
                  %fprintf('Symmetry confirmation %d\n',symmetry);
                    [ iTrimParticle ] = BH_resample3d(iparticle, RotMat,... 
                                                  estPeakCoord,...
                                                  {'Bah',symmetry,'linear',1,volBinary_tmp}, ...
                                                  'GPU', 'inv',inputVectors);



                     if (getInitialCCC)
                     [ iTrimInitial ] =  BH_resample3d(iparticle, ...
                                                    reshape(angles,3,3),... 
                                                    shiftVAL,...
                                                    {'Bah',symmetry,'linear',1,volBinary_tmp}, ...
                                                    'GPU', 'inv',inputVectors);

                    
           
% % %                      powerInitial =  sum(abs(iTrimInitial(volBinary_tmp))).^2;   
% % % 
                
                      iWedgeInitial = BH_resample3d(iMaxWedgeMask, reshape(angles,3,3), [0,0,0], ...
                                               {'Bah',symmetry,'linear',1,wdgBinary_tmp}, ...
                                               'GPU', 'inv',inputWgtVectors);
                      
                                 
                     end           

                        iWedgeMask = BH_resample3d(iMaxWedgeMask, RotMat, [0,0,0], ...
                                                 {'Bah',symmetry,'linear',1,wdgBinary_tmp}, ...
                                                 'GPU', 'inv',inputWgtVectors);
                    




                else
                  symmetry = 1;
                  % Transform the particle, and then trim to motif size

                  [ iTrimParticle ] = BH_resample3d(iparticle, RotMat, ...
                                            estPeakCoord, {'Bah',1,'linear',1,volBinary_tmp}, 'GPU', 'inv',inputVectors);
                       iWedgeMask = BH_resample3d(iMaxWedgeMask, RotMat, [0,0,0], ...
                                       {'Bah',1,'linear',1,wdgBinary_tmp}, 'GPU', 'inv',inputWgtVectors);



                  if (getInitialCCC)
                     
                     iTrimInitial =  BH_resample3d(iparticle, reshape(angles,3,3),... 
                                           shiftVAL,{'Bah',1,'linear',1,volBinary_tmp}, 'GPU', 'inv',inputVectors);
                                         

                                   
                     iWedgeInitial = BH_resample3d(iMaxWedgeMask, reshape(angles,3,3), [0,0,0], ...
                                                  {'Bah',1,'linear',1,wdgBinary_tmp}, 'GPU', 'inv',inputWgtVectors);         
                                                                          
                  end
                end % Symmetry or not + interpolation
                
                 
% % %                  powerOut(angCount) =   sum(abs(iTrimParticle(volBinary_tmp))).^2;
                 
              end

                
              switch flgMultiRefAlignment
                case 0
                  refToAlign = 1;
                case 1
                  refToAlign = 1:max(nReferences(:));
                case 2
                  refToAlign = classIDX;
                otherwise
                  error('flgMultiRefAlignment is not 0,1,2')
              end

              for iRef = refToAlign % 1:max(nReferences(:)) 

                switch alignLoop
                  
                  case 1

                  
                    % use transpose of RotMat
           
                    iRotRef = BH_resample3d(ref_FT2_tmp{iGold}{iRef}, RotMat', ...
                                            estPeakCoord, {'Bah',1,'linear',1,peakBinary_tmp}, 'GPU', 'forward',inputVectors);
                    
                    iRotWdg = BH_resample3d(ref_WGT_rot{iGold}{iRef}, RotMat', ...
                                            [0,0,0], {'Bah',1,'linear',1,wdgBinary_tmp}, 'GPU', 'forward',inputWgtVectors);
                                      
%                     iRotRef = ...
%                               iRotRef(padWindow(1,1) + 1:end - padWindow(2,1) , ...
%                                       padWindow(1,2) + 1:end - padWindow(2,2) , ...
%                                       padWindow(1,3) + 1:end - padWindow(2,3) );                                          

                    % maybe I should be rotating peak mask here in case it has
                    % an odd shape, since we are leaving the proper frame
                        
                    iRotRef = BH_bandLimitCenterNormalize(...
                                                         iRotRef.*peakMask_tmp,...
                                                         bandpassFiltREF_tmp{iRef} ,peakBinary_tmp,...
                                                         padCalc,flgPrecision);
                                                       
                    rotPart_FT = BH_bandLimitCenterNormalize(...
                                 iTrimParticle.*peakMask_tmp,...
                                 bandpassFilt_tmp{iRef} ,peakBinary_tmp,padCalc,flgPrecision); 
                               
                    [ peakCoord ] =  BH_multi_xcf_Translational( ...
                                                            rotPart_FT.*ifftshift(iRotWdg), ...
                                                            conj(iRotRef).*iMaxWedgeIfft,...
                                                            peakMask_tmp, peakCOM,eraseMask_tmp);


                    cccStorageTrans(iRef,:) = [iRef, particleIDX, ...
                                               phi, theta, psi - phi, ...
                                               0, 0, ...
                                               peakCoord + estPeakCoord];                                      
                  case 2
                    
                    % get starting point
                    if (getInitialCCC)
                      
                      initialRotPart_FT = BH_bandLimitCenterNormalize(...
                                              iTrimInitial.*volMask_tmp,...
                                              bandpassFilt_tmp{iRef} ,volBinary_tmp,padCalc,flgPrecision);
                    
                    

                      [ iCCC, ~ ] = ...
                                  BH_multi_xcf_Rotational( initialRotPart_FT, ...
                                                           ref_FT1_tmp{iGold}{iRef}, ...
                                                           ifftshift(iWedgeInitial),...
                                                           ref_WGT_tmp{iGold}{iRef}, ...
                                                           wCCC_tmp{iRef});


                          

                        cccInitial(iRef,:) = [iRef, particleIDX, ...
                                              0,0,0, ...
                                              iCCC, 1, ...
                                              shiftVAL];
                      
                   
                    initialRotPart_FT = [];
                    
                  
                    end
                    
                  rotPart_FT = BH_bandLimitCenterNormalize(...
                                            iTrimParticle.*volMask_tmp,...
                                            bandpassFilt_tmp{iRef} ,volBinary_tmp,padCalc,flgPrecision);
                                          
                                          
                     
                      
                      [ iCCC, ~ ] = ...
                                  BH_multi_xcf_Rotational( rotPart_FT, ...
                                                           ref_FT1_tmp{iGold}{iRef},...
                                                           ifftshift(iWedgeMask),...
                                                           ref_WGT_tmp{iGold}{iRef}, ...
                                                           wCCC_tmp{iRef});
                                                         



                                              

                    % Note that no new translational estimate is made, so no
                    % need to multiply by RotMat
                    cccStorage2(angCount,:) = ...
                                          [iRef, particleIDX, ...
                                           phi, theta, psi , ...
                                           iCCC, 1, ...
                                           estPeakCoord];                                         


                angCount = angCount + 1;                          
                end


              end % loop over references.
              

           
            end
            % This volume won't be needed until the next subTomo is considered, 
            % which is also where getInitialCCC Boolean is set to True again.
            iTrimInitial = [];
            getInitialCCC = 0;
            
            end % in plane angles
          end % azimuth
        end % polar

% % %         fprintf('Power ratio is  %3.3f\n',powerOut./powerInitial);

        cccPreRefineSort =  sortrows(gather(cccStorage2),-6);
     
        if (length(refToAlign) > 1)
          cccInitial = sortrows(gather(cccInitial), -6);
          cccInitial = cccInitial(1,:);
        else
          cccInitial = gather(cccInitial(refToAlign,:));
          
        end
        
        if cccInitial(1,6 ) > cccPreRefineSort(1,6)
          cccPreRefineSort(1,:) = cccInitial(1,:);
        end
        
       
        
        % This only seems to be a problem with cut out volumes.
        % Normalization maybe?
        if ~any(cccPreRefineSort(1,:))
          cccStorageBest{iPeak}(iSubTomo,:) = cccInitial(1,:);
          fprintf('all Zeros in PreRefine search, revert on subtomo %d peak %d\n',iSubTomo,iPeak);
          continue     
        end
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     

        if (flgRefine)
         
          

          % Get the results from just this subTomo and sort on CCC

          rRef  = cccPreRefineSort(1,1);
          rPart = cccPreRefineSort(1,2);
          rPhi  = cccPreRefineSort(1,3);
          rPhiInc = phiInc / 4;
          rTheta= cccPreRefineSort(1,4);
          rTheInc = thetaInc /2;
          rPsi  = cccPreRefineSort(1,5);
          rPsiInc = psiInc /2;
          % Confirm shiftVAL is doing what it should be
          rXYZest  = cccPreRefineSort(1,8:10);

          if (rTheInc)
            % For a larger out of plane step, search a larger range in plane
            psiRefineStep = floor(sqrt(rTheInc));
          else
            psiRefineStep = 1;
          end
          
          thetaRefineStep =1;
          phiRefineStep=2; 
          totalRefineStep = [psiRefineStep, thetaRefineStep, phiRefineStep];
          totalRefineStep = prod((2.*totalRefineStep)+1);
          
          cccStorage3 = zeros(totalRefineStep,10,'gpuArray');
          
          if (rPsiInc == 0)
            inPlaneRefine = rPsi - psiRefineStep*rTheInc./2:rTheInc./2: rPsi+psiRefineStep*rTheInc./2;
          else
            inPlaneRefine  = rPsi-  psiRefineStep*rPsiInc : rPsiInc : rPsi    + psiRefineStep*rPsiInc;
          end
          polarRefine    = rTheta-thetaRefineStep*rTheInc : rTheInc : rTheta  + thetaRefineStep*rTheInc;
          azimuthalRefine= rPhi-phiRefineStep*rPhiInc : rPhiInc : rPhi  + phiRefineStep*rPhiInc;

          searchList = zeros(totalRefineStep,3);
          nSearch = 1; 
          for iPhi = azimuthalRefine
            for iTheta = polarRefine 
              for iPsi = inPlaneRefine
                % best iPsi is origin Psi - Phi, no need to subtract here.

                  searchList(nSearch, :) = [iPhi, iTheta, iPsi-iPhi];

                nSearch = nSearch + 1;
              end
            end
          end % end of building angle list
             
          for iRefine = 1:nSearch-1
            for alignLoop = 1:2
              if alignLoop == 1  
                rXYZ = rXYZest;
              elseif alignLoop == 2
                rXYZ = cccStorage3(iRefine,8:10);
              end
              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
              
              RotMat = BH_defineMatrix(searchList(iRefine,:),rotConvention, 'inv');
              RotMat = reshape(angles,3,3) * RotMat;

       


            switch alignLoop
                % This keeps seperate shifts due to windowing and binning from
                % shifts found in CCC
              case 1
               
                % Estimate the peakshift by rotating the ref not the particle.
                iTrimParticle = ...
                       iparticle(padWindow(1,1) + 1:end - padWindow(2,1) , ...
                                 padWindow(1,2) + 1:end - padWindow(2,2) , ...
                                 padWindow(1,3) + 1:end - padWindow(2,3) );
                               
              case 2
            
               % Assuming if class specific symmetry, then some not just 1
                if (flgSymmetry)
                  symmetry = classSymmetry{iGold}(1, classPosition);

                    [ iTrimParticle ] = BH_resample3d(iparticle, RotMat,... 
                                                  rXYZ,...
                                                  {'Bah',symmetry,'linear',1,volBinary_tmp}, ...
                                                  'GPU', 'inv',inputVectors);


%                     iTrimParticle =  iTrimParticle(...
%                                      padWindow(1,1) + 1:end - padWindow(2,1) , ...
%                                      padWindow(1,2) + 1:end - padWindow(2,2) , ...
%                                      padWindow(1,3) + 1:end - padWindow(2,3) );
                    
                         iWedgeMask =  BH_resample3d(iMaxWedgeMask, RotMat, [0,0,0], ...
                                                         {'Bah',symmetry,'linear',1,wdgBinary_tmp},...
                                                         'GPU', 'inv',inputWgtVectors);
                      

                else
                  symmetry = 1;
                  % Transform the particle, and then trim to motif size

                  [ iTrimParticle ] = BH_resample3d(iparticle, RotMat, ...
                                            rXYZ, {'Bah',1,'linear',1,volBinary_tmp}, 'GPU', 'inv',inputVectors);


                         iWedgeMask = BH_resample3d(iMaxWedgeMask, RotMat, [0,0,0], ...
                                           {'Bah',1,'linear',1,wdgBinary_tmp}, 'GPU', 'inv');
                     
                          
%                     iTrimParticle = ...
%                          iTrimParticle(padWindow(1,1) + 1:end - padWindow(2,1) , ...
%                                      padWindow(1,2) + 1:end - padWindow(2,2) , ...
%                                      padWindow(1,3) + 1:end - padWindow(2,3) );


                end % Symmetry or not + interpolation




              
            end



              if alignLoop == 1

                    % use transpose of RotMat
                    try
                    iRotRef = BH_resample3d(ref_FT2_tmp{iGold}{rRef}, RotMat', ...
                                            rXYZ, {'Bah',1,'linear',1,volBinary_tmp}, 'GPU', 'forward',inputVectors);
                    catch
                    cccPreRefineSort(1,1)
                    end                      
                    iRotWdg = BH_resample3d(ref_WGT_rot{iGold}{rRef}, RotMat', ...
                                [0,0,0], {'Bah',1,'linear',1,wdgBinary_tmp}, 'GPU', 'forward',inputWgtVectors);                                          
                                      


                    iRotRef = BH_bandLimitCenterNormalize(...
                                                         iRotRef.*peakMask_tmp,...
                                                         bandpassFiltREF_tmp{rRef},peakBinary_tmp,...
                                                         padCalc,flgPrecision);
                                                       
                    rotPart_FT = BH_bandLimitCenterNormalize(...
                                 iTrimParticle.*peakMask_tmp,...
                                 bandpassFilt_tmp{rRef} ,peakBinary_tmp,padCalc,flgPrecision); 
                               
                    [ peakCoord ] =  BH_multi_xcf_Translational( ...
                                                            rotPart_FT.*ifftshift(iRotWdg), ...
                                                            conj(iRotRef).*iMaxWedgeIfft,...
                                                            peakMask_tmp, peakCOM,eraseMask_tmp);

                                                          
                % 2016-11-11 also took out (+ rXYZ) 
                cccStorage3(iRefine,:) = [rRef, rPart, ...
                                          searchList(iRefine,:), ...
                                          1, 1, ...
                                          peakCoord+rXYZ]; 
              else
                  rotPart_FT = BH_bandLimitCenterNormalize(...
                                            iTrimParticle.*volMask_tmp,...
                                            bandpassFilt_tmp{rRef},volBinary_tmp,...
                                            padCalc,flgPrecision);
               
                      [ iCCC, ~ ] = ...
                                  BH_multi_xcf_Rotational( rotPart_FT, ...
                                                           ref_FT1_tmp{iGold}{rRef},...
                                                           ifftshift(iWedgeMask),...
                                                           ref_WGT_tmp{iGold}{rRef}, ...
                                                           wCCC_tmp{iRef});


                cccStorage3(iRefine,:) = [rRef, rPart, ...
                                          searchList(iRefine,:), ...
                                          iCCC, 1, ...
                                          rXYZ] ;
              end
            end

          end

            sortRef = sortrows(gather(cccStorage3),-6);
            cccStorageRefine{iPeak}(iSubTomo,:) = sortRef(1,:);

        end % end of refinement loop

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        % Get the final translational shift for the best scoring angular
        % match.
        try
          if (flgRefine) && any(cccStorageRefine{iPeak}(iSubTomo,:))
            bestRotPeak = cccStorageRefine{iPeak}(iSubTomo,:);
% % %             % Get the negative slope of the top ten CCC scores.
% % %             topTen = fit([.1:.1:1]',sortRef(1:10,6),'linear');
% % %             bestRotPeak(1,7) = topTen(100)-topTen(101);

          else
            bestRotPeak = cccPreRefineSort(1,:);
            bestRotPeak(1,5) = bestRotPeak(1,5) - bestRotPeak(1,3);
% % %             rowNum = min(size(cccPreRefineSort,1),10*nPeaks);
% % %             topX = 1- 0.1.*(10-rowNum);
% % %             % Get the negative slope of the top ten CCC scores.
% % %             topTen = fit([.1:.1:topX]',cccPreRefineSort(1:rowNum,6),'linear');
% % %             bestRotPeak(1,7) = topTen(100)-topTen(101);

          end
        catch
          fprintf('\nflgRefine %d, iPeak %d, iSubTomo %d\n',flgRefine,iPeak,iSubTomo);
                cccStorageRefine{iPeak}(iSubTomo,:)
          cccPreRefineSort(1,:)
% % %                 rowNum = min(size(cccPreRefineSort,1),10*nPeaks)
% % %                 topX = 1- 0.1.*(10-rowNum)
% % %           fprintf('\nNow check the fits, first and second clause\n');
% % %                 topTen = fit([.1:.1:1]',sortRef(1:10,6),'linear')
% % %           fprintf('\nSecond\n');
% % %           topTen = fit([.1:.1:topX]',cccPreRefineSort(1:rowNum,6),'linear')
% % %           error('Error in sorting the best peak in alignRaw');
        end

          finalRef  = bestRotPeak(1,1);
          finalPart = bestRotPeak(1,2);
          finalPhi  = bestRotPeak(1,3); 
          finalTheta= bestRotPeak(1,4); 
          finalPsi  = bestRotPeak(1,5); 
          % Confirm shiftVAL is doing what it should be
          finalrXYZest  = bestRotPeak(1,8:10);

          RotMat = BH_defineMatrix([finalPhi, finalTheta, finalPsi],rotConvention, 'inv');
          RotMat = reshape(angles,3,3) * RotMat;



    
                  iTrimParticle = ...
                       iparticle(padWindow(1,1) + 1:end - padWindow(2,1) , ...
                                   padWindow(1,2) + 1:end - padWindow(2,2) , ...
                                   padWindow(1,3) + 1:end - padWindow(2,3) );
                                 
 
                    % use transpose of RotMat
                    %%% 2016-11-11 estPeakCoord should have been finalrXYZest in
                    %%% the last writing, but now switching to zeros
                    try
                    iRotRef = BH_resample3d(ref_FT2_tmp{iGold}{finalRef}, RotMat', ...
                                            finalrXYZest, {'Bah',1,'linear',1,volBinary_tmp}, 'GPU', 'forward',inputVectors);
                    iRotWdg = BH_resample3d(ref_WGT_rot{iGold}{finalRef}, RotMat', ...
                                            [0,0,0], {'Bah',1,'linear',1,wdgBinary_tmp}, 'GPU', 'forward',inputWgtVectors);    
                    catch
                      fprintf('\n\nFinal ref,part,phi,theta,psi %f %f %f %f %f\n\n',...
                        bestRotPeak(:,1:5));
                       bestRotPeak(1,1:5)
                       fprintf('BreakPeak %d\n',breakPeak);
                      error('errrorsoedfsdf')
                    end
                                        
%                     iRotRef = ...
%                               iRotRef(padWindow(1,1) + 1:end - padWindow(2,1) , ...
%                                       padWindow(1,2) + 1:end - padWindow(2,2) , ...
%                                       padWindow(1,3) + 1:end - padWindow(2,3) );                                          


                    iRotRef = BH_bandLimitCenterNormalize(...
                                                         iRotRef.*peakMask_tmp,...
                                                         bandpassFiltREF_tmp{finalRef} ,peakBinary_tmp,...
                                                         padCalc,flgPrecision);
                                                       
                    rotPart_FT = BH_bandLimitCenterNormalize(...
                                 iTrimParticle.*peakMask_tmp,...
                                 bandpassFilt_tmp{finalRef} ,peakBinary_tmp,padCalc,flgPrecision ); 
                               
                    [ peakCoord ] =  BH_multi_xcf_Translational( ...
                                                            rotPart_FT.*ifftshift(iRotWdg), ...
                                                            conj(iRotRef).*iMaxWedgeIfft,...
                                                            peakMask_tmp, peakCOM,eraseMask_tmp);

% % %           end




          % Subtract shiftVAL since this is due to windowing, not the actual
          % position.
          cccStorageBest{iPeak}(iSubTomo,:) = gather([bestRotPeak(1,1:7), ...
                                        peakCoord + finalrXYZest - shiftVAL]) ;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
        % It is probably more useful see the shifts in the particle
        % reference frame vs. the avg which was the original
        if (bh_global_print_shifts_in_particle_basis)
          printShifts = zeros(3,3);
          printShifts(1,:) = RotMat * reshape(cccInitial(1,end-2:end),3,1);
          printShifts(2,:) = RotMat * reshape(cccPreRefineSort(1,end-2:end),3,1);
          printShifts(3,:) = RotMat * reshape(cccStorageBest{iPeak}(iSubTomo,end-2:end),3,1);
        else
          printShifts = [cccInitial(1,end-2:end); ...
                         cccPreRefineSort(1,end-2:end);...
                         cccStorageBest{iPeak}(iSubTomo,end-2:end)];
        end
        
        % Print out in Angstrom
        printShifts = printShifts .* pixelSize;
          

        deltaCCC = cccStorageBest{iPeak}(iSubTomo,6) - cccInitial(1,6);
        if (deltaCCC < 0) && (abs(deltaCCC) > 0.15*cccInitial(1,6))
          fprintf('Drop in CCC greater than 15 pph (%2.3f), reverting to prior.\n', deltaCCC);
          fprintf(['\n%s\t%d, %d,%d,%6.3f,%6.3f,%6.3f,%6.6f,%6.6f,%6.3f,%6.3f,%6.3f\n', ...
                   '%s\t%d, %d,%d,%6.3f,%6.3f,%6.3f,%6.6f,%6.6f,%6.3f,%6.3f,%6.3f\n'], ...
                   'PreInitial',iPeak,cccInitial(1,1:end-3),printShifts(1,:),...
                   'PreRefine', iPeak,cccStorageBest{iPeak}(iSubTomo,1:end-3),printShifts(3,:));          
          cccStorageBest{iPeak}(iSubTomo,:) = cccInitial(1,:);
          
        end

        if (flgRefine)
          fprintf(['\n%s\t%d, %d,%d,%d,%6.3f,%6.3f,%6.3f,%6.6f,%6.6f,%6.3f,%6.3f,%6.3f\n', ...
                   '%s\t%d, %d,%d,%d,%6.3f,%6.3f,%6.3f,%6.6f,%6.6f,%6.3f,%6.3f,%6.3f\n', ...
                   '%s\t%d, %d,%d,%d,%6.3f,%6.3f,%6.3f,%6.6f,%6.6f,%6.3f,%6.3f,%6.3f\n'], ...
                   'PreInitial',iPeak,classIDX,cccInitial(1,1:end-3),printShifts(1,:), ...
                   'PreRefine', iPeak,classIDX,[cccPreRefineSort(1,1:4),cccPreRefineSort(1,5)-...
                   cccPreRefineSort(1,3),cccPreRefineSort(1,6:7),printShifts(2,:)], ...
                   'PostRefine',iPeak,classIDX,cccStorageBest{iPeak}(iSubTomo,1:end-3),printShifts(3,:));
          
        else
          fprintf(['\n%s\t%d, %d,%d,%d,%6.3f,%6.3f,%6.3f,%6.6f,%6.6f,%6.3f,%6.3f,%6.3f\n', ...
                   '%s\t%d, %d,%d,%d,%6.3f,%6.3f,%6.3f,%6.6f,%6.6f,%6.3f,%6.3f,%6.3f\n'], ...
                   'PreInitial',iPeak,classIDX,cccInitial(1,1:end-3),printShifts(1,:),...
                   'PreRefine',iPeak,classIDX,cccStorageBest{iPeak}(iSubTomo,1:end-3),printShifts(3,:));
          
        end

    
      end % if condition on newly ignored particles
      
      end

      if ~(rem(iSubTomo,100))
        timeClass = toc;
        fprintf('\nworking on %d/%d subTomo from %s...%fs\n',...
                                      iSubTomo,nSubTomos,tomoName,timeClass);
        tic;
      end


    iParticle = [];
    iSymParti = [];
    iTrimParticle = [];
    iAsym = [];
    iTrimAsym = [];
    iWedgeMask = [];
    rotPart_FT = [];
    rotParticle = [];
        end % end loop over possible peaks
    end % loop over subTomos

    
    for iPeak = 1:nPeaks

      % Get rid of any zero entries left over from pre-initialization
      if iPeak == 1
        nonZeroInits = ( cccStorageBest{iPeak}(:,2) ~= 0 );
        cccStorageBest{1}=cccStorageBest{1}(nonZeroInits,:);
        sortCCC = zeros(size(cccStorageBest{1},1),10*nPeaks);
      else
        cccStorageBest{iPeak}=cccStorageBest{iPeak}(nonZeroInits,:);
      end

      sortCCC(:,1+10*(iPeak-1):10+10*(iPeak-1)) = cccStorageBest{iPeak};
    end

% % %     % I think this is redundant now, but leaving until I double check.
% % %     save('sortCCC.mat','sortCCC');
    [~,a,~] = unique(sortCCC(:,2), 'stable','rows');

    cccSortedandUnique = sortCCC(a,:);
% % %     save('cccSortedandUnique.mat','cccSortedandUnique');
% % %     g = gather(geometry);
% % %     save('TBL_geom.mat','g');

   bestAngles_tmp.(tomoList{iTomo}) = gather(cccSortedandUnique);
   
  % save doesn't work in a parfor, so write out the results for each tomogram so that a 
  % run may be resumed if cancelled.
  angOut = fopen(sprintf('alignResume/%s/%s.txt',outputPrefix,tomoList{iTomo}),'w');
  
  for iRow = 1:size( bestAngles_tmp.(tomoList{iTomo}),1)
    for iPeak = 1:nPeaks
      fprintf(angOut,'%d %d %6.3f %6.3f %6.3f %6.6f %6.6f %6.3f %6.3f %6.3f ', ...
        bestAngles_tmp.(tomoList{iTomo})(iRow,1+10*(iPeak-1):10+10*(iPeak-1)));
    end
    fprintf(angOut,'\n');
  end
  fclose(angOut);
    
  end % if clause to check for previous alignment 
  end % loop over tomos
  bestAnglesResults{iParProc} = bestAngles_tmp;
  geometryResults{iParProc} = geometry_tmp;
%profile off
%profsave
end % parfor



save('bestAnglesResults.mat', 'bestAnglesResults');
bestAngles = struct();
for iParProc = 1:nParProcesses  
  for iTomo = iterList{iParProc}
    geometry.(tomoList{iTomo}) = geometryResults{iParProc}.(tomoList{iTomo});
    bestAngles.(tomoList{iTomo}) = bestAnglesResults{iParProc}.(tomoList{iTomo});
  end
end
%   save('bestAnglesTemp.mat', 'bestAngles');
 save('bestAngles.mat', 'bestAngles');
 
  [ rawAlign ] = BH_rawAlignmentsApply( gather(geometry), bestAngles, samplingRate, nPeaks,rotConvention );
  masterTM.(cycleNumber).('RawAlign') = rawAlign;
  masterTM.(cycleNumber).('newIgnored_rawAlign') = gather(nIgnored);

clear bestAngles rawAlign
subTomoMeta = masterTM;

if ( flgReverseOrder || flgStartThird )
  fprintf('This reverse run will not write the metaData\n');
else
  save(pBH.('subTomoMeta'), 'subTomoMeta');
end

delete(gcp('nocreate'))
for iGPU = 1:nGPUs
  gpuDevice(iGPU);
end

end % end of alignRaw3d
