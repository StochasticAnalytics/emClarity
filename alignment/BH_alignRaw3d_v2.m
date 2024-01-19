 function [  ] = BH_alignRaw3d_v2(PARAMETER_FILE, CYCLE, varargin)
                                                               
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
  bh_global_zero_lag_score = false;
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
CYCLE = EMC_str2double(CYCLE);
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



emc = BH_parseParameterFile(PARAMETER_FILE);
cycleNumber = sprintf('cycle%0.3u', CYCLE);
load(sprintf('%s.mat', emc.('subTomoMeta')), 'subTomoMeta');
mapBackIter = subTomoMeta.currentTomoCPR; 
reconScaling = 1;


try 
  track_stats = emc.('track_stats');
catch
  track_stats = false;
end

try
  flgCutOutVolumes=emc.('flgCutOutVolumes')
catch
  flgCutOutVolumes=0
end


% TODO decide on a "reasonable" padding based on expected shifts.
try
  CUTPADDING = subTomoMeta.('CUTPADDING')
catch
  CUTPADDING=20
end

try
  use_v2_SF3D = emc.('use_v2_SF3D')
catch
  use_v2_SF3D = true;
end

try
  symmetry_op = emc.('symmetry');
catch
  error('You must now specify a symmetry=X parameter, where symmetry E (C1,C2..CX,O,I)');
end

try 
  use_new_grid_search = emc.('use_new_grid_search');
catch
  use_new_grid_search = true;
end

try
  force_no_symmetry = emc.('force_no_symmetry');
catch
  force_no_symmetry = false;
end
if (force_no_symmetry)
  symmetry_op='C1'
  fprintf('\nWarning, overriding symmetry in the alignment. THis is just for benchmarking\n');
end

maxGoldStandard = subTomoMeta.('maxGoldStandard');


nGPUs = emc.('nGPUs')


flgClassify= emc.('flgClassify');
try
  flgMultiRefAlignment=emc.('flgMultiRefAlignment');
catch
  flgMultiRefAlignment = 0;
end

try 
  updateClassByBestReferenceScore = emc.('updateClassByBestReferenceScore');
catch
  updateClassByBestReferenceScore = false;
end
if (~flgMultiRefAlignment)
  updateClassByBestReferenceScore = false;
end

try
  flgCenterRefCOM = emc.('flgCenterRefCOM');
catch
  flgCenterRefCOM = 1;
end

% FIXME: unused, fix experimental options option
try
  flgSymmetrizeSubTomos = emc.('flgSymmetrizeSubTomos');
catch 
  flgSymmetrizeSubTomos = 0;
end
flgRaw_shapeMask =  0;%= emc.('experimentalOpts')(3)
samplingRate = emc.('Ali_samplingRate');

pixelSize      = emc.('PIXEL_SIZE').*10^10.*samplingRate;
if emc.('SuperResolution')
  pixelSize = pixelSize * 2;
end

flgPrecision = 'single'; %emc.('flgPrecision');
angleSearch  = emc.('Raw_angleSearch');
peakSearch   = (emc.('particleRadius')./pixelSize);
peakCOM      = [1,1,1].*3;
className    = emc.('Raw_className');  

try 
  loadTomo = emc.('loadTomo')
catch
  loadTomo = 0;
end
try 
  eraseMaskType = emc.('Peak_mType');
	eraseMaskRadius = emc.('Peak_mRadius')./pixelSize;
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
  doHelical = emc.('doHelical');
catch
  doHelical = 0;
end
if ( doHelical )
  rotConvention = 'Helical'
end

rotConvention

try
  bFactor = emc.('Fsc_bfactor');
catch
  bFactor = 0;
end
if length(bFactor) > 1
  fprintf('multiple bFactors specified, using the first for alignment.\n');
  bFactor = bFactor(1);
end

try
  scaleCalcSize = emc.('scaleCalcSize');
catch
  scaleCalcSize = 1.5;
end
% % % % if (flgClassify || flgMultiRefAlignment)
if (flgClassify)
  refName      = emc.('Ref_className');
else
  refName = emc.('Raw_className');
end
 
outputPrefix = sprintf('%s_%s', cycleNumber, emc.('subTomoMeta'));



classVector{1}  = emc.('Raw_classes_odd')(1,:);


classVector{2}  = emc.('Raw_classes_eve')(1,:);


% % % % if (flgClassify || flgMultiRefAlignment)
if (flgClassify)
  geometry = subTomoMeta.(cycleNumber).ClassAlignment;
  refVectorFull{1}= [emc.('Ref_references_odd');1]
  refVectorFull{2}= [emc.('Ref_references_eve');1]
elseif (flgMultiRefAlignment)
  geometry = subTomoMeta.(cycleNumber).ClusterRefGeom;
  refVectorFull{1}= [emc.('Raw_classes_odd');classVector{1} ]
  refVectorFull{2}= [emc.('Raw_classes_eve');classVector{2} ]
else
  geometry = subTomoMeta.(cycleNumber).Avg_geometry;
  refVectorFull{1} = [emc.('Raw_classes_odd');1];
  refVectorFull{2} = [emc.('Raw_classes_eve');1];
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
nReferences = nReferences .* [~isempty(refGroup{1}),~isempty(refGroup{2})]''


nRefOut(1:2) = [length(unique(refGroup{1})) + sum(( refSym{1} < 0 )),...
                length(unique(refGroup{2})) + sum(( refSym{2} < 0 ))];


%%%%%%%%%%%%%%%%%%%%%%%

% Get the number of tomograms to process.
tomoList = fieldnames(geometry);
nTomograms = length(tomoList);
tiltList = masterTM.tiltGeometry;
ctfGroupList = masterTM.('ctfGroupSize');

% % Sort the list by number of active subtomos to improve parallelism
% sortedTomoList = zeros(nTomograms,1);
% for iTomo = 1:nTomograms
%   sortedTomoList(iTomo) = sum(geometry.(tomoList{iTomo})(:,26)~=-9999);
% end
% [~, sortedTomoIDX] = sort(sortedTomoList,'descend')

% mask defines area for angular search, peakRADIUS restricts translational


[ maskType, maskSize, maskRadius, maskCenter ] = ...
                                  BH_multi_maskCheck(emc, 'Ali', pixelSize)

[ sizeWindow, sizeCalc, sizeMask, padWindow, padCalc ] = ...
                                       BH_multi_validArea( maskSize, maskRadius, scaleCalcSize  )


try 
  flgLimitToOneProcess = emc.('flgLimitToOneProcess');
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
  limitToOne = emc.('nCpuCores');
end

[ nParProcesses, iterList] = BH_multi_parallelJobs(nTomograms,nGPUs, sizeCalc(1),limitToOne);                                   
if ( flgReverseOrder )
  % fprintf('nCpuCores is %d\n', limitToOne);
  % [ nParProcesses, iterList] = BH_multi_parallelJobs(nTomograms,nGPUs, sizeCalc(1),limitToOne);   
  % for iParProc = 1:nParProcesses
  %   iterList{iParProc} = sortedTomoIDX(iterList{iParProc})'
  % end
  % Flip the order for reverse processing on a second machine. This will also disable saving of 
  % of the metadata so there aren't conflicts.
  for iParProc = 1:nParProcesses
    iterList{iParProc} = flip(iterList{iParProc});
  end
  
elseif ( flgStartThird )
  % fprintf('nCpuCores is %d\n', limitToOne);
  % [ nParProcesses, iterList_full] = BH_multi_parallelJobs(nTomograms,nGPUs*cycle_denominator, sizeCalc(1),limitToOne*cycle_denominator);   

  % for iParProc = 1:nParProcesses
  %   iterList_full{iParProc} = sortedTomoIDX(iterList_full{iParProc})';
  % end

  % % Need to scale this back down
  % nParProcesses = limitToOne;

  % Shift to start at one third through to process on a third machine. This will also disable saving of 
  % of the metadata so there aren't conflicts.
  % iterList = {};
  % for iParProc = 1:nParProcesses
  %   idx = cycle_numerator + (iParProc-1)*cycle_denominator;
  %   if (idx <= length(iterList_full))
  %     iterList{iParProc} = iterList_full{idx};
  %   end
  % end
  for iParProc = 1:nParProcesses
    % Note the use of floor is more like ceiling here (rounds away from
    % zero)
    nParts = ceil(length(iterList{iParProc}) ./ cycle_denominator);
    fIDX = 1+(cycle_numerator - 1)*nParts;
    lIDX = min(cycle_numerator*nParts,length(iterList{iParProc}));
    iterList{iParProc} = iterList{iParProc}(fIDX:lIDX);
  end


else
  % error('not supported run config');
  % fprintf('nCpuCores is %d\n', limitToOne);
  % [ nParProcesses, iterList] = BH_multi_parallelJobs(nTomograms,nGPUs, sizeCalc(1),limitToOne);   
  % for iParProc = 1:nParProcesses
  %   iterList{iParProc} = sortedTomoIDX(iterList{iParProc})'
  % end


end

if any(peakSearch > maskRadius)
  fprintf('\n\n\tpeakRADIUS should be <= maskRADIUS!!\n\n')
  peakSearch( (peakSearch > maskRadius) ) = ...
                                        maskRadius( (peakSearch > maskRadius) );
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

 
    wdgPAD = BH_multi_padVal(size(tWDG{1}), sizeCalc);
    for iWdg = 1:n-1
      tWDG_r{iWdg} = BH_padZeros3d(tWDG{iWdg},wdgPAD(1,:),wdgPAD(2,:),...
                                                             'cpu',flgPrecision);
      tWDG{iWdg} = ifftshift(tWDG_r{iWdg});
    end

    refWGT{iGold} = tWDG; clear tWDG wdgTMP
    refWgtROT{iGold} = tWDG_r; clear tWDG_r
 
  
  refIMG{iGold} = tIMG ; clear tIMG refTMP


  clear comMask
  
end

[ refIMG ] = BH_multi_combineLowResInfo( refIMG, imgCounts, pixelSize, maxGoldStandard );



% optimize the fft for the given size. Padding to the next power of 2 is usually
% slower given the dimensionalityl of the volume data.
fftPlanner = rand(sizeCalc);
fftw('planner', 'exhaustive');
fftn(fftPlanner);
clear fftPlanner



  
 

  stat_mask = [];
  if (eraseMask)
    peakMask = EMC_maskShape(eraseMaskType,sizeCalc,floor(eraseMaskRadius),'cpu',{'kernel',false});

    if track_stats
      stat_mask =  single(find(peakMask > 0.95));
    end
      
  else
     if track_stats
      stat_mask =   EMC_maskShape('sphere', sizeCalc, [1,1,1].*floor(max(peakSearch)), 'cpu', {'shift', maskCenter});
      stat_mask = single(find(stat_mask > 0.95));
     end
    [ peakMask ]  = EMC_maskShape('sphere', sizeCalc, [1,1,1].*floor(max(peakSearch)), 'cpu', {'shift', maskCenter;'kernel',false});
  end
  
   
  if ( flgRaw_shapeMask )
 
    [ volMask ] = gather(sqrt(volMask .* ...            
              EMC_maskReference(refIMG{1}{iRef}+refIMG{2}{iRef}, pixelSize, {'fsc', true})));

  else
% % % % % % %     [ volMask ] = gather(BH_mask3d(maskType, sizeWindow, maskRadius, maskCenter));
    [ volMask ] = gather(EMC_maskShape(maskType, sizeWindow, maskRadius, 'gpu', {'shift', maskCenter}));
    
  end     
                                                   


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
                          refTMP_2, '', '', padCalc, 'single')));
    refOUT{nOut+1} = gather(refOUT{nOut+1}(padCalc(1,1) + 1: end - padCalc(2,1), ...
                            padCalc(1,2) + 1: end - padCalc(2,2), ...
                            padCalc(1,3) + 1: end - padCalc(2,3)) );  
    nOut = nOut + 2;
    
    refOUT{nOut} = refOUT{nOut} - mean(refOUT{nOut}(:));
    refOUT{nOut} = refOUT{nOut} ./ rms(refOUT{nOut}(:));
    
    refOUT{nOut+1} = refOUT{nOut+1} - mean(refOUT{nOut+1}(:));
    refOUT{nOut+1} = refOUT{nOut+1} ./ rms(refOUT{nOut+1}(:));
  end
 

  % Save a montage of the masked reference & shape masks if requested.

%   maskedOUTFILE = sprintf('%s_maskedRef-mont_%s.mrc',outputPrefix,halfSet);
%   [ maskedReferences, ~ ] = BH_montage4d(refOUT, '');
%   SAVE_IMG(MRCImage(single(maskedReferences)), maskedOUTFILE);
  

end

clear refIMG refWDG refOUT iRef

%%%%%%%%%%%%%%%%%%%%% Determine the angular search, if any are zero, don't
%%%%%%%%%%%%%%%%%%%%% search at all in that dimension.

updateWeights = false;                                  
gridSearch = '';                                 
if (use_new_grid_search)
  gridSearch = eulerSearch(symmetry_op, angleSearch(1),...
        angleSearch(2),angleSearch(3),angleSearch(4), 0, 0, true);
  nAngles = sum(gridSearch.number_of_angles_at_each_theta);
  inPlaneSearch = gridSearch.parameter_map.psi

try
  symmetry_constrained_search = emc.('symmetry_constrained_search');
  fprintf('Using symmetry constrained search\n');
catch
  symmetry_constrained_search = false;
end
  
  if (symmetry_constrained_search)
    % symmetry expansion on in-plane search only
    if (gridSearch.symmetry_symbol(1) ~= 'C')
      error('symmetry constrained search only implemented for Cn symmetry');
    else
      symmetry_number = EMC_str2double(gridSearch.symmetry_symbol(2:end));
      for iSym = 1:symmetry_number-1
        inPlaneSearch = [inPlaneSearch,inPlaneSearch + iSym.*(360/symmetry_number)];
      end
    end
  end

  flgRefine=false;

  for i = 1:length(gridSearch.parameter_map.phi)
    if gridSearch.parameter_map.phi{i} > 0
      flgRefine=true;
      break;
    end
  end

    
    angleStep = [];
else
  [  nInPlane, inPlaneSearch, angleStep, nAngles] ...
                                      = BH_multi_gridSearchAngles(angleSearch);
  if any(angleStep(:,1))
    flgRefine = true;
  else
    flgRefine = false;
  end       
  
  if sum(angleStep(:,2) > 0)
    updateWeights = true;
  else
  end
end


% [masterTM] = BH_recordAngularSampling( masterTM, cycleNumber, angleStep, inPlaneSearch);                               
                                   
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
if ~(use_v2_SF3D)
  for iParProc = 1:nParProcesses

    % Caclulating weights takes up a lot of memory, so do all that are necessary
    % prior to the main loop -- CHANGE THE CHECK TO JUST READ THE HEADER NOT LOAD
    % THE WEIGHT INTO GPU MEMORY

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
end
parVect = 1:nParProcesses;
fprintf('Starting main loopwith N references %d\n', nReferences(1));
parfor iParProc = parVect
  symmetry = symmetry_op; % Why TF would this be necessary?
% for iParProc = 1:nParProcesses
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
    peakMaskInterpolator  = '';
    peakMaskInterpolator = interpolator(gpuArray(peakMask),[0,0,0],[0,0,0], rotConvention , 'forward', 'C1', false);
    
    if (track_stats)
      mip = struct();
      mip.('mask') = gpuArray(stat_mask); 
    end
    
    
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
      if ~(use_v2_SF3D)
        wgtName = sprintf('cache/%s_bin%d.wgt',iTiltName,samplingRate);       
%         wgtName = sprintf('cache/%s_bin%d.wgt', tomoList{iTomo},...
%                                                 samplingRate);
        maxWedgeMask = BH_unStackMontage4d(1:nCtfGroups,wgtName,...
                                          ceil(sqrt(nCtfGroups)).*[1,1],'');
        maxWedgeIfft = maxWedgeMask;

        for iWdg = 1:length(maxWedgeMask)
          if ~isempty(maxWedgeMask{iWdg})    
            maxWedgeMask{iWdg} = (maxWedgeMask{iWdg} - min(maxWedgeMask{iWdg}(:))) + 1e-3;
            maxWedgeMask{iWdg} = maxWedgeMask{iWdg}.^softenWeight;
            maxWedgeIfft{iWdg} = ifftshift(maxWedgeMask{iWdg});
            
          end
        end
        fprintf('loaded %s.\n',wgtName);

      end
      

        
              % Can't clear inside the parfor, but make sure we don't have two tomograms
      % in memory at once.
     
     tomoNumber = masterTM.mapBackGeometry.tomoName.(tomoList{iTomo}).tomoNumber;
     tiltName = masterTM.mapBackGeometry.tomoName.(tomoList{iTomo}).tiltName;
     reconCoords = masterTM.mapBackGeometry.(tiltName).coords(tomoNumber,:);
     TLT = masterTM.('tiltGeometry').(tomoList{iTomo});
     
     if (flgCutOutVolumes)
       volumeData = [];
     else
      [ volumeData, reconGeometry ] = BH_multi_loadOrBuild( tomoList{iTomo}, ...
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
    cccStorageBest = cell(emc.nPeaks,1);
    cccStorageRefine = cell(emc.nPeaks,1);
    for iPeak = 1:emc.nPeaks
      cccStorageBest{iPeak} = zeros(nSubTomos,10); 
      cccStorageRefine{iPeak}= zeros(nSubTomos,10);
    end
    % reset for each tomogram
    wdgIDX = 0;
    
    for iSubTomo = 1:nSubTomos
      

      make_SF3D = true;
      breakPeak = 0; % for try catch on cut out vols    
      if (wdgIDX ~= positionList(iSubTomo,9)) && ~(use_v2_SF3D)
        % Geometry is sorted on this value so that tranfers are minimized,
        % as these can take up a lot of mem. For 9 ctf Groups on an 80s
        % ribo at 2 Ang/pix at full sampling ~ 2Gb eache.

        wdgIDX = positionList(iSubTomo,9);
        fprintf('pulling the wedge %d onto the GPU\n',wdgIDX);
        % Avoid temporar

        iMaxWedgeMask = []; iMaxWedgeIfft = [];
        iMaxWedgeMask = gpuArray(maxWedgeMask{wdgIDX});
        iMaxWedgeIfft = gpuArray(maxWedgeIfft{wdgIDX});  
        imgWdgInterpolator = '';
        [imgWdgInterpolator, ~] = interpolator(iMaxWedgeMask,[0,0,0],[0,0,0], 'Bah', 'forward', 'C1', false);


      end
      

%         [~,iw1,iw2,iw3] = BH_resample3d(iMaxWedgeMask, eye(3), [0,0,0], ...
%                               {'Bah',1,'linear',1,wdgBinary_tmp}, ...
%                                                          'GPU', 'inv');  
%         inputWgtVectors = {iw1,iw2,iw3};
%         iw1 = []; iw2 = []; iw3 = [];
                                        

      
      for iPeak = 1:emc.nPeaks
        
        if (track_stats)
          measure_noise = true;
          mip.('x') = {};
          mip.('x2') = {};
          mip.('N') = 0;
%           mip.('X') = zeros(1,3,'single','gpuArray');
%           mip.('X2') = zeros(3,3,'single','gpuArray');
        end
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
        half_set = positionList(iSubTomo, 7);
     

        % if classVector{half_set}(1,:) == 0
        %   classPosition = 1;
        %   flgAllClasses = true;
        % else
        %   classPosition = find(classVector{half_set}(1,:) == classIDX);
        %   flgAllClasses = false;
        % end
        % Align all valid subtomos, even if the do not belong to the classes we've selected as references.
        % To ignore particles, remove them with geometry RemoveClases.m
        flgAllClasses = true;



        if (classIDX ~= -9999) && ... % All previously ignored particles
         ( flgAllClasses ||  ismember(classIDX, classVector{half_set}(1,:)) ) 

            
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
         

        if (make_SF3D)
          make_SF3D = false;
          if use_v2_SF3D
            % For now excluding the soften weight.
            [ iMaxWedgeIfft ] = BH_weightMaskMex(sizeCalc, samplingRate, TLT, center,reconGeometry, emc.wiener_constant);       
             imgWdgInterpolator = '';
             % The unshifted mask is kept in texture mem until no longer
             % needed
             [imgWdgInterpolator, ~] = interpolator(iMaxWedgeIfft,[0,0,0],[0,0,0], 'Bah', 'forward', 'C1', false);
             iMaxWedgeIfft =ifftshift(iMaxWedgeIfft);

          end
          % Just use C1 to initialize, whether or not this is the final
          refInterpolator = '';
          refWdgInterpolator= '';
          particleInterpolator= '';

          [refInterpolator, ~] = interpolator(gpuArray(ref_FT2_tmp{1}{1}),[0,0,0],[0,0,0], 'Bah', 'forward', 'C1', false);
          refWdgInterpolator   = interpolator(gpuArray(ref_WGT_rot{half_set}{iRef}),[0,0,0],[0,0,0],'Bah','forward','C1',false);
          particleInterpolator = interpolator(gpuArray(iparticle),[0,0,0],[0,0,0], 'Bah', 'inv', 'C1', false);          
        end
 
        if (use_new_grid_search)
          theta_search = 1:gridSearch.number_of_out_of_plane_angles;
        else
          theta_search = 1:size(angleStep,1);
        end

        for iAngle = theta_search
       
            if (use_new_grid_search)
              theta = gridSearch.parameter_map.theta(iAngle);
              if length(gridSearch.parameter_map.phi{iAngle}) > 1
                phiInc = gridSearch.parameter_map.phi{iAngle}(2)-gridSearch.parameter_map.phi{iAngle}(1);
              else 
                phiInc = 0;
              end
              thetaInc = gridSearch.theta_step;
              numRefIter = gridSearch.number_of_angles_at_each_theta(iAngle);
            else
              theta = angleStep(iAngle,1);  
              phiInc = angleStep(iAngle,3);
              thetaInc = angleStep(iAngle,4);       
              numRefIter = angleStep(iAngle,2)*length(inPlaneSearch)+1;
            
            end
            
            % To prevent only searching the same increments each time in a limited
            % grid search, radomly offset the azimuthal angle by a random number
            % between 0 and 1/2 the azimuthal increment. 
            
            azimuthalRandomizer = (rand(1)-0.5)*phiInc;
           
            
          % Calculate the increment in phi so that the azimuthal sampling is
          % consistent and equal to the out of plane increment.


          
              if (use_new_grid_search)
                % FIXME randomizer passed as bool to eulerSearch
                phi_search = gridSearch.parameter_map.phi{iAngle};
              else
                phi_search = 0:angleStep(iAngle,2);
              end

              
          for iAzimuth = phi_search
            
            if (use_new_grid_search)
              phi = rem(iAzimuth + azimuthalRandomizer,360);
              psiInc = gridSearch.psi_step;
            else
              phi = rem((phiInc * iAzimuth)+azimuthalRandomizer,360);
              psiInc = angleStep(iAngle,5);

            end
            

           for iInPlane = inPlaneSearch
            psi    = iInPlane;
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
                 


                
               
%                   fprintf('Symmetry confirmation %d\n',symmetry);
%                     [ iTrimParticle ] = BH_resample3d(iparticle, RotMat,... 
%                                                   estPeakCoord,...
%                                                   {'Bah',symmetry,'linear',1,volBinary_tmp}, ...
%                                                   'GPU', 'inv',inputVectors);
                    [ iTrimParticle ] = particleInterpolator.interp3d(...
                                                  RotMat,... 
                                                  estPeakCoord,rotConvention ,...
                                                  'inv',symmetry);
                                                  
                    


                     if (getInitialCCC)
%                      [ iTrimInitial ] =  BH_resample3d(iparticle, ...
%                                                     reshape(angles,3,3),... 
%                                                     shiftVAL,...
%                                                     {rotConvention ,symmetry,'linear',1,volBinary_tmp}, ...
%                                                     'GPU', 'inv',inputVectors);
                    [ iTrimInitial ] = particleInterpolator.interp3d(...
                                                  reshape(angles,3,3),... 
                                                  shiftVAL,rotConvention ,...
                                                  'inv',symmetry);
                    
           
% % %                      powerInitial =  sum(abs(iTrimInitial(volBinary_tmp))).^2;   
% % % 
                
%                       iWedgeInitial = BH_resample3d(iMaxWedgeMask, reshape(angles,3,3), [0,0,0], ...
%                                                {rotConvention ,symmetry,'linear',1,wdgBinary_tmp}, ...
%                                                'GPU', 'inv',inputWgtVectors);
                    [ iWedgeInitial ] = imgWdgInterpolator.interp3d(...
                                                   reshape(angles,3,3),... 
                                                   [0,0,0],rotConvention ,...
                                                  'inv',symmetry);                                             
                      
                                 
                     end           

%                         iWedgeMask = BH_resample3d(iMaxWedgeMask, RotMat, [0,0,0], ...
%                                                  {rotConvention ,symmetry,'linear',1,wdgBinary_tmp}, ...
%                                                  'GPU', 'inv',inputWgtVectors);
                    
                    [ iWedgeMask ] = imgWdgInterpolator.interp3d(...
                                                   RotMat,... 
                                                   [0,0,0],rotConvention ,...
                                                  'inv',symmetry); 





                
                 
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

                    [ iRotRef ] = refInterpolator.interp3d(...
                                                   RotMat',... 
                                                  estPeakCoord,rotConvention ,...
                                                  'forward','C1'); 
                     


                    [ iRotWdg ] = refWdgInterpolator.interp3d(...
                                                   RotMat',... 
                                                   [0,0,0],rotConvention ,...
                                                  'forward','C1');  
                                                
                                                
                    [ iRotMask ] = peakMaskInterpolator.interp3d(...
                                                    RotMat',... 
                                                   [0,0,0],rotConvention ,...
                                                  'forward','C1');  
                                                
                                                
                                    

                    % maybe I should be rotating peak mask here in case it has
                    % an odd shape, since we are leaving the proper frame

                    iRotRef = BH_bandLimitCenterNormalize(...
                                                         iRotRef,...
                                                         bandpassFiltREF_tmp{iRef} ,'',...
                                                         padCalc,flgPrecision);

                    rotPart_FT = BH_bandLimitCenterNormalize(...
                                 iTrimParticle,...
                                 bandpassFilt_tmp{iRef} ,'',padCalc,flgPrecision); 
                      
                    if (track_stats && measure_noise)


                        [ ~, mip ] =  BH_multi_xcf_Translational_2( ...
                                                                rotPart_FT, ...
                                                                conj(iRotRef),...
                                                                ifftshift(iRotWdg),...
                                                                iMaxWedgeIfft,...                                                            
                                                                iRotMask, peakCOM,...
                                                                mip);
                                                              

                    end
                       [ peakCoord ] =  BH_multi_xcf_Translational( ...
                                                              rotPart_FT.*ifftshift(iRotWdg), ...
                                                              conj(iRotRef).*iMaxWedgeIfft,...
                                                              iRotMask, peakCOM);                   
                    

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
                                                           ref_FT1_tmp{half_set}{iRef}, ...
                                                           ifftshift(iWedgeInitial),...
                                                           ref_WGT_tmp{half_set}{iRef}, ...
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
                                                           ref_FT1_tmp{half_set}{iRef},...
                                                           ifftshift(iWedgeMask),...
                                                           ref_WGT_tmp{half_set}{iRef}, ...
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
               if (force_no_symmetry)
                 symmetry = 'C1';
               end

%                     [ iTrimParticle ] = BH_resample3d(iparticle, RotMat,... 
%                                                   rXYZ,...
%                                                   {rotConvention ,symmetry,'linear',1,volBinary_tmp}, ...
%                                                   'GPU', 'inv',inputVectors);
                  [ iTrimParticle ] = particleInterpolator.interp3d(...
                                                 RotMat,... 
                                                rXYZ,rotConvention ,...
                                                'inv',symmetry);                                                


%                     iTrimParticle =  iTrimParticle(...
%                                      padWindow(1,1) + 1:end - padWindow(2,1) , ...
%                                      padWindow(1,2) + 1:end - padWindow(2,2) , ...
%                                      padWindow(1,3) + 1:end - padWindow(2,3) );
% 
%                        iWedgeMask =  BH_resample3d(iMaxWedgeMask, RotMat, [0,0,0], ...
%                                                        {rotConvention ,symmetry,'linear',1,wdgBinary_tmp},...
%                                                        'GPU', 'inv',inputWgtVectors);
                                                     
                    [ iWedgeMask ] = imgWdgInterpolator.interp3d(...
                                                   RotMat,... 
                                                   [0,0,0],rotConvention ,...
                                                  'inv',symmetry);                                                      







              
            end



              if alignLoop == 1

                    % use transpose of RotMat
%                     try
%                     iRotRef = BH_resample3d(ref_FT2_tmp{half_set}{rRef}, RotMat', ...
%                                             rXYZ, {rotConvention ,1,'linear',1,volBinary_tmp}, 'GPU', 'forward',inputVectors);
%                     catch
%                     cccPreRefineSort(1,1)
%                     end                      
%                     iRotWdg = BH_resample3d(ref_WGT_rot{half_set}{rRef}, RotMat', ...
%                                 [0,0,0], {rotConvention ,1,'linear',1,wdgBinary_tmp}, 'GPU', 'forward',inputWgtVectors);                                          
%                                       

                    [ iRotRef ] = refInterpolator.interp3d(...
                                                  RotMat',... 
                                                  rXYZ,rotConvention ,...
                                                  'forward','C1');                                          
                    

                    [ iRotWdg ] = refWdgInterpolator.interp3d(...
                                                  RotMat',... 
                                                   [0,0,0],rotConvention ,...
                                                  'forward','C1'); 
                                                
                    [ iRotMask ] = peakMaskInterpolator.interp3d(...
                                                   RotMat',... 
                                                   [0,0,0],rotConvention ,...
                                                  'forward','C1');                                                 

                    iRotRef = BH_bandLimitCenterNormalize(...
                                                         iRotRef,...
                                                         bandpassFiltREF_tmp{rRef},'',...
                                                         padCalc,flgPrecision);
                                                       
                    rotPart_FT = BH_bandLimitCenterNormalize(...
                                 iTrimParticle,...
                                 bandpassFilt_tmp{rRef} ,'',padCalc,flgPrecision); 
                               
                    [ peakCoord ] =  BH_multi_xcf_Translational( ...
                                                            rotPart_FT.*ifftshift(iRotWdg), ...
                                                            conj(iRotRef).*iMaxWedgeIfft,...
                                                            iRotMask, peakCOM);

                                                          
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
                                                           ref_FT1_tmp{half_set}{rRef},...
                                                           ifftshift(iWedgeMask),...
                                                           ref_WGT_tmp{half_set}{rRef}, ...
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


          else
            bestRotPeak = cccPreRefineSort(1,:);
            bestRotPeak(1,5) = bestRotPeak(1,5) - bestRotPeak(1,3);


          end
        catch
          fprintf('\nflgRefine %d, iPeak %d, iSubTomo %d\n',flgRefine,iPeak,iSubTomo);
                cccStorageRefine{iPeak}(iSubTomo,:)
          cccPreRefineSort(1,:)

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
%               iRotRef = BH_resample3d(ref_FT2_tmp{half_set}{finalRef}, RotMat', ...
%                                       finalrXYZest, {rotConvention ,1,'linear',1,volBinary_tmp}, 'GPU', 'forward',inputVectors);
%               iRotWdg = BH_resample3d(ref_WGT_rot{half_set}{finalRef}, RotMat', ...
%                                       [0,0,0], {rotConvention ,1,'linear',1,wdgBinary_tmp}, 'GPU', 'forward',inputWgtVectors);   

                    [ iRotRef ] = refInterpolator.interp3d(...
                                                  RotMat',... 
                                                  finalrXYZest,rotConvention ,...
                                                  'forward','C1');                                          
                    

                    [ iRotWdg ] = refWdgInterpolator.interp3d(...
                                                   RotMat',... 
                                                   [0,0,0],rotConvention ,...
                                                  'forward','C1');            
                                                
                    [ iRotMask ] = peakMaskInterpolator.interp3d(...
                                                   RotMat',... 
                                                   [0,0,0],rotConvention ,...
                                                  'forward','C1');                                                 
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
                                                         iRotRef,...
                                                         bandpassFiltREF_tmp{finalRef} ,'',...
                                                         padCalc,flgPrecision);
                                                       
                    rotPart_FT = BH_bandLimitCenterNormalize(...
                                 iTrimParticle,...
                                 bandpassFilt_tmp{finalRef} ,'',padCalc,flgPrecision ); 
                  

                      [ peakCoord ] =  BH_multi_xcf_Translational( ...
                                                            rotPart_FT.*ifftshift(iRotWdg), ...
                                                            conj(iRotRef).*iMaxWedgeIfft,...
                                                            iRotMask, peakCOM);
                    
                      


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

        if (track_stats)

          if thetaInc > 0
            cccStorageBest{iPeak}(iSubTomo,end-3) = gather(mean(mip.x , 'all')./std(mip.x,0,'all')./thetaInc);
          else
            cccStorageBest{iPeak}(iSubTomo,end-3) = 0;
          end
          
%           % I'm not sold on what do do with this. The distribution over the
%           % shift parameters doesn't really seem to make sense to me. There
%           % are too many factors that can lead to large shifts (e.g.
%           % tomoCPR) If we were searching the full angular space each
%           % iteration, then this would make sense.
%           mip_mean = mip.X./mip.N;
%           mip_covar = mip.X2./mip.N - transpose(mip_mean)*(mip_mean);
%           mip_covar_inv = mip_covar\eye(3);          
%           gauss_norm = ((2.*pi).^(3/2).*abs(mip_covar)).^-1;
%           gauss_exp = exp(-0.5.*(printShifts(2,:)-mip_mean)*mip_covar_inv*transpose(printShifts(2,:)-mip_mean));
       
        end
        
        
        cccInitial(1,1) = classVector{iGold}(cccInitial(1,1));
        cccStorageBest{iPeak}(iSubTomo,1) = classVector{iGold}(cccStorageBest{iPeak}(iSubTomo,1));
        if (flgRefine)
          cccPreRefineSort(1,1) = classVector{iGold}(cccPreRefineSort(1,1));
          fprintf(['\n%s\t%d, %d,%d,%d,%6.3f,%6.3f,%6.3f,%6.6f,%6.6f,%6.3f,%6.3f,%6.3f\n', ...
                   '%s\t%d, %d,%d,%d,%6.3f,%6.3f,%6.3f,%6.6f,%6.6f,%6.3f,%6.3f,%6.3f\n', ...
                   '%s\t%d, %d,%d,%d,%6.3f,%6.3f,%6.3f,%6.6f,%6.6f,%6.3f,%6.3f,%6.3f\n'], ...
                   'PreInitial',iPeak,classIDX, cccInitial(1,1:end-3),printShifts(1,:), ...
                   'PreRefine', iPeak,classIDX,[cccPreRefineSort(1,1:4),cccPreRefineSort(1,5)-...
                   cccPreRefineSort(1,3),cccPreRefineSort(1,6:7),printShifts(2,:)], ...
                   'PostRefine',iPeak,classIDX,cccStorageBest{iPeak}(iSubTomo,1:end-3),printShifts(3,:));
          
        else
          fprintf(['\n%s\t%d, %d,%d,%d,%6.3f,%6.3f,%6.3f,%6.6f,%6.6f,%6.3f,%6.3f,%6.3f\n', ...
                   '%s\t%d, %d,%d,%d,%6.3f,%6.3f,%6.3f,%6.6f,%6.6f,%6.3f,%6.3f,%6.3f\n'], ...
                   'PreInitial',iPeak,classIDX, cccInitial(1,1:end-3),printShifts(1,:),...
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
        
        if use_v2_SF3D
          iMaxWedgeIfft = [];
        end
    end % loop over subTomos

    
    for iPeak = 1:emc.nPeaks

      % Get rid of any zero entries left over from pre-initialization
      if iPeak == 1
        nonZeroInits = ( cccStorageBest{iPeak}(:,2) ~= 0 );
        cccStorageBest{1}=cccStorageBest{1}(nonZeroInits,:);
        sortCCC = zeros(size(cccStorageBest{1},1),10*emc.nPeaks);
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
    for iPeak = 1:emc.nPeaks
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




if ( flgReverseOrder || flgStartThird )
  fprintf('This reverse run will not write the metaData\n');
else
  
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
 
  [ rawAlign ] = BH_rawAlignmentsApply( gather(geometry), bestAngles, samplingRate, emc.nPeaks, rotConvention, updateWeights, updateClassByBestReferenceScore);
  masterTM.(cycleNumber).('RawAlign') = rawAlign;
  masterTM.(cycleNumber).('newIgnored_rawAlign') = gather(nIgnored);
  masterTM.('updatedWeights') = true;

  clear bestAngles rawAlign
  subTomoMeta = masterTM;
  save(emc.('subTomoMeta'), 'subTomoMeta');
  


end

delete(gcp('nocreate'))
for iGPU = 1:nGPUs
  gpuDevice(iGPU);
end

end % end of alignRaw3d
