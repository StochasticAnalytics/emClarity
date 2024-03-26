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

startTime =  datetime("now");
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





maxGoldStandard = subTomoMeta.('maxGoldStandard');


nGPUs = emc.('nGPUs');




flgRaw_shapeMask =  0;%= emc.('experimentalOpts')(3)
samplingRate = emc.('Ali_samplingRate');

emc.pixel_size_angstroms = emc.pixel_size_angstroms.*samplingRate;


flgPrecision = 'single'; %emc.('flgPrecision');
angleSearch  = emc.('Raw_angleSearch');
peakSearch   = (emc.('particleRadius')./emc.pixel_size_angstroms);
peakCOM      = [1,1,1].*3;
className    = emc.('Raw_className');


try
  eraseMaskType = emc.('Peak_mType');
  eraseMaskRadius = emc.('Peak_mRadius')./emc.pixel_size_angstroms;
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

if ( emc.doHelical )
  rotConvention = 'Helical';
end

if (emc.classification)
  refName      = emc.('Ref_className');
else
  refName = emc.('Raw_className');
end

outputPrefix = sprintf('%s_%s', cycleNumber, emc.('subTomoMeta'));



classVector{1}  = emc.('Raw_classes_odd')(1,:);


classVector{2}  = emc.('Raw_classes_eve')(1,:);


% % % % if (emc.classification || emc.multi_reference_alignment)
if (emc.classification)
  geometry = subTomoMeta.(cycleNumber).ClassAlignment;
  refVectorFull{1}= [emc.('Ref_references_odd');1]
  refVectorFull{2}= [emc.('Ref_references_eve');1]
elseif (emc.multi_reference_alignment)
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

refVector = cell(2,1);
refGroup = cell(2,1);
refSym = cell(2,1);

for iGold = 1:2
  % Sort low to high, because order is rearranged as such unstack
  refVectorFull{iGold} = sortrows(refVectorFull{iGold}', 1)';
  % class id corresponding to membership in ???_refName
  refVector{iGold} = refVectorFull{iGold}(1,:);
  % reference id, so multiple classes can be merged into one
  refGroup{iGold}  = refVectorFull{iGold}(3,:);
  % axial symmetry to apply, negative value indicates creating a mirrored ref
  % accros the corresponding axis
  refSym{iGold}    = refVectorFull{iGold}(2,:);
end

% make sure the number of references match the unique groups in the classVector
% and also that the class/group pairs match the class/ref pairs.
nReferences(1:2) = [length(unique(refGroup{1})),length(unique(refGroup{1}))];
nReferences = nReferences .* [~isempty(refGroup{1}),~isempty(refGroup{2})];

if (nReferences(1) ~= length(refVector{1}))
  error('Number of references does not match the number of unique groups in the classVector')
end
nRefOut(1:2) = [length(unique(refGroup{1})) + sum(( refSym{1} < 0 )),...
  length(unique(refGroup{2})) + sum(( refSym{2} < 0 ))];

% Get the number of tomograms to process.
tomoList = fieldnames(geometry);
nTomograms = length(tomoList);
tiltList = subTomoMeta.tiltGeometry;

% Sort the list by number of active subtomos to improve parallelism
sortedTomoList = zeros(nTomograms,1);
for iTomo = 1:nTomograms
  sortedTomoList(iTomo) = sum(geometry.(tomoList{iTomo})(:,26)~=-9999);
end
[~, sortedTomoIDX] = sort(sortedTomoList,'descend');
tomoList = tomoList(sortedTomoIDX);
% mask defines area for angular search, peakRADIUS restricts translational


[ maskType, maskSize, maskRadius, maskCenter ] = ...
  BH_multi_maskCheck(emc, 'Ali', emc.pixel_size_angstroms);

[ sizeWindow, sizeCalc, sizeMask, padWindow, padCalc ] = ...
  BH_multi_validArea( maskSize, maskRadius, emc.scale_calc_size  )



if (flgStartThird)
  [ nParProcesses, iterList] = BH_multi_parallelJobs(nTomograms, nGPUs, sizeCalc(1), emc.nCpuCores, [cycle_numerator,cycle_denominator]);
else
  [ nParProcesses, iterList] = BH_multi_parallelJobs(nTomograms, nGPUs, sizeCalc(1), emc.nCpuCores);
end
if ( flgReverseOrder )
  % Flip the order for reverse processing on a second machine. This will also disable saving of
  % of the metadata so there aren't conflicts.
  for iParProc = 1:nParProcesses
    iterList{iParProc} = flip(iterList{iParProc});
  end
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
  
  
  imgNAME = sprintf('class_%d_Locations_Ref_%s', refName, halfSet)
  
  
  weightNAME = sprintf('class_%d_Locations_Ref_%s_Wgt', refName, halfSet);
  imgCounts{iGold} = subTomoMeta.(cycleNumber).(imgNAME){3};
  
  
  [ refTMP ] = BH_unStackMontage4d(1:nReferences(iGold), ...
    subTomoMeta.(cycleNumber).(imgNAME){1}, ...
    subTomoMeta.(cycleNumber).(imgNAME){2},...
    sizeWindow);
  
  [ wdgTMP ] = BH_unStackMontage4d(1:nReferences(iGold), ...
    subTomoMeta.(cycleNumber).(weightNAME){1},...
    subTomoMeta.(cycleNumber).(weightNAME){2},...
    sizeCalc);
  
  sizeREF = subTomoMeta.(cycleNumber).(imgNAME){2}{1}(2:2:6)';
  
  if (emc.move_reference_by_com)
    % % % % % % %    [ comMask ] = BH_mask3d(maskType, sizeMask, maskRadius, maskCenter);
    [ comMask ]  = EMC_maskShape(maskType, sizeMask, maskRadius, 'gpu', {'shift', maskCenter});
  end
  
  % get boxSize
  n = 1 ; tIMG = cell(numel(refVector{iGold})); tWDG = cell(numel(refVector{iGold}));tWDG_r = tWDG;
  for iP = 1:numel(refTMP)
    if ~isempty(refTMP{iP})
      tIMG{n} = refTMP{iP}; refTMP{iP} = [];
      if (emc.move_reference_by_com)
        
        [~, ~, ~,iCOM] = EMC_maskReference(gpuArray(tIMG{n}).*comMask, emc.pixel_size_angstroms, {'fsc',true; 'com', true});
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

[ refIMG ] = BH_multi_combineLowResInfo( refIMG, imgCounts, emc.pixel_size_angstroms, maxGoldStandard );



% optimize the fft for the given size. Padding to the next power of 2 is usually
% slower given the dimensionalityl of the volume data.
fftPlanner = rand(sizeCalc);
fftw('planner', 'exhaustive');
fftn(fftPlanner);
clear fftPlanner






stat_mask = [];
if (eraseMask)
  peakMask = EMC_maskShape(eraseMaskType,sizeCalc,floor(eraseMaskRadius),'cpu',{'kernel',false});
  
  if ( emc.track_stats )
    stat_mask =  single(find(peakMask > 0.95));
  end
  
else
  if ( emc.track_stats )
    stat_mask =   EMC_maskShape('sphere', sizeCalc, [1,1,1].*floor(max(peakSearch)), 'cpu', {'shift', maskCenter});
    stat_mask = single(find(stat_mask > 0.95));
  end
  [ peakMask ]  = EMC_maskShape('sphere', sizeCalc, [1,1,1].*floor(max(peakSearch)), 'cpu', {'shift', maskCenter;'kernel',false});
end


if ( flgRaw_shapeMask )
  
  [ volMask ] = gather(sqrt(volMask .* ...
    EMC_maskReference(refIMG{1}{iRef}+refIMG{2}{iRef}, emc.pixel_size_angstroms, {'fsc', true})));
  
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
if (emc.classification || emc.multi_reference_alignment)
  for iRef = 1:nReferences(1)

    fscINFO = subTomoMeta.(cycleNumber).('fitFSC').(sprintf('Ref%d',iRef));
  
    [radialGrid,~,~,~,~,~ ] = BH_multi_gridCoordinates(sizeCalc, 'Cartesian', ...
      'GPU', {'none'}, 1, 0, 1 );
    radialGrid = single(radialGrid./emc.pixel_size_angstroms);
    % returns a cpu array
    if (flgWeightCCC)
      [ bandpassFilt{iRef}, ~,wCCC] =  BH_multi_cRef( fscINFO, radialGrid, emc.Fsc_bfactor(1), 1, 1);
    else
      [ bandpassFilt{iRef}, ~] =  BH_multi_cRef( fscINFO, radialGrid, emc.Fsc_bfactor(1), 1);
    end
    
    
    bandpassFiltREF{iRef} = 1;
    
  end
else
  
  
  
  for iRef = 1
    fscINFO = subTomoMeta.(cycleNumber).('fitFSC').('Ref1');
    [radialGrid,~,~,~,~,~ ] = BH_multi_gridCoordinates(sizeCalc, 'Cartesian', ...
      'GPU', {'none'}, 1, 0, 1 );
    radialGrid = single(radialGrid./emc.pixel_size_angstroms);
    % returns a cpu array
    if (flgWeightCCC)
      [ bandpassFilt{iRef},~,wCCC{iRef} ] =  BH_multi_cRef( fscINFO, radialGrid, emc.Fsc_bfactor(1), 1, 1 );
    else
      [ bandpassFilt{iRef},~ ] =  BH_multi_cRef( fscINFO, radialGrid, emc.Fsc_bfactor(1), 1 );
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


    refTMP = refIMG{iGold}{iRef}(padWindow(1,1) + 1: end - padWindow(2,1), ...
                                padWindow(1,2) + 1: end - padWindow(2,2), ...
                                padWindow(1,3) + 1: end - padWindow(2,3));
    
    
    % if not using a weighted average (adapted SPW filter), apply an
    % approximation the cRef from Rosenthal/Henderson. This is currently always set to one
    % and is just doing the masking and normalization. It should be okay to just apply the mask
    % and rely on the normalization during the CCC calc. TODO
    ref_FT1{iGold}{iRef} = gather(conj(BH_bandLimitCenterNormalize(...
      refTMP.*volMask, bandpassFiltREF{iRef}, (volMask>0.01), padCalc, flgPrecision)));
    
    
    
    
    % ref_FT2{iGold}{iRef} = gather(refTMP_2);
    % Trim for output reference
    % refTMP_2 = refTMP_2(padWindow(1,1) + 1: end - padWindow(2,1), ...
    %   padWindow(1,2) + 1: end - padWindow(2,2), ...
    %   padWindow(1,3) + 1: end - padWindow(2,3));
    
    % Overwrite a copy of the filtered, bandpassed ref for output
    % refOUT{nOut} = real(ifftn(conj(ref_FT1{iGold}{iRef})));
    % refOUT{nOut} = gather(refOUT{nOut}(padCalc(1,1) + 1: end - padCalc(2,1), ...
    %   padCalc(1,2) + 1: end - padCalc(2,2), ...
    %   padCalc(1,3) + 1: end - padCalc(2,3)) .* volMask);
    
    
    
    % refOUT{nOut} = refOUT{nOut}.*volMask;
    
    % refOUT{nOut+1} = real(ifftn(BH_bandLimitCenterNormalize(...
    %   refTMP_2, '', '', padCalc, 'single')));
    % refOUT{nOut+1} = gather(refOUT{nOut+1}(padCalc(1,1) + 1: end - padCalc(2,1), ...
    %   padCalc(1,2) + 1: end - padCalc(2,2), ...
    %   padCalc(1,3) + 1: end - padCalc(2,3)) );
    % nOut = nOut + 2;
    
    % refOUT{nOut} = refOUT{nOut} - mean(refOUT{nOut}(:));
    % refOUT{nOut} = refOUT{nOut} ./ rms(refOUT{nOut}(:));
    
    % refOUT{nOut+1} = refOUT{nOut+1} - mean(refOUT{nOut+1}(:));
    % refOUT{nOut+1} = refOUT{nOut+1} ./ rms(refOUT{nOut+1}(:));
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
if (emc.use_new_grid_search)
  gridSearch = eulerSearch(emc.symmetry, angleSearch(1),...
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


% [subTomoMeta] = BH_recordAngularSampling( subTomoMeta, cycleNumber, angleStep, inPlaneSearch);

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


system('mkdir -p alignResume');

system(sprintf('mkdir -p alignResume/%s',outputPrefix));

parVect = 1:nParProcesses;
fprintf('Starting main loopwith N references %d\n', nReferences(1));

% This may be modified in the parfor loop (tho that prob isn't really necessary)
particle_symmetry = emc.symmetry;
if (emc.force_no_symmetry)
  particle_symmetry = 'C1';
end
% for iParProc = parVect

parfor iParProc = parVect
  symmetry = emc.symmetry;

  bestAngles_tmp = struct();
  geometry_tmp = geometry;
  
  gpuIDXList = mod(parVect+nGPUs,nGPUs)+1;
  iGPUidx = gpuIDXList(iParProc);
  gpuDevice(iGPUidx);
  fprintf('parProc %d/%d assigned to GPU %d\n',iParProc,nParProcesses,iGPUidx);
  
  for iTomo = iterList{iParProc}
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
      
      bandpassFilt_tmp = cell(1,1);
      bandpassFiltREF_tmp = cell(1,1);

      
  
      volMask_tmp = gpuArray(volMask);
      volBinary_tmp = single(find( volMask_tmp > 0.01 ));
      peakMaskInterpolator  = '';
      peakMaskInterpolator = interpolator(gpuArray(peakMask),[0,0,0],[0,0,0], rotConvention , 'forward', 'C1', false);
      
      if (emc.track_stats)
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
        end % iWccc
      end % iRef
      
      % sprintf('\nWorking on %d/%d volumes',iTomo,nTomograms)
      tic;
      
      % Load the tomo into gpu
      tomoName = tomoList{iTomo};
      
      tiltGeometry = subTomoMeta.tiltGeometry.(tomoList{iTomo});
      % Load in the geometry for the tomogram, and get number of subTomos.
      positionList = geometry_tmp.(tomoList{iTomo});
      includeList = positionList(:,26:26:26*emc.nPeaks ) ~= -9999;
      includeList = any(includeList,2);
      position_list = positionList(includeList,:);

      binShift = [0,0,0];
      nSubTomos = size(positionList,1);
      
      
      iTiltName = subTomoMeta.mapBackGeometry.tomoName.(tomoName).tiltName;
      
      % Can't clear inside the parfor, but make sure we don't have two tomograms
      % in memory at once.
      
      tomoIdx = subTomoMeta.mapBackGeometry.tomoName.(tomoList{iTomo}).tomoIdx;
      tiltName = subTomoMeta.mapBackGeometry.tomoName.(tomoList{iTomo}).tiltName;
      reconCoords = subTomoMeta.mapBackGeometry.tomoCoords.(tomoList{iTomo});
      
      TLT = subTomoMeta.('tiltGeometry').(tomoList{iTomo});
      
      if (emc.flgCutOutVolumes)
        volumeData = [];
      else
        do_load = false;
        [ volumeData ] = BH_multi_loadOrBuild(tomoList{iTomo}, ...
                                              mapBackIter, ...
                                              samplingRate,...
                                              iGPUidx, ...
                                              do_load);
         volHeader = getHeader(volumeData);
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
      cccStorage2 = cell(emc.nPeaks,1);
      cccInitial_arr = cell(emc.nPeaks,1);
      cccStorageTrans= zeros(1,10, 'single', 'gpuArray');

      for iPeak = 1:emc.nPeaks
        cccStorageBest{iPeak} = zeros(nSubTomos,10, 'single');
        cccStorageRefine{iPeak}= zeros(nSubTomos,10, 'single');
        cccStorage2{iPeak} = zeros(nSubTomos,10, 'single', 'gpuArray');
        cccInitial_arr{iPeak} = zeros(nSubTomos,10, 'single', 'gpuArray');
        cccInitial_arr{iPeak}(:,6) = -9999;
        cccStorageRefine{iPeak}(:,6) = -9999;

      end
      % reset for each tomogram
      wdgIDX = 0;
      for iSubTomo = 1:nSubTomos
        make_SF3D = true;
        breakPeak = 0; % for try catch on cut out vols

 
        for iPeak = 1:emc.nPeaks
          if (emc.track_stats)
            measure_noise = true;
            mip.('x') = {};
            mip.('x2') = {};
            mip.('N') = 0;
          end
          if (breakPeak)
            continue;
          end
          getInitialCCC = 1;
          
          % Check that the given subTomo is not to be ignored
          classIDX = positionList(iSubTomo, 26+26*(iPeak-1));
          particleIDX = positionList(iSubTomo, 4);
          half_set = positionList(iSubTomo, 7);

          if (classIDX == -9999)
            continue;
          end

          center = positionList(iSubTomo,[11:13]+26*(iPeak-1))./samplingRate + binShift;
          angles = positionList(iSubTomo,[17:25]+26*(iPeak-1));
          
          % Find range to extract, and check for domain error.
          if (emc.flgCutOutVolumes)
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
            geometry_tmp.(tomoList{iTomo})(geometry_tmp.(tomoList{iTomo})(:,4) == particleIDX, 26) = -9999;
            continue;
          else
            if (emc.flgCutOutVolumes)
              % Test with some generic padding , only to be used on bin 1 at
              % first!!! TODO add a flag to check this.
              try
                particleOUT_name = sprintf('cache/subtomo_%0.7d_%d.mrc',positionList(iSubTomo,4),iPeak);
                iparticle = gpuArray(OPEN_IMG('single',particleOUT_name,[indVAL(1,1),indVAL(2,1)], ...
                  [indVAL(1,2),indVAL(2,2)], ...
                  [indVAL(1,3),indVAL(2,3)],'keep'));
              catch
                fprintf('\n\nDid not load cut out vol. on subTomo %d FixMEEEEEE\n\n',iSubTomo);
                geometry_tmp.(tomoList{iTomo})(geometry_tmp.(tomoList{iTomo})(:,4) == particleIDX, 26) = -9999;
                breakPeak = 1;
                continue;
              end
            else
  
              iparticle = gpuArray(OPEN_IMG('single', volumeData, [indVAL(1,1),indVAL(2,1)], ...
                                                                  [indVAL(1,2),indVAL(2,2)], ...
                                                                  [indVAL(1,3),indVAL(2,3)],'keep'));
            end
            [ iparticle ] = BH_padZeros3d(iparticle,  padVAL(1,1:3), ...
              padVAL(2,1:3), 'GPU', 'singleTaper');
            
            if (make_SF3D)
              make_SF3D = false;
              % For now excluding the soften weight.
              [ iMaxWedgeIfft ] = BH_weightMaskMex(sizeCalc, samplingRate, TLT, center, reconCoords, emc.wiener_constant);
              imgWdgInterpolator = '';
              % The unshifted mask is kept in texture mem until no longer
              % needed
              [imgWdgInterpolator, ~] = interpolator(iMaxWedgeIfft,[0,0,0],[0,0,0], 'Bah', 'forward', 'C1', false);
              iMaxWedgeIfft =ifftshift(iMaxWedgeIfft);
  
              particleInterpolator= '';
              
              particleInterpolator = interpolator(gpuArray(iparticle),[0,0,0],[0,0,0], 'Bah', 'inv', 'C1', false);
            end
          end % Else clause on windowing (if is a continue)
          
          for iRefIdx = 1:nReferences(1);
            % Just use C1 to initialize, whether or not this is the final
            refInterpolator = '';
            refWdgInterpolator= '';
            ref_FT1_thread_local = gpuArray(ref_FT1{half_set}{iRefIdx});
            refWgtROT_thread_local = gpuArray(refWgtROT{half_set}{iRefIdx});
            refInterpolator = interpolator(ifftn(ref_FT1_thread_local),[0,0,0],[0,0,0], 'Bah', 'forward', 'C1', false);
            refWdgInterpolator = interpolator(refWgtROT_thread_local,[0,0,0],[0,0,0],'Bah','forward','C1',false);
            bandpassFilt_tmp{1} = gpuArray(bandpassFilt{iRefIdx});
            bandpassFiltREF_tmp{1} = gpuArray(bandpassFiltREF{iRefIdx});
            
            
            % In the case we are only aligning refs to those classes they came from, we
            % need to check to see if we duck out here.
            switch emc.multi_reference_alignment
              case 0
                refToAlign = 1;
              case 1
                refToAlign = classVector{half_set}(1,iRefIdx);
              case 2
                refToAlign = find(classIDX ==  classVector{half_set}(1,:));
                if isempty(refToAlign)
                  fprintf('WARNING: wanted class %d not found in classVector\n',classIDX);
                  for i = 1:size(classVector{half_set},2)
                    fprintf('%d ',classVector{half_set}(1,i));
                  end
                  error('classIDX not found in classVector');
                end
                if (refToAlign ~= iRefIdx)
                  continue;
                end
            end
            
              
            if (emc.use_new_grid_search)
              theta_search = 1:gridSearch.number_of_out_of_plane_angles;
            else
              theta_search = 1:size(angleStep,1);
            end
            
            for iAngle = theta_search
              
              if (emc.use_new_grid_search)
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
              
              if (emc.use_new_grid_search)
                % FIXME randomizer passed as bool to eulerSearch
                phi_search = gridSearch.parameter_map.phi{iAngle};
              else
                phi_search = 0:angleStep(iAngle,2);
              end
              
              for iAzimuth = phi_search
                
                if (emc.use_new_grid_search)
                  phi = rem(iAzimuth + azimuthalRandomizer,360);
                  psiInc = gridSearch.psi_step;
                else
                  phi = rem((phiInc * iAzimuth)+azimuthalRandomizer,360);
                  psiInc = angleStep(iAngle,5);
                end
                
                for iInPlane = inPlaneSearch
                  psi    = iInPlane;
                  
                  RotMat = BH_defineMatrix([phi, theta, psi - phi],rotConvention, 'inv');
                  RotMat = reshape(angles,3,3) * RotMat;
                                      
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
                        
                        estPeakCoord = gather(cccStorageTrans(1,8:10));
                      
                        [ iTrimParticle ] = particleInterpolator.interp3d(...
                          RotMat,...
                          estPeakCoord,rotConvention ,...
                          'inv',particle_symmetry);
                        
                        if (getInitialCCC)
                          [ iTrimInitial ] = particleInterpolator.interp3d(...
                            reshape(angles,3,3),...
                            shiftVAL,rotConvention ,...
                            'inv',particle_symmetry);
                          
                          [ iWedgeInitial ] = imgWdgInterpolator.interp3d(...
                            reshape(angles,3,3),...
                            [0,0,0],rotConvention ,...
                            'inv',particle_symmetry);
                        end
                        
                        [ iWedgeMask ] = imgWdgInterpolator.interp3d(...
                          RotMat,...
                          [0,0,0],rotConvention ,...
                          'inv',particle_symmetry);
                    end % switch on align loop
                    
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
                          bandpassFiltREF_tmp{1} ,'',...
                          [0,0,0;0,0,0],flgPrecision);
                        
                        rotPart_FT = BH_bandLimitCenterNormalize(...
                          iTrimParticle,...
                          bandpassFilt_tmp{1} ,'',padCalc,flgPrecision);
                        
                        if (emc.track_stats && measure_noise)
                          
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
                        
                        cccStorageTrans(1,:) = [iRefIdx, particleIDX, ...
                                                phi, theta, psi - phi, ...
                                                0, 0, ...
                                                peakCoord + estPeakCoord];
                      case 2
                        % get starting point
                        if (getInitialCCC)
                          initialRotPart_FT = BH_bandLimitCenterNormalize(...
                            iTrimInitial.*volMask_tmp,...
                            bandpassFilt_tmp{1} ,volBinary_tmp,padCalc,flgPrecision);
                          
                          [ iCCC, ~ ] = ...
                            BH_multi_xcf_Rotational( initialRotPart_FT, ...
                            ref_FT1_thread_local, ...
                            ifftshift(iWedgeInitial),...
                            refWGT{half_set}{iRefIdx}, ...
                            wCCC_tmp{iRefIdx});
                          
                            if ( iCCC > cccInitial_arr{iPeak}(iSubTomo,6))
                              cccInitial_arr{iPeak}(iSubTomo,:) =  [iRefIdx, particleIDX, ...
                                                        0,0,0, ...
                                                        iCCC, 1, ...
                                                        shiftVAL];
                                                      initialRotPart_FT = [];
                            end
                        end
                        
                        rotPart_FT = BH_bandLimitCenterNormalize(...
                          iTrimParticle.*volMask_tmp,...
                          bandpassFilt_tmp{1} ,volBinary_tmp,padCalc,flgPrecision);
                        
                        [ iCCC, ~ ] = ...
                          BH_multi_xcf_Rotational( rotPart_FT, ...
                          ref_FT1_thread_local,...
                          ifftshift(iWedgeMask),...
                          refWGT{half_set}{iRefIdx}, ...
                          wCCC_tmp{iRefIdx});
                        
                        % Note that no new translational estimate is made, so no
                        % need to multiply by RotMat
                        if (iCCC > cccStorage2{iPeak}(iSubTomo,6) )

                          cccStorage2{iPeak}(iSubTomo,:) = ...
                                                      [iRefIdx, particleIDX, ...
                                                        phi, theta, psi , ...
                                                        iCCC, 1, ...
                                                        estPeakCoord];
                        end
                      
                    end % switch
                    % This volume won't be needed until the next subTomo is considered,
                    % which is also where getInitialCCC Boolean is set to True again.
                    iTrimInitial = [];
                  end % alignLoop (trans then rotation)
                  getInitialCCC = 0; % needs to be calculated on the first pass through on the second iter of align loop
                end % in plane angles (psi)
              end % azimuth (phi)
            end % polar (theta)
          end % loop over references

          
          cccInitial = gather( cccInitial_arr{iPeak}(iSubTomo,:) );
          if (cccInitial(1,6) == -9999)
            continue;
          end
          cccPreRefineSort =  gather(cccStorage2{iPeak}(iSubTomo,:));

          % This only seems to be a problem with cut out volumes.
          % Normalization maybe?
          if ~any(cccPreRefineSort(1,:))
            cccStorageBest{iPeak}(iSubTomo,:) = cccInitial(1,:);
            fprintf('all Zeros in PreRefine search, revert on subtomo %d peak %d\n',iSubTomo,iPeak);
            continue;
          end
            %%

          if cccInitial(1,6) > cccPreRefineSort(1,6)
            cccPreRefineSort(1,:) = cccInitial(1,:);
          end
          
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

            % Host to Device for best reference
            bandpassFilt_tmp{1} = gpuArray(bandpassFilt{rRef});
            bandpassFiltREF_tmp{1} = gpuArray(bandpassFiltREF{rRef});
            refInterpolator = '';
            refWdgInterpolator= '';
            ref_FT1_thread_local = gpuArray(ref_FT1{half_set}{rRef});
            refWgtROT_thread_local = gpuArray(refWgtROT{half_set}{rRef});
            refInterpolator = interpolator(ifftn(ref_FT1_thread_local),[0,0,0],[0,0,0], 'Bah', 'forward', 'C1', false);
            refWdgInterpolator = interpolator(gpuArray(refWgtROT_thread_local),[0,0,0],[0,0,0],'Bah','forward','C1',false);
                
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
               
                RotMat = BH_defineMatrix(searchList(iRefine,:),rotConvention, 'inv');
                RotMat = reshape(angles,3,3) * RotMat;
                  
                switch alignLoop
                  % This keeps seperate shifts due to windowing and binning from
                  % shifts found in CCC
                  case 1
                    rXYZ = rXYZest;
                    % Estimate the peakshift by rotating the ref not the particle.
                    iTrimParticle = ...
                      iparticle(padWindow(1,1) + 1:end - padWindow(2,1) , ...
                      padWindow(1,2) + 1:end - padWindow(2,2) , ...
                      padWindow(1,3) + 1:end - padWindow(2,3) );

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
                        bandpassFiltREF_tmp{1},'',...
                          [0,0,0;0,0,0],flgPrecision);
                      
                      rotPart_FT = BH_bandLimitCenterNormalize(...
                        iTrimParticle,...
                        bandpassFilt_tmp{1} ,'',padCalc,flgPrecision);
      
                        [ peakCoord ] =  BH_multi_xcf_Translational( ...
                          rotPart_FT.*ifftshift(iRotWdg), ...
                          conj(iRotRef).*iMaxWedgeIfft,...
                          iRotMask, peakCOM);
                          
                      % 2016-11-11 also took out (+ rXYZ)
                      cccStorage3(iRefine,:) = [rRef, rPart, ...
                        searchList(iRefine,:), ...
                        -9999, 1, ...
                        peakCoord+rXYZ];     
                  case 2
                    rXYZ = cccStorage3(iRefine,8:10);
                  
                    [ iTrimParticle ] = particleInterpolator.interp3d(...
                      RotMat,...
                      rXYZ,rotConvention ,...
                      'inv',particle_symmetry);
                    
                    [ iWedgeMask ] = imgWdgInterpolator.interp3d(...
                      RotMat,...
                      [0,0,0],rotConvention ,...
                      'inv',particle_symmetry);

                      rotPart_FT = BH_bandLimitCenterNormalize(...
                        iTrimParticle.*volMask_tmp,...
                        bandpassFilt_tmp{1},volBinary_tmp,...
                        padCalc,flgPrecision);
                      
                      [ iCCC, ~ ] = ...
                        BH_multi_xcf_Rotational( rotPart_FT, ...
                        ref_FT1_thread_local,...
                        ifftshift(iWedgeMask),...
                        refWGT{half_set}{rRef}, ...
                        wCCC_tmp{rRef});
                      
                      
                      cccStorage3(iRefine,:) = [rRef, rPart, ...
                        searchList(iRefine,:), ...
                        iCCC, 1, ...
                        rXYZ] ;
                end % switch align loop
         
              end % end of alignLoop
            end % end of iRefine loop
            sortRef = sortrows(gather(cccStorage3),-6);
            
            if (sortRef(1,6) > cccStorageRefine{iPeak}(iSubTomo,6))
              cccStorageRefine{iPeak}(iSubTomo,:) = sortRef(1,:);
            end

          end % end of if flgRegine
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

          if (~flgRefine)
            % We can use the same reference and interpolators if we've reloaded them for the refine step, otherwise
            % we need to transfer them host to device here to get the final shifts.
            bandpassFilt_tmp{1} = gpuArray(bandpassFilt{finalRef});
            bandpassFiltREF_tmp{1} = gpuArray(bandpassFiltREF{finalRef});
            refInterpolator = '';
            refWdgInterpolator= '';
            ref_FT1_thread_local = gpuArray(ref_FT1{half_set}{finalRef});
            refWgtROT_thread_local = gpuArray(refWgtROT{half_set}{finalRef});
            refInterpolator = interpolator(ifftn(ref_FT1_thread_local),[0,0,0],[0,0,0], 'Bah', 'forward', 'C1', false);
            refWdgInterpolator = interpolator(gpuArray(refWgtROT_thread_local),[0,0,0],[0,0,0],'Bah','forward','C1',false);
          end
            
          iTrimParticle = ...
            iparticle(padWindow(1,1) + 1:end - padWindow(2,1) , ...
            padWindow(1,2) + 1:end - padWindow(2,2) , ...
            padWindow(1,3) + 1:end - padWindow(2,3) );
            
          % use transpose of RotMat
          %%% 2016-11-11 estPeakCoord should have been finalrXYZest in
          %%% the last writing, but now switching to zeros
          try
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
            
          iRotRef = BH_bandLimitCenterNormalize(...
            iRotRef,...
            bandpassFiltREF_tmp{1} ,'',...
              [0,0,0;0,0,0],flgPrecision);
            
          rotPart_FT = BH_bandLimitCenterNormalize(...
          iTrimParticle,...
          bandpassFilt_tmp{1} ,'',padCalc,flgPrecision );
          
          
          [ peakCoord ] =  BH_multi_xcf_Translational( ...
                                                    rotPart_FT.*ifftshift(iRotWdg), ...
                                                    conj(iRotRef).*iMaxWedgeIfft,...
                                                    iRotMask, peakCOM);
          
          
          % Subtract shiftVAL since this is due to windowing, not the actual
          % position.
          cccStorageBest{iPeak}(iSubTomo,:) = gather([bestRotPeak(1,1:7), peakCoord + finalrXYZest - shiftVAL]) ;
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          
          % It is probably more useful see the shifts in the particle
          % reference frame vs. the avg which was the original
          if (emc.printShiftsInParticleBasis)
            printShifts = zeros(3,3);
            printShifts(1,:) = RotMat * reshape(cccInitial(1,end-2:end),3,1);
            printShifts(2,:) = RotMat * reshape(cccPreRefineSort(1,end-2:end),3,1);
            printShifts(3,:) = RotMat * reshape(cccStorageBest{iPeak}(iSubTomo,end-2:end),3,1);
          else
            printShifts =  [cccInitial(1,end-2:end); ...
                            cccPreRefineSort(1,end-2:end);...
                            cccStorageBest{iPeak}(iSubTomo,end-2:end)];
          end
              
          % Print out in Angstrom
          printShifts = printShifts .* emc.pixel_size_angstroms;
          
          
          deltaCCC = cccStorageBest{iPeak}(iSubTomo,6) - cccInitial(1,6);
          if (emc.print_alignment_stats && deltaCCC < 0 && abs(deltaCCC) > 0.15*cccInitial(1,6))
            fprintf('Drop in CCC greater than 15 pph (%2.3f), reverting to prior.\n', deltaCCC);
            fprintf(['\n%s\t%d, %d,%d,%6.3f,%6.3f,%6.3f,%6.6f,%6.6f,%6.3f,%6.3f,%6.3f\n', ...
              '%s\t%d, %d,%d,%6.3f,%6.3f,%6.3f,%6.6f,%6.6f,%6.3f,%6.3f,%6.3f\n'], ...
              'PreInitial',iPeak,cccInitial(1,1:end-3),printShifts(1,:),...
              'PreRefine', iPeak,cccStorageBest{iPeak}(iSubTomo,1:end-3),printShifts(3,:));
            cccStorageBest{iPeak}(iSubTomo,:) = cccInitial(1,:);
            
          end
          
          if (emc.track_stats)
            
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
          if (emc.print_alignment_stats && flgRefine)
            cccPreRefineSort(1,1) = classVector{iGold}(cccPreRefineSort(1,1));
            fprintf(['\n%s\t%d, %d,%d,%d,%6.3f,%6.3f,%6.3f,%6.6f,%6.6f,%6.3f,%6.3f,%6.3f\n', ...
              '%s\t%d, %d,%d,%d,%6.3f,%6.3f,%6.3f,%6.6f,%6.6f,%6.3f,%6.3f,%6.3f\n', ...
              '%s\t%d, %d,%d,%d,%6.3f,%6.3f,%6.3f,%6.6f,%6.6f,%6.3f,%6.3f,%6.3f\n'], ...
              'PreInitial',iPeak,classIDX, cccInitial(1,1:end-3),printShifts(1,:), ...
              'PreRefine', iPeak,classIDX,[cccPreRefineSort(1,1:4),cccPreRefineSort(1,5)-...
              cccPreRefineSort(1,3),cccPreRefineSort(1,6:7),printShifts(2,:)], ...
              'PostRefine',iPeak,classIDX,cccStorageBest{iPeak}(iSubTomo,1:end-3),printShifts(3,:));
            
          else 
            if (emc.print_alignment_stats)
              fprintf(['\n%s\t%d, %d,%d,%d,%6.3f,%6.3f,%6.3f,%6.6f,%6.6f,%6.3f,%6.3f,%6.3f\n', ...
                '%s\t%d, %d,%d,%d,%6.3f,%6.3f,%6.3f,%6.6f,%6.6f,%6.3f,%6.3f,%6.3f\n'], ...
                'PreInitial',iPeak,classIDX, cccInitial(1,1:end-3),printShifts(1,:),...
                'PreRefine',iPeak,classIDX,cccStorageBest{iPeak}(iSubTomo,1:end-3),printShifts(3,:));
            end
          end

          
          if ~(rem(iSubTomo,100))
            timeClass = toc;
            fprintf('Refining %d/%d subTomo from %s...%fs\n', iSubTomo, nSubTomos, tomoName, timeClass);
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
      iMaxWedgeIfft = [];
      
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
    end % if/ -> else clause to check for previous alignment, return to loop on references
  end % loop over tomos
  bestAnglesResults{iParProc} = bestAngles_tmp;
  geometryResults{iParProc} = geometry_tmp;
  %profile off
  %profsave
end % parfor




if ( flgReverseOrder || flgStartThird )
  fprintf('This multi-node run will not write the metaData\n');
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
  
  [ rawAlign ] = BH_rawAlignmentsApply( gather(geometry), bestAngles, samplingRate, emc.nPeaks, rotConvention, updateWeights, emc.update_class_by_ccc);
  subTomoMeta.(cycleNumber).('RawAlign') = rawAlign;
  subTomoMeta.(cycleNumber).('newIgnored_rawAlign') = gather(nIgnored);
  subTomoMeta.('updatedWeights') = true;
  
  clear bestAngles rawAlign
  save(emc.('subTomoMeta'), 'subTomoMeta');
   
end

delete(gcp('nocreate'))
for iGPU = 1:nGPUs
  gpuDevice(iGPU);
end

end % end of alignRaw3d
