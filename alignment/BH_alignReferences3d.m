efunction [  ] = BH_alignReferences3d( PARAMETER_FILE, CYCLE)
%Align an symmetrize references
%   


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (nargin ~= 2)
  error('args = PARAMETER_FILE, CYCLE')
end
startTime =  clock;
CYCLE = str2num(CYCLE);

% not for normal use, but under certain circumstances allow override of 0.5 cutoff for 
% ref and allow to 0.143
if (CYCLE < 0)
  CYCLE = abs(CYCLE);
  flgAlignCutoff = 0;
else
  flgAlignCutoff = 1;
end

cycleNumber = sprintf('cycle%0.3u', CYCLE);

pBH = BH_parseParameterFile(PARAMETER_FILE);
load(sprintf('%s.mat', pBH.('subTomoMeta')), 'subTomoMeta');

maxGoldStandard = subTomoMeta.('maxGoldStandard');
flgClassify = pBH.('flgClassify');
try
  flgMultiRefAlignment = pBH.('flgMultiRefAlignment');
catch
  flgMultiRefAlignment = 0;
end


memoryLevel = pBH.('avgMemory');

% specify an angular shift to apply to the first reference, "" translational,
% identify the ref to align all other refs to, and whether to use bandpass from
% param (based on est from FSC on RawAlign average, or re-extract and calc FSC
% to design filter. The second option will only work in later iterations as the
% alignment of the raw images gets better, as symmetry will magnify small
% geometrical errors.)
bFactor = pBH.('Fsc_bfactor');
samplingRate   = pBH.('Ali_samplingRate');
pixelSize      = pBH.('PIXEL_SIZE').*10^10.*samplingRate;
if pBH.('SuperResolution')
  pixelSize = pixelSize * 2;
end
try
  scaleCalcSize = pBH.('scaleCalcSize');
catch
  scaleCalcSize = 1.5;
end
angleSearch  = pBH.('Ref_angleSearch');

refName    = pBH.('Ref_className');

peakSearch   = floor(pBH.('particleRadius')./pixelSize);
peakCOM      = [1,1,1].*peakCOM;
outputPrefix = sprintf('%s_%s', cycleNumber, pBH.('subTomoMeta'));



flgAngleShift{1}= pBH.('ref_AngleShift_odd');
flgTransShift{1}= pBH.('ref_TransShift_odd');
flgRefRef{1}    = pBH.('ref_Ref_odd');
refVectorFull{1}= pBH.('Ref_references_odd');
features{1}     = pBH.('Pca_coeffs_odd');
% % % vol_geometry{1} = subTomoMeta.(cycleNumber).('ClusterResults').( ...
% % %                                          sprintf('%s_%d_%d_nClass_%d_ODD', ...
% % %                                          outputPrefix,features{1}(1,1), ...
% % %                                          features{1}(1,end), refName));

flgAngleShift{2}= pBH.('ref_AngleShift_eve');
flgTransShift{2}= pBH.('ref_TransShift_eve');
flgRefRef{2}    = pBH.('ref_Ref_eve');
refVectorFull{2}= pBH.('Ref_references_eve');
features{2}     = pBH.('Pca_coeffs_eve');
% % % vol_geometry{2} = subTomoMeta.(cycleNumber).('ClusterResults').( ...
% % %                                          sprintf('%s_%d_%d_nClass_%d_EVE', ...
% % %                                          outputPrefix,features{2}(1,1), ...
% % %                                          features{2}(1,end), refName));
% Merge alignments that could come from using different feature vectors in the
% classification, back into one metadata.
% % % vol_geometry = BH_mergeClassGeometry(vol_geometry{1},vol_geometry{2});
  % Copy to reset after initial extraction of primary references.
vol_geometry = subTomoMeta.(cycleNumber).('ClusterRefGeom'); 

vol_geometry_clean = vol_geometry;









masterTM = subTomoMeta; clear subTomoMeta



refVector = cell(2,1);
refGroup = cell(2,1);
refSym = cell(2,1);

for iGold = 1:2
  refVectorFull{iGold}
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
nReferences = nReferences .* [~isempty(refGroup{1}),~isempty(refGroup{2})];

uniqueSym = cell(2,1);
for iGold = 1:2
  [~,uniqueGroup,~] = unique(refGroup{iGold});
  uniqueSym{iGold} = refSym{iGold}(uniqueGroup)
end

    



[ maskType, maskSize, maskRadius, maskCenter ] = ...
                                  BH_multi_maskCheck(pBH, 'Cls', samplingRate);
                                
[ sizeWindow, sizeCalc, sizeMask, padWindow, padCalc ] = ...
                                       BH_multi_validArea(maskSize, maskRadius, scaleCalcSize)
padREF = [0,0,0;0,0,0];                                     


if any(peakSearch > maskRadius)
  fprintf('\n\n\tpeakRADIUS should be <= maskRADIUS!!\n\n')

  peakSearch( (peakSearch > maskRadius) ) = ...
                                        maskRadius( (peakSearch > maskRadius) );
end



% To properly average averages without re-extracting, the original number
% in each average must be taken in to account.


% Read in the references.
refIMG = cell(2,1);
imgCounts = cell(2,1);
for iGold = 1:2

  if iGold == 1
    halfSet = 'ODD';
  else
    halfSet = 'EVE';
  end

  
  % To allow more flexibility ref/ cls field prefixes are used, and rather than
  % making something special for the naming of cls/ref locations which SHOULD be
  % interchangeable, just try one then the other since this impacts nothing
  % else. Use the unweighted average for the first pass

    try
      imgNAME = sprintf('class_%d_Locations_Ref_%s', refName, halfSet);
      imgCounts{iGold} = masterTM.(cycleNumber).(imgNAME){3};

      [ refIMG{iGold} ] = BH_unStackMontage4d(refVector{iGold}, ...
                                       masterTM.(cycleNumber).(imgNAME){1}, ...
                                       masterTM.(cycleNumber).(imgNAME){2},...
                                       sizeWindow);
    catch
      imgNAME = sprintf('class_%d_Locations_Cls_%s', refName, halfSet);
      imgCounts{iGold} = masterTM.(cycleNumber).(imgNAME){3};

      [ refIMG{iGold} ] = BH_unStackMontage4d(refVector{iGold}, ...
                                       masterTM.(cycleNumber).(imgNAME){1}, ...
                                       masterTM.(cycleNumber).(imgNAME){2},...
                                       sizeWindow);
    end

  occCell = BH_multi_isCell( refIMG{iGold} );
  nRefs = length(occCell);
  tIMG = cell(nRefs,1);
  for iRef = 1:nRefs
    tIMG{iRef} = refIMG{iGold}{occCell(iRef)};
  end
  refIMG{iGold} = tIMG; clear tIMG
  
  sizeREF = masterTM.(cycleNumber).(imgNAME){2}{1};
  sizeREF = sizeREF(2:2:6)'                                      
                            
% % %   % clear out empty cell contents, and return  
% % %   n = 1 ; tIMG = cell(nReferences(iGold)); 
% % %   for iP = 1:length(refTMP)
% % %     if ~isempty(refTMP{iP})
% % %       tIMG{n} = refTMP{iP};
% % % 
% % %       n = n + 1;
% % %     end
% % %   end
% % % 
% % %   refIMG{iGold} = tIMG ; clear tIMG refTMP
% % %   

  
end

% % Prevent divergent orientations
% [ refIMG ] = BH_multi_combineLowResInfo( refIMG, imgCounts, pixelSize, 40 );

% % % [ sizeWindow, sizeCalc, sizeMask, padWindow, padCalc, padREF ] = ...
% % %                                        BH_multi_validArea( maskRadius, sizeREF )
                                     

% optimize the fft for the given size. Padding to the next power of 2 is usually
% slower given the dimensionality of the volume data.
fftPlanner = rand(sizeCalc);
fftw('planner', 'exhaustive');
fftn(fftPlanner);
clear fftPlanner


% Make a mask, and apply to the average motif && save a masked, binned copy of
% the average for inspection. 
mask = struct();


  [ volMask ] = BH_mask3d(maskType, sizeMask, maskRadius, maskCenter);
  mask.('volMask')  = gather(volMask);
                                                     
  [ peakMask ] = gpuArray(BH_mask3d(maskType, sizeMask, peakSearch, maskCenter)); 
  mask.('peakMask') = gather(peakMask);
  
  peakBinary = (peakMask >= 0.01);
 
% The bandpass here is from the full average - assuming each class is a
% substantial portion of the total population.
  
  bandpassFilt = cell(1,1);
  [radialGrid,~,~,~,~,~ ] = BH_multi_gridCoordinates(sizeCalc, 'Cartesian', ...
                                                       'cpu', {'none'}, 1, 0, 1 );
  radialGrid = single(radialGrid./pixelSize);
  for iRef = 1:1
    fscINFO = masterTM.(cycleNumber).('fitFSC').('Raw1');

    % The class averages have roughly the same SNR as the references  so apply any
    % bFactor to them as well.    
    [ bandpassFilt{iRef}, ~ ] =  BH_multi_cRef( fscINFO, radialGrid , bFactor,flgAlignCutoff); 
  end
  clear radialGrid
  mask.('bandpassFilt') = (bandpassFilt); 
  bandpassFilt{1} = gpuArray(bandpassFilt{1});
%   [maskMontage,~] = BH_montage4d({gather(volMask), gather(peakMask)},'');
%   SAVE_IMG(MRCImage(maskMontage), sprintf('%s_maskMontage.mrc', outputPrefix));

%Improve the estimated center for the reference used to align any other
% optional references or sub references, applying any symmetry as well. Then
% re-extract these refs with the new transformation.
estShifts = cell(2,1);

for iGold = 1:2

  if iGold == 1
    halfSet = 'ODD';
  else
   halfSet = 'EVE';
  end
  halfNUM = iGold;

  % The first class in the refRef vector is used to align the primary reference
  % for each reference group, to which the members of each group are
  % subsequently aligned.
  refRefIDX = find(refVector{iGold} == flgRefRef{iGold}(1));
  refRef = refIMG{iGold}{refRefIDX};
  % Only apply trans shifts from estimate, then later use angles as well.
  refRef = BH_axialSymmetry(refRef,1, 0,'GPU',flgTransShift{iGold}(1,:));

  refRef = refRef(padWindow(1,1)+1 : end - padWindow(2,1), ...
                  padWindow(1,2)+1 : end - padWindow(2,2), ...
                  padWindow(1,3)+1 : end - padWindow(2,3) );

  refAsym= refIMG{iGold}{refRefIDX}(padWindow(1,1)+1 : end - padWindow(2,1), ...
                                    padWindow(1,2)+1 : end - padWindow(2,2), ...
                                    padWindow(1,3)+1 : end - padWindow(2,3) );

                                  
  
  refAsym= gpuArray(refAsym);

  
  
  
  refAsym_FT = BH_bandLimitCenterNormalize(refAsym.*peakMask, bandpassFilt{1}, ...
                                                   peakBinary,padCalc,'double');
  refRef_FT =conj(BH_bandLimitCenterNormalize(refRef.*peakMask, bandpassFilt{1},...
                                                  peakBinary,padCalc,'double'));
                                          

  [ estPeakCoord ] =  BH_multi_xcf_Translational( refAsym_FT, refRef_FT, ...
                                                              peakMask, peakCOM);
                                          
  estShifts{iGold} =  gather(estPeakCoord); 
  fprintf('estPeak at %f %f %f\n estShifts now %f %f %f \n', estPeakCoord', estShifts{iGold}');
  clear refRef refAsym refRef_FT refAsym_FT


  refRef = BH_axialSymmetry(refIMG{iGold}{refRefIDX}, refSym{iGold}(refRefIDX),...
                           flgAngleShift{iGold}(1), ...
                           'GPU',estShifts{iGold});
  refRefRotAvg = BH_axialSymmetry(refIMG{iGold}{refRefIDX}, 120,...
                           flgAngleShift{iGold}(1), ...
                           'GPU',estShifts{iGold});

  if isa(refRef,'cell')
    refRef = refRef{2};
  end

  SAVE_IMG(MRCImage(gather(refRef)),sprintf('initialAxialOffsetCheck_%s.mrc',halfSet));
  refWDG = NaN;
 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % % 
  % To allow more flexibility ref/ cls field prefixes are used, and rather than
  % making something special for the naming of cls/ref locations which SHOULD be
  % interchangeable, just try one then the other since this impacts nothing
  % else.

    try
      imgClassNAME = sprintf('class_%d_Locations_%s_%s', refName, 'Ref',halfSet);


      [ classIMG ] = BH_unStackMontage4d(1:refName, ...
                                        masterTM.(cycleNumber).(imgClassNAME){1},...
                                        masterTM.(cycleNumber).(imgClassNAME){2},...
                                        sizeWindow);
    catch
      imgClassNAME = sprintf('class_%d_Locations_%s_%s', refName, 'Cls',halfSet);


      [ classIMG ] = BH_unStackMontage4d(1:refName, ...
                                        masterTM.(cycleNumber).(imgClassNAME){1},...
                                        masterTM.(cycleNumber).(imgClassNAME){2},...
                                        sizeWindow);
    end




  %%%%%%%%%%%%%%%%%%%%% Determine the angular search, if any are zero, don't
  %%%%%%%%%%%%%%%%%%%%% search at all in that dimension.

  [  nInPlane, inPlaneSearch, angleStep, nAngles] ...
                                        = BH_multi_gridSearchAngles(angleSearch(1,:))



  % First align the reference reference to any sub-references in its group. After
  % finding the best alignment extract from the tomograms, with any symmetry.

  % Second align the first member of each ref group to the new reference, and use
  % alignment to resample the head of each ref group.

  % Third align each sub-ref to the head of each ref group, and finally extract
  % these.



%   % Store the cross correlation score, peak location, and wedge weight
  bestAnglesTotal = zeros(nReferences(iGold),10);
  nCount = 1;

  for iClass = flgRefRef{iGold}(1,:)
    tic;
   % Load the class onto gpu
    iClassImg = gpuArray(classIMG{iClass});

    iClassWdg = NaN;

   

    [ alignOUT ] = run_alignment(angleStep, inPlaneSearch, ...
                                 iClassImg, iClassWdg, ...                                          
                                 refRef, refWDG, refRefRotAvg, ...
                                 volMask, bandpassFilt, padCalc, padWindow, ...
                                 peakMask, peakCOM, iClass, ...
                                 uniqueSym{iGold});


    bestAnglesTotal(nCount,:) = [gather(alignOUT(1,1:10))]; 

    timeClass = toc;
    fprintf('finished working on %d/%d classes...%fs\n',nCount, ...
                                              length(flgRefRef{iGold}),timeClass);
    nCount = nCount + 1;

  end % loop over (classes)

  bestAnglesTotal = gather(bestAnglesTotal);
  

  [ vol_geometry ] = BH_refAlignmentsApply( vol_geometry, bestAnglesTotal,...
                                        samplingRate,  ...
                                        flgRefRef{iGold}(1,:), flgRefRef{iGold}(2,:),halfNUM );



  
end 

masterTM.(cycleNumber).('RefAlignment') = vol_geometry;
subTomoMeta = masterTM;
save(pBH.('subTomoMeta'), 'subTomoMeta');
  
 
gpuDevice(1);
BH_average3d(PARAMETER_FILE,num2str(CYCLE),'RefAlignment');

% average3d puts info about the ref locations in the montage, so reload the
% metadata

load(pBH.('subTomoMeta'));
masterTM = subTomoMeta;

% reload masks onto gpu
volMask = gpuArray(mask.('volMask'));
peakMask = gpuArray(mask.('peakMask'));
for iRef = 1:length(mask.('bandpassFilt'))
  bandpassFilt{iRef} = gpuArray(mask.('bandpassFilt'){iRef});
end

clear refRef refRefRotAvg initPeakMask wdgIMG refIMG


% Reset the geometry so that the  original classes persist.
vol_geometry = vol_geometry_clean;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now we have two montages of references (if gold standard) which later I will
% add in an fsc comparision, which will be used to derive an ideal bandpass
% filter. For now, continue on to use these to further align "sub references"
% using the supplied bandpass filter.
refIMG = cell(2,1);
imgCounts = cell(2,1);

for iGold = 1:2

  if iGold == 1
    halfSet = 'ODD';
  else
    halfSet = 'EVE';
  end
  halfNUM = iGold;

  % Read in the new refs we just created.
  newRefNAME = sprintf('class_%d_Locations_%s_%s', refName, 'REF',halfSet);

  fprintf('reading in new ref %s.\n',halfSet);
  [ refIMG{iGold} ] = BH_unStackMontage4d(1:nRefs, ...
                                    masterTM.(cycleNumber).(newRefNAME){1},...
                                    masterTM.(cycleNumber).(newRefNAME){2},...
                                    sizeWindow);
  
  imgCounts{iGold} = masterTM.(cycleNumber).(newRefNAME){3};
  
  occCell = BH_multi_isCell( refIMG{iGold} );
  nRefs = length(occCell);
  tIMG = cell(nRefs,1);
  for iRef = 1:nRefs
    tIMG{iRef} = refIMG{iGold}{occCell(iRef)};
  end
  refIMG{iGold} = tIMG; clear tIMG
  
  % Save the intermediate "primaryRef" since this will make trouble shooting
  % errors in symmetry/grouping until explicit error checks can be added to
  % BH_parseParameterFile
  CMD = sprintf('mv %s_filtered%d_REF_%s.mrc %s_primaryRef_%d_%s.mrc ', ...
                 outputPrefix, refName, halfSet,outputPrefix, refName, halfSet);
  system(CMD, '-echo');
end


[ refIMG ] = BH_multi_combineLowResInfo( refIMG, imgCounts, pixelSize, maxGoldStandard );

[  nInPlane, inPlaneSearch, angleStep, nAngles] ...
                                        = BH_multi_gridSearchAngles(angleSearch(2,:))

for iGold = 1:2

  if iGold == 1
    halfSet = 'ODD';
  else
    halfSet = 'EVE';
  end
  halfNUM = iGold;
  
  refImg = zeros([size(refIMG{iGold}{1}),nRefs],'single','gpuArray');
  refRotAvg = zeros([size(refIMG{iGold}{1}),nRefs], 'single','gpuArray');
  %refWdg = zeros([size(refIMG{iGold}{1}),nRefs],'single','gpuArray');
  refWdg = NaN
  for iRef = 1:nRefs;
    refRotAvg(:,:,:,iRef) = BH_axialSymmetry(refIMG{iGold}{iRef},120,0,'GPU',[0,0,0]);
    refImg(:,:,:,iRef) = refIMG{iGold}{iRef};
  end

  

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  % Store the cross correlation score, peak location, and wedge weight
  bestAnglesTotal = zeros(nReferences(iGold),10);
  nCount = 1;


    imgClassNAME = sprintf('class_%d_Locations_%s_%s', refName, 'Ref',halfSet)
    refName

    [ classIMG ] = BH_unStackMontage4d(1:refName, ...
                                      masterTM.(cycleNumber).(imgClassNAME){1},...
                                      masterTM.(cycleNumber).(imgClassNAME){2},...
                                      sizeWindow);




  refVector{iGold}(1,:)
  for iClass = refVector{iGold}(1,:)
    tic;
    iClass
   % Load the class onto gpu
    iClassImg = gpuArray(classIMG{iClass});

    iClassWdg = NaN;

 

    % only use the relevant class
    iGroup = refGroup{iGold}(ismember(refVector{iGold},iClass))


    [ alignOUT ] = run_alignment(angleStep, inPlaneSearch, ...
                                 iClassImg, iClassWdg, ...                                          
                                 refImg(:,:,:,iGroup), refWdg, ...
                                 refRotAvg(:,:,:,iGroup), ...
                                 volMask, bandpassFilt, padCalc, padWindow, ...
                                 peakMask, peakCOM, iClass,uniqueSym{iGold});


    bestAnglesTotal(nCount,:) = [gather(iGroup),gather(alignOUT(1,2:10))]; 
 
    timeClass = toc;
    fprintf('finished working on %d/%d classes %s...%fs\n',nCount, ...
                                              length(refVector{iGold}),halfSet,timeClass);
    nCount = nCount + 1;

  end % loop over (classes)

  bestAnglesTotal = gather(bestAnglesTotal)


  [ vol_geometry ] = BH_refAlignmentsApply( vol_geometry, bestAnglesTotal,...
                                        samplingRate,  ...
                                        refVector{iGold}(1,:), refGroup{iGold}(1,:),halfNUM );
  
 
end

  masterTM.(cycleNumber).('RefAlignment') = vol_geometry;
  subTomoMeta = masterTM;
  save(pBH.('subTomoMeta'), 'subTomoMeta');
  
 BH_average3d(PARAMETER_FILE,num2str(CYCLE),'RefAlignment');
 

fprintf('total execution time : %f seconds\n', etime(clock, startTime));



function [ alignOUT ] = run_alignment(angleStep, inPlaneSearch, ...
                               iClassImg, iClassWdg, ...                                          
                               refRef, ref_WDG, refRefRotAvg, ...
                               volMask, bandpassFilt, padCalc, padWindow, ...
                               peakMask, peakCOM, iClass,uniqueSym)

  clear  cccStorage1 cccStorage2 cccStorage3 cccStorage4 
 
  % Out of plane search
  if any(angleStep(2:4))
    flgSearchDepth = 3;

  elseif any(angleStep(5))
    % in plane
    flgSearchDepth = 2;
    peakListTop10(angleStep(5).*nReferences,6) = gpuArray(0);
    nPeak = 1;
    for iRef = 1
      for iPsi = inPlaneSearch
        peakListTop10(nPeak,1) = iRef;
        nPeak = nPeak + 1;
      end
    end  
  else
    error('specify at least an in plane search')
  end

  if (flgSearchDepth == 3)

    [ cccStorage1 ] =  BH_multi_angularSearch( angleStep, 0, 0, ...
                                              iClassImg, iClassWdg, ...                                          
                                              refRefRotAvg, ref_WDG, ...
                                              refRefRotAvg, ...
                                              volMask, bandpassFilt, ...
                                              padCalc, padWindow, ...
                                              peakMask, peakCOM, iClass, ...
                                              uniqueSym);


    cccStorage1(1:10,:)


    % Second loop over top 10 peaks now using the non-rotationally averaged
    % reference and including out of plane angles.

    % peakList is # rows = top peaks
    % reference, phi, theta
    % put zero peak at top of list
    zeroPeak = sortrows(cccStorage1,[3, 4, 5]);
    zeroPeak = [zeroPeak(1,:) ; cccStorage1 ];
    % if zero peak was already there, remove it so no duplicate
    zeroPeak = unique(zeroPeak, 'stable', 'rows');
    peakListTop10 = [zeroPeak(1:10,1),zeroPeak(1:10,3:4),zeroPeak(1:10,8:10)]
  end
  

  [ cccStorage2 ] = BH_multi_angularSearch( angleStep, peakListTop10, ...
                                            inPlaneSearch, ...
                                            iClassImg, iClassWdg, ...                                          
                                            refRef, ref_WDG, ...
                                            refRefRotAvg, ...
                                            volMask, bandpassFilt, ...
                                            padCalc, padWindow, ...
                                            peakMask, peakCOM, iClass, ...
                                            uniqueSym);


  cccStorage2 = unique(cccStorage2((cccStorage2(:,6) ~= 0),:), 'stable','rows');
  if size(cccStorage2, 1) > 9
     cccStorage2(1:10,:)
  else
    cccStorage2
  end
  
  % This is to save time assuming that we can get a good estimate of the
  % particles shift by taking the average of the higher ranking alignments. In
  % testing this was always within ~ half a pixel. If CCC scores are strangely
  % low, suspect this as a break point.


  % Get the top three peaks with unique phi, and theta
  % Return [ref,phi,theta,psi,phistep,thetastep,psistep]
  % Stable prevents any sorting
  if (flgSearchDepth== 3)
    [~,ia,~] = unique(cccStorage2(:,3:4),'stable' ,'rows');
    peakListTop3 = zeros(3,10);
  else
    [~,ia,~] = unique(cccStorage2(:,3:5),'stable' ,'rows');
    peakListTop3 = zeros(3,10);
  end
  
  if numel(ia) >= 10
    TOP = 10;
  else
    TOP = numel(ia);
  end


  for top3 = 1:TOP
    outOfPlaneAngle = cccStorage2(ia(top3),4);
    angleIndex = find(angleStep(:,1)==outOfPlaneAngle,1,'first');
    if (flgSearchDepth == 3 )
      % Search around top 3 +/- 0.5 the original out of plane angular increment
      peakListTop3(top3,:) = [cccStorage2(ia(top3),1), ...
                              cccStorage2(ia(top3),3:5),...
                              angleStep(angleIndex,3)./4,...
                              angleStep(angleIndex,4)./2,...
                              angleStep(angleIndex,5)./2, cccStorage2(ia(top3),8:10)];
    else
      % Search around top 3 +/- 0.5 the original out of plane angular increment
      peakListTop3(top3,:) = [cccStorage2(ia(top3),1), ...
                              cccStorage2(ia(top3),3:5),...
                              0,...
                              0,...
                              angleStep(1,5)./2, cccStorage2(ia(top3),8:10)];
    end
  end
    peakListTop3
  [ cccStorage3 ] = BH_multi_angularSearch( angleStep, peakListTop3, ...
                                            0, ...
                                            iClassImg, iClassWdg, ...                                          
                                            refRef, ref_WDG, ...
                                            refRefRotAvg, ...
                                            volMask, bandpassFilt, ...
                                            padCalc, padWindow, ...
                                            peakMask, peakCOM,iClass, ...
                                            uniqueSym);
    
                                 
  cccStorage3 = unique(cccStorage3((cccStorage3(:,6) ~= 0),:), 'stable','rows');
  if size(cccStorage3, 1) > 9
     cccStorage3(1:10,:)
  else
    cccStorage3
  end
  
  
  if (flgSearchDepth == 3 )
    % Use previous increments/2
    
    % Search around the top peak +/- 0.25 the orginal angular increment
    topPeak = [cccStorage3(1,1), ...
               cccStorage3(1,3:5), ...
               cccStorage3(1,11:13)./2, ...
               cccStorage3(1,8:10)]
  else
    topPeak = [cccStorage3(1,1), ...
               cccStorage3(1,3:5), ...
               0, ...
               0,...
               angleStep(1,5)./3, cccStorage3(1,8:10)]
  end
  
  [ cccStorage4 ] = BH_multi_angularSearch( angleStep, topPeak, ...
                                            0, ...
                                            iClassImg, iClassWdg, ...                                          
                                            refRef, ref_WDG, ...
                                            refRefRotAvg, ...
                                            volMask, bandpassFilt, ...
                                            padCalc, padWindow, ...
                                            peakMask, peakCOM, iClass, ...
                                            uniqueSym);
                                          
  alignOUT = cccStorage4(1,:) 
 


end

