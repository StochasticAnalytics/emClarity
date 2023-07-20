function [ ] = BH_alignClassRotAvg3d(PARAMETER_FILE, CYCLE)
                               
                                
%Extract and align class averages and references from 4D montages derived.
%
%   Input variables:
%
%   IMAGE = 4d volume, or a string specifing a volume to read in.
%
%   CLASSES = Align subset of class averages. [1, 2, 5, 6]
%                                                      
%   CLASS_NAME = a number (e.g. 64) that refers to the class/montage to use.
%
%   REFERENCES a list same as CLASSES with the class id, and symmetry to apply.
%
%   REF_NAME = a number (e.g. 8) that refers to the class/montage to draw the
%              references from.
%
%   REAL_MASK = {maskType, maskSize, maskRadius, maskCenter}
%
%   BANDPASS =  [HIGH_THRESH, HIGH_CUT, LOW_CUT, PIXEL_SIZE] applied to the
%   particle prior to interpolation.
%
%   ANGLE_SEARCH = [a1 a2 a3 a4 ] the angular search is a grid searched
%                  designed to exhaustively sample a unit sphere over the range
%                  you specify. Out of plane +/- a1 in steps of size a2, in
%                  plane +/- a3 in steps of a4. These together define a set of
%                  "latitudes" if you will, and the in plane sampling at each
%                  point sampled on that latitude, while the longitudinal
%                  sampling is calculated to be evenly sampled at the same rate
%                  as the latitude. The smaller the out of plane step size, the
%                  more longitudinal sampling points, and also the larger the
%                  absolute out of plane angle, the greater number of steps it
%                  takes to make a full revolution.
%
%   PEAK = [x y z] = *RADIUS* of peak search
%   PEAK_MASS = [x y z] = *RADIUS* for center of mass search around max peak
%         Set to a zero to ignore either option.
%
%   REFERENCES = a cell with list of images to use as references.
%
%   GEOMETRY = A structure with tomogram names as the field names, and geometry
%              information in a 26 column array.
%              Additionally, a field called 'source_path' has a value with the
%              absolute path to the location of the tomograms.
%
%              The input is a string 'Geometry_templatematching.mat' for
%              example, and it is expected that the structure is saved as the
%              variable named geometry.
%
%
%   OUTPUT_PREFIX = String to prepend to output volumes.
%
%
%   Output variables:
%
%   None = files are written to disk in the current directory.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Goals & Limitations:
%
%   Align class averages to a reference. Here it is implicitly assumed that the
%   references and the volumes they will be aligned against are the same
%   dimension. I am not going to remove all of the information and steps related
%   to binning, as these will be needed when writing the function to handle
%   alignment of raw subTomos, but later I will clean this up.
%
%   Assumed to run on GPU.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   TODO

%     - Update geometry to record % sampling
%     - deal with assumption that 256 256 256 is sufficient for volumes.
%     - Add error check for mask size/ v radius
%     - Change angular searches to allow for a translational only search.
%     - Work through angular sampling to be sure it is doing what you think it
%     is.
%     - Read sampling from metadata to deal with bandpass.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (nargin ~= 2)
  error('args = PARAMETER_FILE, CYCLE')
end

PRECISION = 'single';

startTime =  clock;
CYCLE = EMC_str2double(CYCLE); 


cycleNumber = sprintf('cycle%0.3u', CYCLE);

pBH = BH_parseParameterFile(PARAMETER_FILE);
load(sprintf('%s.mat', pBH.('subTomoMeta')), 'subTomoMeta');

maxGoldStandard = subTomoMeta.('maxGoldStandard');

bFactor = pBH.('Fsc_bfactor');

try
  flgRotAvgRef    = pBH.('flgRotAvgRef');
catch
  flgRotAvgRef = 0;
end
samplingRate   = pBH.('Ali_samplingRate');
pixelSize      = pBH.('PIXEL_SIZE').*samplingRate.*10^10;
if pBH.('SuperResolution')
  pixelSize = pixelSize * 2;
end
try
  scaleCalcSize = pBH.('scaleCalcSize');
catch
  scaleCalcSize = 1.5;
end

angleSearch  = pBH.('Cls_angleSearch');
refName    = pBH.('Cls_className'); %pBH.('Ref_className');
className = pBH.('Cls_className');
peakSearch   = floor(pBH.('particleRadius')./pixelSize)
peakCOM      = [1,1,1].*peakCOM;
outputPrefix = sprintf('%s_%s', cycleNumber, pBH.('subTomoMeta'));



flgAngleShift{1}= pBH.('ref_AngleShift_odd');
flgTransShift{1}= pBH.('ref_TransShift_odd');
flgRefRef{1}    = pBH.('ref_Ref_odd');
refVectorFull{1}= pBH.('Ref_references_odd');
classVector{1}  = pBH.('Cls_classes_odd')(1,:);
features{1}     = pBH.('Pca_coeffs_odd');
vol_geometry{1} = subTomoMeta.(cycleNumber).('ClusterResults').( ...
                                         sprintf('%s_%d_%d_nClass_%d_ODD', ...
                                         outputPrefix,features{1}(1,1), ...
                                         features{1}(1,end), className));

flgAngleShift{2}= pBH.('ref_AngleShift_eve');
flgTransShift{2}= pBH.('ref_TransShift_eve');
flgRefRef{2}    = pBH.('ref_Ref_eve');
refVectorFull{2}= pBH.('Ref_references_eve');
classVector{2}  = pBH.('Cls_classes_eve')(1,:);
features{2}     = pBH.('Pca_coeffs_eve');
vol_geometry{2} = subTomoMeta.(cycleNumber).('ClusterResults').( ...
                                         sprintf('%s_%d_%d_nClass_%d_EVE', ...
                                         outputPrefix,features{2}(1,1), ...
                                         features{2}(1,end), className));
% Merge alignments that could come from using different feature vectors in the
% classification, back into one metadata.
vol_geometry = BH_mergeClassGeometry(vol_geometry{1},vol_geometry{2});

                                         




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


pathList= subTomoMeta.mapPath;
extList = subTomoMeta.mapExt;
masterTM = subTomoMeta; clear subTomoMeta

% make sure the number of references match the unique groups in the classVector
% and also that the class/group pairs match the class/ref pairs.
nReferences(1:2) = [length(unique(refGroup{1})),length(unique(refGroup{1}))];
nReferences = nReferences .* [~isempty(refGroup{1}),~isempty(refGroup{2})]

uniqueSym = cell(2,1);
for iGold = 1:2
  [~,uniqueGroup,~] = unique(refGroup{iGold});
  uniqueSym{iGold} = refSym{iGold}(uniqueGroup)
end






[ maskType, maskSize, maskRadius, maskCenter ] = ...
                                  BH_multi_maskCheck(pBH, 'Cls', samplingRate);
[ sizeWindow, sizeCalc, sizeMask, padWindow, padCalc ] = ...
                                       BH_multi_validArea(  maskSize, maskRadius, scaleCalcSize )
padREF = [0,0,0;0,0,0];                                     
if any(peakSearch > maskRadius)
  fprintf('\n\n\tpeakRADIUS should be <= maskRADIUS!!\n\n')
  peakSearch( (peakSearch > maskRadius) ) = ...
                                        maskRadius( (peakSearch > maskRadius) );
end


% Read in the references.
refIMG = cell(2,1);
imgCounts = cell(2,1);

for iGold = 1:2

  if iGold == 1
    halfSet = 'ODD';
  else
    halfSet = 'EVE';
  end

  imgNAME = sprintf('class_%d_Locations_REF_%s', refName, halfSet);
  weightNAME = sprintf('class_%d_Locations_REF_%s_Wgt', refName, halfSet);
  
  imgCounts{iGold} = masterTM.(cycleNumber).(imgNAME){3};

  [ refIMG{iGold} ] = BH_unStackMontage4d(1:nReferences(iGold), ...
                                   masterTM.(cycleNumber).(imgNAME){1}, ...
                                   masterTM.(cycleNumber).(imgNAME){2},...
                                   sizeWindow);
                                 
  [ refWDG{iGold} ] = BH_unStackMontage4d(1:nReferences(iGold), ...
                                masterTM.(cycleNumber).(weightNAME){1},...
                                masterTM.(cycleNumber).(weightNAME){2},...
                                sizeCalc);                                 


  sizeREF = masterTM.(cycleNumber).(imgNAME){2}{1};
  sizeREF = sizeREF(2:2:6)'     
  
   
end

[ refIMG ] = BH_multi_combineLowResInfo( refIMG, imgCounts, pixelSize, maxGoldStandard);


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


[ volMask ] = gpuArray(BH_mask3d(maskType, sizeMask, maskRadius, maskCenter));
volBinary = (volMask > 0.01);

[ peakMask] = gpuArray(BH_mask3d(maskType, sizeMask, peakSearch, maskCenter));                                                   
peakBinary = (peakMask > 0.01);

bandpassFilt = cell(nReferences(1),1);
  [radialGrid,~,~,~,~,~ ] = BH_multi_gridCoordinates(sizeCalc, 'Cartesian', ...
                                                     'cpu', {'none'}, 1, 0, 1 );
  radialGrid = single(radialGrid./pixelSize);
for iRef = 1:nReferences(1)
  

  fscINFO = masterTM.(cycleNumber).('fitFSC').(sprintf('REF%d',iRef));

  % The class averages have roughly the same SNR as the references  so apply any
  % bFactor to them as well.
  [ ~, bandpassFilt{iRef} ] =  BH_multi_cRef( fscINFO, radialGrid, bFactor,1 ); 
  bandpassFilt{iRef} = gpuArray(bandpassFilt{iRef});
end 

bestAnglesResults = cell(2,1);


try
  EMC_parpool(2)
catch
  delete(gcp('nocreate'))
  EMC_parpool(2)
end

parfor iGold = 1:2

  if iGold == 1
    halfSet = 'ODD';
  else
    halfSet = 'EVE';
  end

    
    ref_FT       = zeros([sizeWindow,nReferences(iGold)], PRECISION, 'gpuArray');
    refRotAvg_FT = zeros([sizeWindow,nReferences(iGold)], PRECISION, 'gpuArray');


    refRotAvg = refIMG{iGold};
    refTrans = cell(length(refIMG{iGold}).*2);
    nRefOut = 1;
    for iRef = 1:nReferences(iGold)
      

      if (flgRotAvgRef)
        refRotAvg{iRef} = gather( ...
          BH_axialSymmetry(gpuArray(refIMG{iGold}{iRef}),120, 0, ...
                                                                'GPU', [0,0,0]));
      
      % For later stages where the angles are very small, still scan only
      % azimuthal and out of plane, but don't rotationally average the in-plane
      % angles.        
      else
        refRotAvg{iRef} = refIMG{iGold}{iRef};
  
        
      end
      
      padTransTrim = padWindow + padREF;
      refTransTrim = refIMG{iGold}{iRef}(padTransTrim(1,1)+1 : end - padTransTrim(2,1), ...
                                          padTransTrim(1,2)+1 : end - padTransTrim(2,2), ...
                                          padTransTrim(1,3)+1 : end - padTransTrim(2,3) );

      refTrans{nRefOut} = real(ifftn( BH_bandLimitCenterNormalize(refTransTrim.*volMask, ...
        bandpassFilt{iRef}, volBinary,padCalc,'double')));
      
      refTrans{nRefOut} = gather(refTrans{nRefOut}(padCalc(1,1) + 1: end - padCalc(2,1),...
        padCalc(1,2) + 1: end - padCalc(2,2),...
        padCalc(1,3) + 1: end - padCalc(2,3)) .*volMask);
      
      refTrans{nRefOut+1} = real(ifftn( BH_bandLimitCenterNormalize(refTransTrim.*peakMask, ...
        bandpassFilt{iRef}, peakBinary,padCalc,'double')));
      
      refTrans{nRefOut+1} = gather(refTrans{nRefOut+1}(padCalc(1,1) + 1: end - padCalc(2,1),...
        padCalc(1,2) + 1: end - padCalc(2,2),...
        padCalc(1,3) + 1: end - padCalc(2,3)) .*peakMask);
      
      nRefOut = nRefOut +2;
      

iGold
iRef
size(refIMG{iGold}{iRef})
size(ref_FT)
      %winRefDiff = padREF - padWindow;
      ref_FT(:,:,:,iRef) = refIMG{iGold}{iRef}( ...
                                  padREF(1,1)+1 : end - padREF(2,1), ...
                                  padREF(1,2)+1 : end - padREF(2,2), ...
                                  padREF(1,3)+1 : end - padREF(2,3) );
      refRotAvg_FT(:,:,:,iRef) = refRotAvg{iRef}( ...
                                  padREF(1,1)+1 : end - padREF(2,1), ...
                                  padREF(1,2)+1 : end - padREF(2,2), ...
                                  padREF(1,3)+1 : end - padREF(2,3) );
                                
      

     

      
    end
    refMontage = BH_montage4d(refTrans,'');
    SAVE_IMG(MRCImage(refMontage),sprintf('%s_class_refFiltered_%s.mrc',cycleNumber,halfSet));
    % Read in the class averages.
%     refMontage{iGold} = BH_montage4d(refTrans,'');   
    imgClassNAME = sprintf('class_%d_Locations_%s_%s_NoWgt', className, 'Cls', halfSet);
    wdgClassNAME = sprintf('class_%d_Locations_%s_%s_Wgt', className, 'Cls', halfSet);
    
    [ classIMG ] = BH_unStackMontage4d(classVector{iGold}, ...
      masterTM.(cycleNumber).(imgClassNAME){1}, ...
      masterTM.(cycleNumber).(imgClassNAME){2},sizeWindow);
    
    [ classWDG ] = BH_unStackMontage4d(classVector{iGold}, ...
      masterTM.(cycleNumber).(imgClassNAME){1}, ...
      masterTM.(cycleNumber).(imgClassNAME){2},sizeWindow);
  

    nClassesPossible = length(classVector{iGold})
    nClasses = length(classIMG)
    
    %%%%%%%%%%%%%%%%%%%%% Determine the angular search, if any are zero, don't
    %%%%%%%%%%%%%%%%%%%%% search at all in that dimension.
    [  nInPlane, inPlaneSearch, angleStep, nAngles] ...
      = BH_multi_gridSearchAngles(angleSearch)
    
    fprintf('%d ',inPlaneSearch);
    fprintf('\n');
    
    % Store the cross correlation score, peak location, and wedge weight
    bestAnglesTotal = zeros(nClassesPossible,10);
    nCount = 1;
    
    
    for iClass = classVector{iGold}
      
      tic;
      
      
      % Load the class into gpu, center and normalize
      iClassImg = gpuArray(classIMG{iClass}( ...
                                  padREF(1,1)+1 : end - padREF(2,1), ...
                                  padREF(1,2)+1 : end - padREF(2,2), ...
                                  padREF(1,3)+1 : end - padREF(2,3) ));
                                

      iClassWdg = ifftshift(gpuArray(classWDG{iClass}( ...
                                  padREF(1,1)+1 : end - padREF(2,1), ...
                                  padREF(1,2)+1 : end - padREF(2,2), ...
                                  padREF(1,3)+1 : end - padREF(2,3) )));
  
   
      
      %   [ iClassImg ] = BH_padZeros3d(iClassImg, padPre, padPost, 'GPU', 'single');
      %[ iClassImg ] = BH_bandLimitCenterNormalize(iClassImg, bandpassFilt, volMask);
      %  iClassImg   = real(ifftn(iClassImg));
      % First loop over all out of plane, no in plane, with rotationally averaged
      % reference.
      cccStorage1 = [];
      cccStorage2 = [];
      cccStorage3 = [];
      cccStorage4 = [];
      % Out of plane search
      if any(angleStep(2:4))
        flgSearchDepth = 3;
      elseif any(angleStep(5))
        % in plane
        flgSearchDepth = 2;
        peakListTop10(angleStep(5).*nReferences(iGold),6) = gpuArray(0);
        
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
        % Note the rotationally averaged ref is passed as main ref
        % 
        [ cccStorage1 ] =  BH_multi_angularSearch( angleStep, 0, 0, ...
                                                  iClassImg, iClassWdg, ...
                                                  refRotAvg_FT, NaN, ...
                                                  refRotAvg_FT, ...
                                                  volMask, bandpassFilt, ...
                                                  padCalc, padWindow,...
                                                  peakMask, peakCOM, iClass, ...
                                                  uniqueSym{iGold});
        
        
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
                                                ref_FT, refWDG{iGold}, ...
                                                refRotAvg_FT, ...
                                                volMask, bandpassFilt, ...
                                                padCalc, padWindow, ...
                                                peakMask, peakCOM,iClass, ...
                                                uniqueSym{iGold});
      
      
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
                                                ref_FT, refWDG{iGold}, ...
                                                refRotAvg_FT, ...
                                                volMask, bandpassFilt, ...
                                                padCalc, padWindow, ...
                                                peakMask,peakCOM,iClass, ...
                                                uniqueSym{iGold});

      
      cccStorage3 = unique(cccStorage3((cccStorage3(:,6) ~= 0),:), 'stable','rows');
      cccStorage3(1:10,:)
      
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
                                                ref_FT, refWDG{iGold}, ...
                                                refRotAvg_FT, ...
                                                volMask, bandpassFilt, ...
                                                padCalc, padWindow,...
                                                peakMask, peakCOM,iClass, ...
                                                uniqueSym{iGold});

      cccStorage4(1,:)
      
      
      bestAnglesTotal(iClass,:) = gather(cccStorage4(1,1:10));
      
      timeClass = toc;
      fprintf('finished working on %d/%d classes...%fs\n',nCount,nClassesPossible,timeClass)
      nCount = nCount + 1;
      
      % Save alignment incase of crash, long runs can be resumed.
      
%%%      save('bestAnglesTotalClass.mat', 'bestAnglesTotal');
      


      
    end % loop over (classes)
   
    bestAnglesResults{iGold} = gather(bestAnglesTotal);
    % Save a text copy of results.
    resultsOut = fopen(sprintf('%s_bestAngles_%s.txt', cycleNumber, halfSet),'w');
    fprintf(resultsOut,'%d %d %3.3f %3.3f %3.3f %1.6f %d %4.4f %4.4f %4.4f\n',bestAnglesTotal');
    fclose(resultsOut);
end

for iGold = 1:2

    if iGold == 1
    halfSet = 'ODD';
    else
      halfSet = 'EVE';
    end
      imgClassNAME = sprintf('class_%d_Locations_%s_%s', className, 'Cls', halfSet);
      [ classIMG ] = BH_unStackMontage4d(classVector{iGold}, ...
      masterTM.(cycleNumber).(imgClassNAME){1}, ...
      masterTM.(cycleNumber).(imgClassNAME){2},sizeWindow);
  
    [ vol_geometry ] = BH_classAlignmentsApply( vol_geometry, bestAnglesResults{iGold},...
                                                 samplingRate,classVector{iGold},iGold,halfSet);

 
 
   % SAVE_IMG(MRCImage(refMontage{iGold}),sprintf('%s_class_refFiltered_%s.mrc',cycleNumber,halfSet));
  % Max a montage of the applied corrections for viewing.
  % Insert blank images so subsets are still corresponding.
  classInterp = cell(className,1);
  emptyImg = zeros(sizeWindow, 'single');
  for iClass = 1:className
    if ismember(iClass, classVector{iGold})
      iClass
      
      iClassImg = single(gpuArray(classIMG{iClass}));
      sizeMask
      size(iClassImg)
      %figure, imshow3D(gather(iClassImg)) ; pause(2) ; close(gcf)
      % Note 'forward' affects the shifts but not the angles, a value of 10  means the
      % class was 10 away from the ref
      idx = find(bestAnglesResults{iGold}(:,2) == iClass);
      iClassImg = BH_resample3d(iClassImg, bestAnglesResults{iGold}(idx, 3:5), ...
                                           bestAnglesResults{iGold}(idx, 8:10), ...
                                                           'Bah', 'GPU', 'inv');
    else

      iClassImg = emptyImg;
    end

    classInterp{iClass} = gather(single(iClassImg(padWindow(1,1)+1 : end - padWindow(2,1), ...
                                                  padWindow(1,2)+1 : end - padWindow(2,2), ...
                                                  padWindow(1,3)+1 : end - padWindow(2,3)) ));

   
  end


  try
    [montOUT] = BH_montage4d(classInterp, '');
    imout = sprintf('%s_class%d_%s_aligned.mrc',outputPrefix, className,halfSet);
    SAVE_IMG(MRCImage(gather(montOUT)), imout);
  catch
    fprintf('error in montaging the aligned classes line 603.\n')
  end
end

masterTM.(cycleNumber).('ClassAlignment') = vol_geometry;
subTomoMeta = masterTM;
save(pBH.('subTomoMeta'), 'subTomoMeta');   

fprintf('Total execution time : %f seconds\n', etime(clock, startTime));
gpuDevice(1);

try
  delete(gcp('nocreate'))
catch
end
end % end of average3d function


