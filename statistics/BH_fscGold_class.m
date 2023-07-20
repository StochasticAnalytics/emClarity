function [ refWGT ] = BH_fscGold_class( PARAMETER_FILE, CYCLE, STAGEofALIGNMENT,varargin)
%Calculate the fsc and relative power distribution for the two input vol.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Goals & limitations:
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   TODO:
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (nargin < 3 && nargin > 4)
  error('args = PARAMETER_FILE, CYCLE, STAGEofALIGNMENT') 
end


% Explicit reference to location of variables in main memory, or on the GPU.
cpuVar = struct();
GPUVar = struct();

% Should add a max size img line


CYCLE = str2num(CYCLE);

% Put all out put in a subdirectory.
system('mkdir -p FSC');

pBH = BH_parseParameterFile(PARAMETER_FILE);
load(sprintf('%s.mat', pBH.('subTomoMeta')), 'subTomoMeta');
masterTM = subTomoMeta;

cycleNumber = sprintf('cycle%0.3u', CYCLE);
prevCycleNumber = sprintf('cycle%0.3u',CYCLE-1);


flgCones = pBH.('flgCones');
flgClassify= pBH.('flgClassify')
try
  flgMultiRefAlignment = pBH.('flgMultiRefAlignment');
catch
  flgMultiRefAlignment = 0;
end

try
  scaleCalcSize = pBH.('scaleCalcSize');
catch
  scaleCalcSize = 1.5;
end

try
  flgFscShapeMask = pBH.('flgFscShapeMask');
catch
  flgFscShapeMask = 1;
end

try
  shape_mask_lowpass = pBH.('shape_mask_lowpass');
catch
  shape_mask_lowpass = 14; 
end

try
  shape_mask_threshold = pBH.('shape_mask_threshold');
catch
  shape_mask_threshold = 2.4;
end

try
  % Apply the mask with the given parameters, save and exit.
  shape_mask_test = pBH.('shape_mask_test');
catch
  shape_mask_test = false;
end

% Estimating the particle volume still occasionaly goes awry. Place a cap and return a cautionary message.

try 
  minimumParticleVolume = pBH.('minimumparticleVolume');
catch
  minimumParticleVolume = 0.1;
end

try fscWithChimera = pBH.('fscWithChimera');
catch fscWithChimera = 0;
end

outputPrefix = sprintf('./FSC/%s_%s', cycleNumber, pBH.('subTomoMeta')); 
samplingRate = pBH.('Ali_samplingRate');

pixelSize = pBH.('PIXEL_SIZE').*10^10.*samplingRate;
if pBH.('SuperResolution')
  pixelSize = pixelSize * 2;
end



if ( flgCones )
  coneInc = 30;
  nCones = 37;
  calcCones = 1;
else
  coneList = 0;
  nCones = 0;
  halfAngle = 0;
  calcCones = 0;  
end

if ( calcCones )
  n=2;
  coneList = cell(nCones,1);
  halfAngle= cell(nCones,1);
  halfAngle{1} = 1.2*coneInc; % coneInc should be afactor of 90 and 360 so 
                                 % this should return an integer
  coneList{1} = [0,0,0];
  for j = coneInc:coneInc:90;
    for i = 0:coneInc:360-coneInc
      coneList{n} = [i,j,0];
      halfAngle{n}= halfAngle{1} ;
      n=n+1;
    end 
  end
end




% Note this is taken from the class section, not Fsc
refName    = pBH.('Cls_className');% pBH.('Ref_className');

peakSearch   = floor(pBH.('particleRadius')./pixelSize);
peakCOM      =3;

global bh_global_MTF
if isempty(bh_global_MTF)
  bh_global_MTF = 2;
end





% The default is fsc-Gold Standard so the two images should need some degree of
% alignment prior to calculating the fsc. 
flgAlignImages = 1;
flgJustFSC=0;
% check to see if images supplied, or to be read in.
flgEstSNR = 0;
if iscell(STAGEofALIGNMENT)
  if ~isnumeric(STAGEofALIGNMENT{1})
    error('Cell contents are not images.')
  elseif length(STAGEofALIGNMENT) ~= 2
    error('Cell must have two images only.')
  elseif size(STAGEofALIGNMENT{1}) ~= size(STAGEofALIGNMENT{2})
    error('Size of img1 and img2 are inconsistent.')
  else
    %flgAlignImages = 0;
    IMG1 = STAGEofALIGNMENT{1};
    IMG2 = STAGEofALIGNMENT{2}; 
    flgJustFSC=1;
    nReferences=2;
      refVector{1} =1;
      refVector{2}= 1;
      STAGEofALIGNMENT = 'NoAlignment';
      fieldPrefix = 'REF'
  end
else
  
  switch STAGEofALIGNMENT
    case 'RawAlignment'
      savePrefix = 'Raw';
% % % %      if (flgClassify || flgMultiRefAlignment)
       if (flgClassify) 
          fieldPrefix = 'Raw';
       
         className = 0;
         classVector = [0;1];
       else
        className    = pBH.(sprintf('Raw_className'));
        classVector   = pBH.(sprintf('Raw_classes_odd'));
        fieldPrefix = 'REF';      
      end
      


      nReferences = length(classVector(1,:))
      
% % % % 
      imageName{1} =  sprintf('class_%d_Locations_%s_ODD_NoWgt', className,fieldPrefix);
      imageName{2} =  sprintf('class_%d_Locations_%s_EVE_NoWgt', className,fieldPrefix);
      weightName{1} = sprintf('class_%d_Locations_%s_ODD_Wgt', className,fieldPrefix);
      weightName{2} = sprintf('class_%d_Locations_%s_EVE_Wgt', className,fieldPrefix);
      imageName{1}

      refVector{1} =1;
      refVector{2}= 1;
% % % %       nReferences = 1;
      outputPrefix = sprintf('%s_Raw', outputPrefix);

     
    case 'ClassAlignment'
      error('Fsc calculation for class averages is not implemented.')
    case 'NoAlignment'
      savePrefix = 'Raw';
      if (flgClassify)
        
        fieldPrefix = 'NoA';
      else
        fieldPrefix = 'REF';                
      end
      imageName{1} = sprintf('class_0_Locations_%s_ODD_NoWgt', fieldPrefix);
      imageName{2} =  sprintf('class_0_Locations_%s_EVE_NoWgt', fieldPrefix);
      weightName{1} = sprintf('class_0_Locations_%s_ODD_Wgt', fieldPrefix);
      weightName{2} = sprintf('class_0_Locations_%s_EVE_Wgt', fieldPrefix);
      
      nReferences = 1;
      refVector{1} =1;
      refVector{2}= 1;
     
      outputPrefix = sprintf('%s_NoA', outputPrefix);
    case 'Cluster'
      error('Fsc calculation for cluster results is not implemented.')
    case 'RefAlignment'
      savePrefix = 'REF';
      imageName{1} = sprintf('class_%d_Locations_REF_ODD_NoWgt', refName);
      imageName{2} = sprintf('class_%d_Locations_REF_EVE_NoWgt', refName);
      weightName{1} = sprintf('class_%d_Locations_REF_ODD_Wgt', refName);
      weightName{2} = sprintf('class_%d_Locations_REF_EVE_Wgt', refName);
      
      outputPrefix = sprintf('%s_Ref', outputPrefix);
      refVector{1} = 1; %pBH.('ref_Ref_odd')(1,:);
      refVector{2} = 1; %pBH.('ref_Ref_eve')(1,:);
      
      fieldPrefix = 'REF';
      if length(refVector{1}) ~= length(refVector{2})
        error('Fsc ref vectors are not the same length.')
      else
        nReferences = length(refVector{1});
      end
    case 'SnrEstimate'
      savePrefix = 'SNR';
      flgEstSNR = 1;
      if ( CYCLE )
        fieldPrefix = 'Raw'
      else
        fieldPrefix = 'NoA'
      end
      imageName{1} = sprintf('class_%d_Locations_%s_ODD_NoWgt', 25,fieldPrefix);
      imageName{2} = sprintf('class_%d_Locations_%s_EVE_NoWgt', 25,fieldPrefix);  
      outputPrefix = sprintf('%s_Snr', outputPrefix);
      refVector{1} = [1:25];
      refVector{2} = [1:25];
      nReferences = length(refVector{1});
    otherwise
      error('STAGEofALIGNMENT incorrect')
  end
end

for iGold = 1:2
  refSym{iGold} = ones(length(refVector{iGold}));
end


% [ maskType, maskSize, maskRadius, maskCenter ] = ...
%                                   BH_multi_maskCheck(pBH, 'Ali', pixelSize,'FSC')
                               
[ maskType, maskSize, maskRadius, maskCenter ] = ...
                                  BH_multi_maskCheck(pBH, 'Ali', pixelSize)
                                
[ sizeWindow, sizeCalc, sizeMask, padWindow, padCalc] = ...
                                       BH_multi_validArea(  maskSize, maskRadius, scaleCalcSize )
                                     
padDIM = max(max(sizeWindow),384);
padREF = [0,0,0;0,0,0];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (flgAlignImages) && ~(flgJustFSC)
  % Read in the refs
  refIMG = cell(2,1);
  refWGT = cell(2,1);
  for iGold = 1:2

    if iGold == 1
      halfSet = 'ODD';
    else
      halfSet = 'EVE';
    end

    [ refIMG{iGold} ] = BH_unStackMontage4d(1:nReferences, ...
                                masterTM.(cycleNumber).(imageName{iGold}){1},...
                                masterTM.(cycleNumber).(imageName{iGold}){2},...
                                sizeWindow);

    [ refWGT{iGold} ] = BH_unStackMontage4d(1:nReferences, ...
                                masterTM.(cycleNumber).(weightName{iGold}){1},...
                                masterTM.(cycleNumber).(weightName{iGold}){2},...
                                sizeCalc);

    padLSQ = BH_multi_padVal(sizeWindow,sizeCalc);
    trimLSQ = BH_multi_padVal(sizeCalc,sizeWindow);
    for iLSQ = 1:nReferences
      padIMG = fftn(BH_padZeros3d(refIMG{iGold}{iLSQ},padLSQ(1,:),padLSQ(2,:),'GPU','single'));
      padWGT = ifftshift(gpuArray(refWGT{iGold}{iLSQ}));

      [padWGT, wienerThreshold] = BH_multi_cRef_wgtCritical(padWGT);
     
%       padIMG = real(ifftn(padIMG./(padWGT+wienerThreshold)));
      padIMG = real(ifftn(padIMG./(padWGT+100)));

      clear padWGT
      refIMG{iGold}{iLSQ} = BH_padZeros3d(gather(padIMG),trimLSQ(1,:),trimLSQ(2,:),'cpu','single');
      clear padIMG 
    end                          
  end
  clear refWGT
else
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Make a mask, and apply to the average motif && save a masked, binned copy of
% the average for inspection. 

% % % % % % % [ tmpMask ] = BH_mask3d(maskType, sizeMask, maskRadius-7, maskCenter);

[ tmpMask ]  = EMC_maskShape(maskType, sizeMask, maskRadius, 'gpu', {'shift', maskCenter});

volMask{1} = gather(BH_multi_randomizeTaper(tmpMask));
volMask{2} = gather(BH_multi_randomizeTaper(tmpMask)); clear tmpMask

% % % % % % % [ tmpMask] = BH_mask3d(maskType, sizeMask, peakSearch, maskCenter);  
[ tmpMask ]  = EMC_maskShape(maskType, sizeMask, peakSearch, 'gpu', {'shift', maskCenter});

peakMask{1} = gather(BH_multi_randomizeTaper(tmpMask));
peakMask{2} = gather(BH_multi_randomizeTaper(tmpMask)); clear tmpMask

bandpassFilt = cell(2,1);
    bandpassFilt{1} = BH_bandpass3d(sizeCalc,0,0,0,'cpu','nyquist');
    bandpassFilt{2} = BH_bandpass3d(sizeCalc,0,0,0,'cpu','nyquist');

clear radialGrid
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 % Default on.
% if ( flgEstSolvent )
   particleVolume = zeros(nReferences,1);
% else
%   particleVolume = ones(nReferences,1);
% end
fprintf('%d\n',nargin);
if nargin > 3
  sLow = EMC_str2double(varargin{1})
  sTop = EMC_str2double(varargin{2})
else
  sLow = 1e-6;
  sTop = 0.999;
end

if (flgAlignImages) && ~(flgJustFSC) && ~(flgEstSNR)
  refRotAvg = cell(nReferences,1);
  for iGold = 1:2

    for iRef = 1:nReferences
     % iHalf = refVector{iGold}(iRef);
      if iGold == 1
        refRotAvg{iRef} = refIMG{iGold}{iRef};%BH_axialSymmetry(refIMG{iGold}{iRef}, 120,...
                                                                %0, 'GPU',[0,0,0]);
      end


    end
  end
 

    bestAnglesTotal = zeros(nReferences,12);

  nCount = 1;
  
  for iRef = 1:nReferences
    
      fprintf('working on %d/ %d references FscGold\n', iRef, nReferences);


      [shapeMask_1, pV1, particleFraction1, ~] = EMC_maskReference(gpuArray(refIMG{1}{iRef}), pixelSize, {'fsc', true; 'lowpass', shape_mask_lowpass; 'threshold', shape_mask_threshold});  
      [shapeMask_2, pV2, particleFraction2, ~] = EMC_maskReference(gpuArray(refIMG{2}{iRef}), pixelSize, {'fsc', true; 'lowpass', shape_mask_lowpass; 'threshold', shape_mask_threshold});
        
      if (shape_mask_test)
        fprintf('\nSaving your masks and exiting!\n');
        SAVE_IMG(shapeMask_1,sprintf('%s-shape_mask_%2.2f_lowpass_%2.2f_threshold.mrc', ...
                                     outputPrefix, shape_mask_lowpass,shape_mask_threshold),pixelSize);
        return;
      end

      if (flgFscShapeMask)
        shapeMask_1 = gather((shapeMask_1.*volMask{1}).^flgFscShapeMask);
        shapeMask_2 = gather((shapeMask_2.*volMask{2}).^flgFscShapeMask);
      else
        shapeMask_1 = 1;
        shapeMask_2 = 1;
      end

      particleVolume(iRef) = gather(mean([particleFraction1,particleFraction2]));
      pV1 = gather(pV1.*volMask{1});
      pV2 = gather(pV2.*volMask{2});

      % Save a copy with headers set (headers only important for chimera)
      % TODO: check to see if these are used anywhere else, otherwise, they can probably be
      % written only in the if fscWithChimera block.
      eveName = sprintf('%s-eveAli.mrc', outputPrefix);
      oddName = sprintf('%s-oddAli.mrc', outputPrefix);
      SAVE_IMG(MRCImage(gather(refIMG{2}{iRef}.*pV2)), ...
               eveName, 1.0, 1);
      SAVE_IMG(MRCImage(gather(refIMG{1}{iRef}.*pV1)), ...
               oddName, 1.0,1);

      if (fscWithChimera)
        [whereIsChimera, ~] = system('which chimera');
        if (whereIsChimera)
          error('fscWithChimera is called, but "chimera" not in system PATH.')
        end
        writeOutPyAli()
        system(sprintf('chimera --nogui --script "FSC/fitInMap.py %s %s %s-fitInMap.txt" ',...
                       oddName,eveName,outputPrefix));
                             % Read in the results
        % 1:9 rotation matrix 10:12 = dXYZ
        bestAnglesTotal(nCount,:) = load(sprintf('%s-fitInMap.txt',outputPrefix));
        % The results from fit in map are the transpose of Bah, forward rotmat
        bestAnglesTotal(nCount,1:9) = bestAnglesTotal(nCount,[1,4,7,2,5,8,3,6,9]);
        % The results are in Angstrom
        bestAnglesTotal(nCount,(10:12)) = bestAnglesTotal(nCount,(10:12));
      else
        % If we did not align the two halfsets together, we just need to return dummy values
        % for an identity matrix and zero translation.
        bestAnglesTotal(nCount,:) = [1,0,0,0,1,0,0,0,1,0,0,0];
      end

      nCount = nCount+1;
  end
  
elseif (flgAlignImages) && ~(flgJustFSC) && (flgEstSNR)
 error('this block in FSC alignment is slated for removal');
end



gpuDevice(1);

if ( flgEstSNR )
  % Expand cell for all permutations - inefficient memory (2x) but for now makes
  % for easier code...simply treating all as unique classes.
  tmpRef = cell(2,100);
  [rx,ry,rz] = size(refIMG{1}{1});
  % write a continuos block of memory, but don't explicitly fill in all of the
  % zeros.
  for i = 1:2
    for j = 1:100
      tmpRef{i}{j}(rx,ry,rz) = single(0);
    end
  end
  
  t = zeros(2,100);
  
  nRep = 1;
  for iRef = 1:5
    for iRep = 1:4     
      for iPerm = iRep+1:5
        tmpRef{1}{nRep}   = refIMG{1}{iRep  + 5*(iRef-1)};
        tmpRef{2}{nRep}   = refIMG{2}{iPerm + 5*(iRef-1)};
        tmpRef{1}{nRep+1} = refIMG{1}{iPerm + 5*(iRef-1)};
        tmpRef{2}{nRep+1} = refIMG{2}{iRep  + 5*(iRef-1)};
        t(1,nRep) = iRep  + 5*(iRef-1);
        t(2,nRep) = iPerm + 5*(iRef-1);
        t(1,nRep+1) = iPerm + 5*(iRef-1);
        t(2,nRep+1) = iRep  + 5*(iRef-1);
        nRep = nRep +2;
        
      end
    end
  end 
refSymmetry = [1:100;ones(1,100)];
refIMG = tmpRef; clear tmpRef
nReferences = 100;
end

for iRef = 1:nReferences
  
  if (flgAlignImages) && ~(flgJustFSC)


    img2 = refIMG{2}{iRef};

    if (flgEstSNR)
     img1 = BH_resample3d(gather(refIMG{1}{iRef}), ...
                                  bestAnglesTotal(1,3:5), ...
                                  bestAnglesTotal(1,8:10), ...
                                  {'Bah',1,'spline'}, 'cpu', 'forward');
    else

       rotMat = bestAnglesTotal(iRef,1:9);
       dXYZ = bestAnglesTotal(iRef,10:12);

     if abs(3-sum(rotMat([1,5,9]))) < .005
       % Should I pad this? The shifts are small enough it prob is okay.
       fprintf('\nRotation xform is very small, just applying phase shifts\n')
      [ dU, dV, dZ ] = BH_multi_gridCoordinates(size(refIMG{1}{iRef}),'Cartesian','cpu', ...
                                                      {'none'},1,0,0);

        shiftVect = exp((-2i*pi).*(dU.*dXYZ(1) + dV.*dXYZ(2)+ dZ.*dXYZ(3)));
        clear dU dV dZ
        img1 = real(ifftn(fftn(gather(refIMG{1}{iRef})).*shiftVect));
     else
       img1 = BH_resample3d(gather(refIMG{1}{iRef}), ...
                                  rotMat, ...
                                  dXYZ, ...
                                  {'Bah',1,'spline'}, 'cpu', 'forward');
     end
    end
                         
                      
    halfSet = 'GLD';
    fprintf('resampling ref %d.\n', iRef);
    img1 = BH_padZeros3d(img1,[0,0,0],[0,0,0],'cpu','singleTaper');
    img2 = BH_padZeros3d(img2,[0,0,0],[0,0,0],'cpu','singleTaper');
  elseif (flgJustFSC)
    img1=IMG1;
    img2=IMG2;
    
       
      [shapeMask_1, pV1, particleFraction1, ~] = EMC_maskReference(gpuArray(img1), pixelSize, {'fsc', true; 'lowpass', mask_lowpass; 'threshold', mask_threshold});
      [shapeMask_2, pV2, particleFraction2, ~] = EMC_maskReference(gpuArray(img2), pixelSize, {'fsc', true; 'lowpass', mask_lowpass; 'threshold', mask_threshold});
      
      if (shape_mask_test)
        fprintf('\nSaving your masks and exiting!\n');
        SAVE_IMG(shapeMask_1,sprintf('%s-shape_mask_%2.2f_lowpass_%2.2f_threshold.mrc', ...
                                     outputPrefix, shape_mask_lowpass,shape_mask_threshold),pixelSize);
        return;
      end
      
        
      if (flgFscShapeMask)
        shapeMask_1 = gather((shapeMask_1.*volMask{1}).^flgFscShapeMask);
        shapeMask_2 = gather((shapeMask_2.*volMask{2}).^flgFscShapeMask);
      else
        shapeMask_1 = 1;
        shapeMask_2 = 1;
      end
      
      
      particleVolume(iRef) = gather(mean([particleFraction1,particleFraction2]));
      pV1 = gather(pV1.*volMask{1});
      pV2 = gather(pV2.*volMask{2});
     
 

    halfSet = 'OUT';
  else
    halfSet = 'STD';
  end
  
  
    
    
    imgFilt1 = BH_bandLimitCenterNormalize(img1.*volMask{1}.*shapeMask_1,bandpassFilt{iGold}, ...
                                             (volMask{1} > 0.01),padCalc, 'single');
    imgFilt1 = real(ifftn(imgFilt1));
    imgFilt1 = single(gather(imgFilt1(padCalc(1,1)+1 : end - padCalc(2,1), ...
                                      padCalc(1,2)+1 : end - padCalc(2,2), ...
                                      padCalc(1,3)+1 : end - padCalc(2,3) )));
    imgFilt2 = BH_bandLimitCenterNormalize(img2.*volMask{2}.*shapeMask_2,bandpassFilt{iGold}, ...
                                             (volMask{2} > 0.01),padCalc, 'single');
    imgFilt2 = real(ifftn(imgFilt2));
    imgFilt2 = single(gather(imgFilt2(padCalc(1,1)+1 : end - padCalc(2,1), ...
                                      padCalc(1,2)+1 : end - padCalc(2,2), ...
                                      padCalc(1,3)+1 : end - padCalc(2,3) )));                      
      % Should replace these with a call to padZeros3d
    img1 = img1(padWindow(1,1)+1 : end - padWindow(2,1), ...
                padWindow(1,2)+1 : end - padWindow(2,2), ...
                padWindow(1,3)+1 : end - padWindow(2,3) );

    img2 = img2(padWindow(1,1)+1 : end - padWindow(2,1), ...
                padWindow(1,2)+1 : end - padWindow(2,2), ...
                padWindow(1,3)+1 : end - padWindow(2,3) );

  % Save a temp copy of the aligned references to visualy inspect the results
  SAVE_IMG(MRCImage(imgFilt1), sprintf('./FSC/fscTmp_%d_ODD.mrc',iRef));
  SAVE_IMG(MRCImage(imgFilt2), sprintf('./FSC/fscTmp_%d_EVE.mrc',iRef)); 
  
  clear imgFilt1 imgFilt2
  
  [ fscPAD ] = BH_multi_padVal(size(img1), padDIM(1));
  [ rad,~,~,~,~,~ ] = BH_multi_gridCoordinates(padDIM.*[1,1,1], 'Cartesian', 'GPU', ...
                                   {'none'}, 1, 0, 1 );
  rad = single(rad)./pixelSize;
  
  if (flgFscShapeMask)
    
    %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    % TODO add a flag since the phase randomized is not used in practice
    %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
   % Only need to calculate phase randomized masks if the fscShapeMask is
   % applied during the FSC calculation. If instead it is used to estimate
   % the particle volume (flgEstSolvent) then no mask is directly applied.
    fscRandCutoffRes = 3*masterTM.currentResForDefocusError(1);
    lowResShift = pixelSize*2 - 10;
    if lowResShift <= 0
      lowResShift = 0
    else
      lowResShift = 1/(fscRandCutoffRes+lowResShift);
    end
    % Randomize beyond ~ 20A -- calc so that the cutoff is exactly where a FSC
    % shell is bound. the 10 in the divisor is set in the calc_shells function,.
    binDiv = ceil(1.5*padDIM(1)^(1/3));
    shellInc = 0.5/(floor(padDIM(1)/binDiv)*pixelSize);
    randCutoff = floor((1/(fscRandCutoffRes)-lowResShift)/ shellInc) * shellInc;
    fscTcutoff = (floor((1/(fscRandCutoffRes*.95)-lowResShift)/ shellInc)) * shellInc;
    %Calculate the fsc on phase randomized masked volumes.
    [randGrid,~,~,~,~,~ ] = BH_multi_gridCoordinates(padDIM.*[1,1,1], 'Cartesian', ...
                                                       'cpu', {'none'}, 1, 0, 1 );
    randGrid = single(randGrid./pixelSize);
    randLowRES  = (randGrid <  randCutoff);
    randHighRES = (randGrid >= randCutoff);


    clear randGrid
    [ fou1 ] = fftn(BH_padZeros3d(img1, fscPAD(1,:), fscPAD(2,:), 'cpu', 'singleTaper'));
       
    for iShuffle = 1    
      rng('shuffle')
      fou1 = single(real(ifftn(abs(fou1) .* exp(1i .* (...
                        randLowRES  .* angle(fou1) + ...
                        randHighRES .* pi.*(rand(size(fou1),'single')))))));       
    end

   
    
    [ fou2 ] = fftn(BH_padZeros3d(img2 , fscPAD(1,:), fscPAD(2,:), 'cpu', 'singleTaper'));
    for iShuffle = 1
      rng('shuffle')
      fou2 = single(real(ifftn(abs(fou2) .* exp(1i .* (...                    
                        randLowRES  .* angle(fou2) + ...
                        randHighRES .* pi.*(rand(size(fou2),'single'))))))); 
         
    end
 

    
    
    clear randLowRES randHighRES
    fou1 = fftn(fou1.*BH_padZeros3d(pV1, fscPAD(1,:), fscPAD(2,:), 'GPU', 'single'));
    fou2 = fftn(fou2.*BH_padZeros3d(pV2 , fscPAD(1,:), fscPAD(2,:), 'GPU', 'single'));

    [shellsRandFreq, shellsRandFSC, ~,~] = ...
                  calc_shells(fou1, fou2, rad, pixelSize, coneList,'rand');
    clear fou1 fou2
    
  else
    fprintf('masking here\n')
    shapeMask_1 = 1.*volMask{1};
    shapeMask_2 = 1.*volMask{2};
  end

  img1 = img1 - mean(img1(shapeMask_1>0.01));
  img1 = img1 ./ rms(img1(shapeMask_1>0.01));
  img1 = img1 .* shapeMask_1;
  
  img2 = img2 - mean(img2(shapeMask_2>0.01));
  img2 = img2 ./ rms(img2(shapeMask_2>0.01));
  img2 = img2 .* shapeMask_2;
  
  SAVE_IMG(MRCImage(single(gather(img1))), sprintf('./FSC/fscTmp_%d_noFilt_ODD.mrc',iRef));
  SAVE_IMG(MRCImage(single(gather(img2))), sprintf('./FSC/fscTmp_%d_noFilt_EVE.mrc',iRef));
    
  [ img1 ] = BH_padZeros3d(img1 , fscPAD(1,:), fscPAD(2,:), 'GPU', 'singleTaper');
  fou1 = fftn(img1); %clear img1
  [ img2 ] = BH_padZeros3d(img2 , fscPAD(1,:), fscPAD(2,:), 'GPU', 'singleTaper');
  fou2 = fftn(img2); %clear img2
  [shellsFreq, shellsFSC, shellsNUM,shellsPOWER] = ...
                  calc_shells(fou1, fou2, rad, pixelSize,coneList, halfAngle);
  clear fou1 fou2
  if (flgFscShapeMask)
     fou1 = fftn(img1.*BH_padZeros3d(pV1, fscPAD(1,:), fscPAD(2,:), 'GPU', 'single'));
     clear pv1
     fou2 = fftn(img2.*BH_padZeros3d(pV2, fscPAD(1,:), fscPAD(2,:), 'GPU', 'single'));
    [tightFreq, tightFSC,~,~] = ...
                  calc_shells(fou1, fou2, rad, pixelSize,coneList, halfAngle); 
     fitTightFSC = csape(tightFreq(:,1),tightFSC(:,1),'variational');                
     clear pv2
  end
  clear img1 img2

  
  for iFSC = 1:size(shellsFSC,2)
    f = 2.* particleVolume(iRef);
    fscUnMasked = shellsFSC(:,iFSC);
    firstZero = find(shellsFSC(:,iFSC)<0.001,1,'first')-1;
    % The abs() isn't in the paper, but is in the cisTEM code.
    fscParticle = f.*fscUnMasked ./ (1 + (f-1).*abs(fscUnMasked));
    shellsFSC(1:firstZero,iFSC) = fscParticle(1:firstZero);
  end

   % Oversampled curve
  osX = [0:0.001:0.5]'./pixelSize; 
  


  forceMask = cell(nCones+1,1);
  forceMaskAlign = cell(nCones+1,1);
  fitFSC = cell(nCones+1,1);
  fitNUM = cell(nCones+1,1);
  oneBitCut  = zeros(nCones+1,1);
  halfBitCut = zeros(nCones+1,1);


 
    for iCone = 1:nCones+1
      iCone
      fitFSC{iCone} = csape(shellsFreq(:,iCone),shellsFSC(:,iCone),'variational');
      fitNUM{iCone} = csape(shellsFreq(:,iCone),shellsNUM(:,iCone));
    end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% test save masking curve
  if (flgFscShapeMask)


    transitionFreq = find(osX > fscTcutoff,1,'first');
    
    fscRand = csape( shellsRandFreq(:,1), shellsRandFSC(:,1), 'variational');
   
    fscDiff = csape( osX, (fnval(fitTightFSC,osX)-fnval(fscRand,osX))./(1-fnval(fscRand,osX)), 'variational');
    fscTrue = csape( osX, [fnval(fitTightFSC,osX(1:transitionFreq));fnval(fscDiff,osX(transitionFreq+1:end))],'variational');

  SAVE_IMG(MRCImage(single(gather(shapeMask_1))), ...
                sprintf('%s-%d-shapeMask_%d.mrc', outputPrefix, iRef, 1));
   
   SAVE_IMG(MRCImage(single(gather(shapeMask_2))), ...
                 sprintf('%s-%d-shapeMask_%d.mrc', outputPrefix, iRef, 2));
   if (minimumParticleVolume < 1)
     % Only save if used.
     SAVE_IMG(MRCImage(single(gather(pV1))), ...
                sprintf('%s-%d-particleVolEst_%d.mrc', outputPrefix, iRef, 1));

     SAVE_IMG(MRCImage(single(gather(pV2))), ...
                 sprintf('%s-%d-particleVolEst_%d.mrc', outputPrefix, iRef, 2));
   end
           
  end           
 %%%%%%%%%%%%%%%%%%%%%%%%%%%% 


  whiteningFilter = 1;
  % Initial pass at the 1bit and 1/2bit curves, for now set the D/L parameter ast
  % 2/3 as planned for first cycle. Later use (volumeEst)^1/3 as saved from
  % cRef_Vnorm.

  DbyL = 2/3;
  % This should be the same for both half sets
  if ~(flgEstSNR)
    % Otherwise refSymmetry is set in the block where the permutations are
    % created.
    try
      refSymmetry = masterTM.(cycleNumber).('SymmetryApplied').(STAGEofALIGNMENT).(sprintf('%s_className',fieldPrefix)){1}
    catch
      fprintf('\nCould not load the symmetry applied, just using 1\n');
      refSymmetry = ones(2,nReferences);
    end
  end
  nCones
  nEffective = fnval(fitNUM{1},osX).*(3/2.*DbyL).^2 ./ (2.*refSymmetry(2,iRef));

  oneBIT = ( 0.5+2.4142./sqrt(nEffective) ) ./ ...
           ( 1.5 + 1.4142./sqrt(nEffective) );
        

  halfBIT= ( 0.2077+1.9102./sqrt(nEffective) ) ./ ...
           ( 1.2071 + 0.9102./sqrt(nEffective) );
         
  aliBIT = mean([oneBIT,halfBIT]);          
  % Find the two common cutoff values -- need a better way to determine the second
  % value that handles non-monotonic curves and is still smooth/gentle without
  % falling off too slowly.

  
  lowCut1 = find(fnval(fitFSC{1},osX) <= 0.143 & osX > 1/100, 1, 'first');
  try
    oneBitCut(1) = find(fnval(fitFSC{1},osX)-aliBIT < 0 & osX > 1/100, 1, 'first');
  catch
    oneBitCut(1) = find(osX .* pixelSize > 0.425, 1, 'first');
  end
  try
    halfBitCut(1)= find(fnval(fitFSC{1},osX)-halfBIT < 0 & osX > 1/100, 1, 'first');
  catch
    halfBitCut(1) = find(osX .* pixelSize > 0.425, 1, 'first');
  end

% Particularly for working at higher binning, this allows using the full
% frequency range, which is the most information/calc.
if isempty(lowCut1)
  lowCut1 = find(osX .* pixelSize > 0.425, 1, 'first');
end



% The value 0.005 looked steep but not too steep when I did a first test. This
% should be chosen with some more consideration.
forceMask{1} = exp(-0.005.^-2 .* (osX-osX(halfBitCut(1))).^2);
forceMask{1}(1:halfBitCut(1)) = 1;

% Make a more conservative taper for alignment, forcing the cRef curve to zero
% starting from ssnr 2
lowCutAlign = find(fnval(fitFSC{1},osX) <= 1/2 & osX > 1/100, 1, 'first');
% Particularly for working at higher binning, this allows using the full
% frequency range, which is the most information/calc.
if isempty(lowCutAlign)
  lowCutAlign = find(osX .* pixelSize > 0.425, 1, 'first');
end

forceMaskAlign{1} = exp(-0.005.^-2 .* (osX-osX(oneBitCut(1))).^2);
forceMaskAlign{1}(1:oneBitCut(1)) = 1;

lowestRes = 0;
highestRes = 0;

if (flgCones)
  lowestRes  = 1./osX(lowCut1);
  highestRes = 1./osX(lowCut1);
  for iCone = 1:nCones
 
    forceMask{iCone+1} = ones(size(osX));
    lowCut1 = find(fnval(fitFSC{iCone+1},osX) <= 0.143 & osX > 1/100, 1, 'first');
    try
      oneBitCut(iCone+1) = find(fnval(fitFSC{iCone+1},osX)- aliBIT  < 0 & osX > 1/100, 1, 'first');
    catch
      oneBitCut(iCone+1) = find(osX .* pixelSize > 0.425, 1, 'first');
    end
    try
      halfBitCut(iCone+1)= find(fnval(fitFSC{iCone+1},osX)-halfBIT < 0 & osX > 1/100, 1, 'first');
    catch
      halfBitCut(iCone+1) = find(osX .* pixelSize > 0.425, 1, 'first');
    end
    
    % first try to use 0.5, if not default to 0.425 cyc/pix
    if isempty(lowCut1)
      lowCut1 = find(fnval(fitFSC{iCone+1},osX) <= 0.5 & osX > 1/100, 1, 'first');
    end
    if isempty(lowCut1)
      lowCut1 = find(osX .* pixelSize > 0.425, 1, 'first');
    end
 
    if lowestRes < 1./osX(lowCut1)
      lowestRes = 1./osX(lowCut1);
    elseif highestRes > 1./ osX(lowCut1)
      highestRes = 1./osX(lowCut1);
    end
        
    forceMask{iCone+1} = exp(-0.005.^-2 .* (osX-osX(lowCut1)).^2);
    forceMask{iCone+1}(1:lowCut1) = 1;
          
    forceMaskAlign{iCone+1} = exp(-0.005.^-2 .* (osX-osX(oneBitCut(iCone+1))).^2);
    forceMaskAlign{iCone+1}(1:oneBitCut(iCone+1)) = 1;
  end
end

% Calculate cRef as in Rosenthal & Henderson 2003

cRef = cell(nCones+1,1);
cRefAli= cell(nCones+1,1);

cRefCurve = sqrt( abs(2.*fnval(fitFSC{1},osX)./(1+fnval(fitFSC{1},osX))) ) ...
              .* forceMask{1} .* whiteningFilter;
            
cRef{1} = fit(osX, cRefCurve./max(cRefCurve(:)),'cubicSpline');

cRefCurve = sqrt( abs(2.*fnval(fitFSC{1},osX)./(1+fnval(fitFSC{1},osX))) ) ...
              .* forceMaskAlign{1} .* whiteningFilter;
cRefAli{1} = fit(osX,cRefCurve./max(cRefCurve(:)),'cubicSpline');     
if (flgCones)
  for iCone = 1:nCones

    cRefCurve = sqrt( abs(2.*fnval(fitFSC{iCone+1},osX) ./ ...
                    (1+fnval(fitFSC{iCone+1},osX))) ) .* ...
                    forceMask{iCone+1} .* whiteningFilter;
                    
    cRef{iCone+1} = fit(osX,cRefCurve./max(cRefCurve(:)),'cubicSpline');
    
    cRefCurve = sqrt( abs(2.*fnval(fitFSC{iCone+1},osX) ./ ...
                    (1+fnval(fitFSC{iCone+1},osX))) ) .* ...
                    forceMaskAlign{iCone+1} .* whiteningFilter;
                    
    cRefAli{iCone+1} = fit(osX,cRefCurve./max(cRefCurve(:)),'cubicSpline');
                                     
  end
end


% Save the fsc info, for calculating cRef. For some reason, when compiled, fit
% function can run, but storing a cFit object in a struct or cell doesn't seem
% to work?
if ~(flgJustFSC)
  masterTM.(cycleNumber).('fitFSC').(sprintf('%s%d',savePrefix,iRef)) = ...
  {shellsFreq,shellsFSC,{cRef,cRefAli,bh_global_MTF},osX,forceMaskAlign,forceMask,nCones,coneList,halfAngle,samplingRate};

  sprintf('Resample_%s%d',savePrefix,iRef)

  if (flgEstSNR)
    masterTM.(cycleNumber).('fitFSC').(sprintf('fit%s%d',savePrefix,1))=...
                                                   [bestAnglesTotal(1,1:9); ...
                                                    bestAnglesTotal(1,10:12),0,0,0,0,0,0];
  else
    masterTM.(cycleNumber).('fitFSC').(sprintf('Resample%s%d',savePrefix,iRef))=...
                                                   [bestAnglesTotal(iRef,1:9); ...
                                                    bestAnglesTotal(iRef,10:12),0,0,0,0,0,0];
  end

  masterTM.(cycleNumber).('fitFSC').(sprintf('Mask%s%d',savePrefix,iRef)) =  ...
                                               {maskType, sizeMask, ...
                                                maskRadius, maskCenter, ...
                                                1, shape_mask_lowpass, shape_mask_threshold};
                                              
  masterTM.('currentResForDefocusError') = osX(oneBitCut).^-1;
end


imgOUT = STAGEofALIGNMENT;

clear fout famp1 famp2 fphase1 fphase2


    fmid = osX(find(fnval(fitFSC{1},osX) < 0.5 & osX > 1/100, 1, 'first'));
    fgold = osX(find(fnval(fitFSC{1},osX) < 0.143 & osX > 1/100, 1, 'first'));
    
    fprintf('\n0.5 = 1/%f\n0.143 = 1/%f\n', fmid, fgold)

  fscOUT = fopen(sprintf('%s-%d-fsc_%s.txt', outputPrefix, iRef, halfSet),'w');
  
  % This seems wildly unecessary, but adding the new fscFull print section
  % seems to not print right, so using this there as well.
  nCol = length(fitFSC);
  fscMat = zeros(length(osX),nCol+1);
  fscMat(:,1) = osX;
  for iCol = 1:nCol
    fscMat(:,iCol+1) = fnval(fitFSC{iCol},osX);
  end
  fscMat = fscMat';
  nRow = 1; 
  nTot = 1;
  while nTot < numel(fscMat)
    while nRow <= nCol + 1
      fprintf(fscOUT,'%0.6f  ',fscMat(nTot));
      nRow = nRow+1;
      nTot = nTot+1;
    end
    nRow = 1;
    fprintf(fscOUT,'\n');
  end
  fclose(fscOUT);


  if (flgCones)
      
     figure('Visible','off'), plot(osX,fnval(fitFSC{1},osX),'kd','MarkerSize',2.5); hold on;
     plot(osX, oneBIT,'c');     
     plot(osX, halfBIT,'b'); 
     plot(osX, 0.*osX,'k');
     plot(osX,fnval(fitFSC{2},osX),'k--');
      for iCone = 3:length(fitFSC)
        plot(osX,fnval(fitFSC{iCone},osX),'k--');
      end

%     title({'FSC',sprintf('0.5 %3.2f\n0.143 %3.2f',1./fmid,1./fgold)}); 
    title({'FSC',sprintf('0.5 - %3.2f\n0.143 - %3.2f  (%3.2f-%3.2f)\noneBit %3.2f\n halfBit %3.2f\n',1./fmid,1./fgold,lowestRes,highestRes,osX(oneBitCut(1)).^-1,osX(halfBitCut(1)).^-1)}); 

    xlabel('Spatial Freq'); ylabel('fsc');
    legend('FSC','oneBit','halfBit','Location', ...
         'northeast','Orientation','vertical');  
    ylim([-.05 1.025])       
  else
    figure('Visible','off'), plot(osX,fnval(fitFSC{1},osX),'k',...
                                  osX, oneBIT, 'c', ...
                                  osX, halfBIT, 'b', ...                                 
                                  osX, 0.*osX,'k');
    title({'FSC',sprintf('0.5 %3.2f\n0.143 %3.2f\noneBit %3.2f\n halfBit %3.2f\n',1./fmid,1./fgold,osX(oneBitCut(1)).^-1,osX(halfBitCut(1)).^-1)}); 
    xlabel('Spatial Freq'); ylabel('fsc');
    legend('FSC','oneBit','halfBit','Location', ...
               'northeast','Orientation','vertical');
    ylim([-.05 1.025])             
   
    
  end

  file_out = sprintf('%s-%d-fsc_%s', outputPrefix, iRef, halfSet);
  saveas(gcf, file_out,'pdf')
  
  %%%%%%%%%%%%%%%%%%%%%%%%%
  % For comparing results from a tight mask and solvent normalized fsc
  % The tight mask calculation is not so reliable (which is by it isn't used in the first place.
  % This is used for the figure S3 in the Nature Methods paper, but not in regular us.
    figure('Visible','off'), plot(osX,fnval(fitTightFSC,osX),'k-.',...
                                  osX,fnval(fscTrue,osX),'k',...
                                  osX,fnval(fitFSC{1},osX),'b',...
                                  osX, 0.*osX+0.143, 'k--',...
                                  osX,zeros(length(osX)),'k');
                               
    fTightGold = osX(find(fnval(fscTrue,osX) < 0.143 & osX > 1/100, 1, 'first'));  
    
    title({'FSC',sprintf('0.143 - %3.2f, %3.2f\n(paritcleVolume,tightMask)',1./fgold,1./fTightGold)}); 
    xlabel('Spatial Freq'); ylabel('fsc');
    ylim([-.05 1.025])  
    file_out = sprintf('%s-%d-fscFull_%s', outputPrefix, iRef, halfSet);
   % saveas(gcf, file_out,'pdf')
    
   % fscRandOUT = fopen(sprintf('%s-%d-fscFull_%s.txt', outputPrefix, iRef, halfSet),'w');
   % fprintf(fscRandOUT,'%4.4f\t%4.4f\t%4.4f\t%4.4f\n',[osX,fnval(fitTightFSC,osX),fnval(fscTrue,osX),fnval(fitFSC{1},osX)]');
   % fclose(fscRandOUT); 
    figure('Visible','off'), plot(osX,cRef{1}(osX),'kd','MarkerSize',3); hold on;
    if (flgCones)

      plot(osX,cRef{2}(osX),'k--'); 
      for iCone = 3:length(cRef)
        iCone
        plot(osX,cRef{iCone}(osX),'k--');
      end
    end

    title({'cRef',sprintf('0.5 %3.2f\n0.143 %3.2f',1./fmid,1./fgold)}); 
    xlabel('Spatial Freq'); ylabel('fsc');
    legend('cRef','Location', ...
         'northeastoutside','Orientation','vertical');
    ylim([-.05 1.025])

      file_out = sprintf('%s-%d-cRef_%s', outputPrefix, iRef, halfSet);
      saveas(gcf, file_out,'pdf')
  
      figure('Visible','off'), plot(osX,cRefAli{1}(osX),'kd','MarkerSize',3); hold on;
      if (flgCones)

        plot(osX,cRefAli{2}(osX),'k--');
      for iCone = 3:length(cRefAli)
          plot(osX,cRefAli{iCone}(osX),'k--');
      end
      end
      title({'cRefAli',sprintf('0.5 %3.2f\n0.143 %3.2f',1./fmid,1./fgold)}); 
      xlabel('Spatial Freq'); ylabel('fsc');
      legend('cRefAli','Location', ...
           'northeastoutside','Orientation','vertical');
      ylim([-.05 1.025])

    file_out = sprintf('%s-%d-cRefAli_%s', outputPrefix, iRef, halfSet);
    saveas(gcf, file_out,'pdf')  
  try
    fitPower = fit(shellsFreq(:,1).^2,log(shellsPOWER(:,1)),'cubicSpline');
    LR = 10;
    MR = 7;
    HR = min((1.05*fgold)^2, shellsFreq(end-1,1).^2);
    lowRes = find(shellsFreq(:,1).^2 > (1/LR)^2, 1,'first');
    midRes = find(shellsFreq(:,1).^2 > (1/MR)^2, 1,'first');
    endRes = find(shellsFreq(:,1).^2 > HR, 1,'first');
    bFactorFIT1 = fit(shellsFreq(lowRes:midRes,1).^2 , ...
                     log(shellsPOWER(lowRes:midRes,1)),'poly1');
    bFactorFIT2 = fit(shellsFreq(midRes:endRes,1).^2 , ...
                     log(shellsPOWER(midRes:endRes,1)),'poly1');                 

    figure('Visible','off'), plot(osX.^2,fitPower(osX.^2),'k'); hold on
                             plot(osX.^2,bFactorFIT1(osX.^2),'b--');
                             plot(osX.^2,bFactorFIT2(osX.^2),'b--'); 
                             line([(1/LR)^2,(1/LR)^2], ...
                                  [min(log(shellsPOWER(:,1))), ...
                                   max(log(shellsPOWER(:,1)))], ...
                                  'Color','k','LineStyle','--'); 
                             line([(1/MR)^2,(1/MR)^2], ...
                                  [min(log(shellsPOWER(:,1))), ...
                                   max(log(shellsPOWER(:,1)))], ...
                                  'Color','k','LineStyle','--');                            
                             line([HR,HR], ...
                                  [min(log(shellsPOWER(:,1))), ...
                                   max(log(shellsPOWER(:,1)))], ...
                                  'Color','k','LineStyle','--');
  %                            plot(osX.^2,bFactorFIT2(osX.^2),'b--');
              %   outCurve(:,1).^2,outCurve(:,8),'g');

             title({'Guinier Plot',sprintf('\nbFactor(%2.1f-%2.1f-%2.1f)\n %d,%d', ...
                                   LR,MR,sqrt(1./HR),round(bFactorFIT1.p1*-4),...
                                   round(bFactorFIT2.p1*-4))}); ...
                                   xlabel('1/Ang^2'),...                               
                                   ylabel('log(F)');
                                   ylim([0.95*min(log(shellsPOWER(:,1))),...
                                         1.05*max(log(shellsPOWER(:,1)))])

             legend('uncorrected','corrected','Location','northeast',...
                    'Orientation', 'vertical');
    file_out = sprintf('%s-%d-guinier_%s', outputPrefix, iRef, halfSet);
    savefig(gcf,file_out);
    saveas(gcf, file_out,'pdf') 
  catch
    fprintf('\nRan into some error in the guinier analysis.\n');
    fprintf('\nSince this is not critical, skipping and continue.\n');
  end
end

subTomoMeta = masterTM;
save(sprintf('%s.mat', pBH.('subTomoMeta')), 'subTomoMeta');
clearvars -except refWGT
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Make cones - 

function [ coneMask ] = make_cones( coneAngles, padDIM, halfAngle)

% The extra position is for the average of the cone mask, used to consider
% overlapping regions.
% % nCones = length(coneList) + 1;
% % coneMask = cell(nCones);

coneShift = [0;0;0];

coneOrientation = BH_defineMatrix(coneAngles,'Bah','invVector');
[ radius,~,height,~,~,~ ] = ...
                                BH_multi_gridCoordinates( padDIM.*[1,1,1], ...
                                                          'Cylindrical', ...
                                                          'GPU', ...
                                                          {'single',...
                                                          coneOrientation,...
                                                          coneShift, ...
                                                          'invVector',...
                                                          1,1},...
                                                          0, 0, 0 );

coneMask = (rad2deg(atan2(radius,abs(height))) < halfAngle);

  
  

clear radius height 
end % end of the make cones function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calc FSC 

function [shellsFreq, shellsFSC, shellsNsamples, shellsPOWER] = ...
                    calc_shells(fou1, fou2, rad, pixelSize, coneList, halfAngle)
                       

                       
padDIM = size(fou1);

if isempty(fou2)
  flgPowerSpectrum = 1;
else
  flgPowerSpectrum = 0;
end

if isa(coneList,'cell') && ~(flgPowerSpectrum)
  % The last value is just the average, but we calculate the global (spherical)
  % fsc, so we never get to this index.
  nIters = length(coneList)+1;
  %bandLimit = (BH_bandpass3d(padDIM,0,0,2,'GPU',1) == 1);
else
  nIters = 1;
end

if ( flgPowerSpectrum )
  cross = abs(fftn(fou1)); clear fou1
else
  cross = fou1.*conj(fou2);
  auto1 = abs(fou1).^2; clear fou1
  auto2 = abs(fou2).^2; clear fou2
end



% Keep same binDiv but override to calc only the spherical FSC for the
% randomized set.
if ~isa(halfAngle,'cell') 
  if strcmpi(halfAngle, 'rand')
    nIters = 1;
  end
end


binDiv = ceil(1.5*padDIM(1)^(1/3));
if ( flgPowerSpectrum )
  bin = floor(size(cross,1)/3);
else
  bin = floor(size(cross,1)/binDiv);
end
inc = 0.5 / (bin*pixelSize);
shellsFreq = zeros(bin, nIters, 'gpuArray');
shellsFSC = zeros(bin, nIters, 'gpuArray');
shellsNsamples = zeros(bin, nIters, 'gpuArray');
shellsPOWER = zeros(bin, nIters, 'gpuArray');

for iCalc = 1:nIters
  
  if iCalc ==1 
    iConeMask = 1;
  else
%     iConeMask = principleAxes{iCalc -1}; 
    iConeMask = make_cones( coneList{iCalc-1}, padDIM, halfAngle{iCalc-1});
  end
     
  for q = 1:bin
  
  qdep=0;
    iMask = gpuArray((q-1)*inc <= rad & rad < (q+qdep)*inc & iConeMask);
    shellsFreq(q, iCalc) = inc.*(q-1/2);
   
    a = real(sum((cross(iMask))));

    if (flgPowerSpectrum)
      shellsFSC(q, iCalc) = (a ./ sum(iMask(:)));
    else
      b = sum((auto1(iMask)));
      c = sum((auto2(iMask)));
      shellsFSC(q, iCalc) = a ./ sqrt(b * c);
      shellsNsamples(q, iCalc) = sum(iMask(:));   
      shellsPOWER(q,iCalc) = sqrt(b*c);
    end

    
    
  end
  clear iConeMask
end 
shellsFreq = gather(shellsFreq);
shellsFSC = gather(shellsFSC);
shellsNsamples = gather(shellsNsamples);
shellsPOWER = gather(shellsPOWER);
clear cross auto1 auto2 iMask rad principleAxes bandLimit
end % end of the calc_shells function/


function writeOutPyAli()

fOUT = fopen(sprintf('FSC/fitInMap.py'),'w');

fprintf(fOUT,['\n'...
'from sys import argv\n\n',...
'def fit_map_in_map(map1_path, map2_path, xformName,\n',...
'                   initial_map1_transform = None,\n',...
'                   map1_threshold = 3.0,\n',...
'                   ijk_step_size_min = 0.01,\n',...    
'                   ijk_step_size_max = 1.5,\n',...     
'                   max_steps = 5000,\n',...
'                   optimize_translation = True,\n',...
'                   optimize_rotation = True):\n',...
'  from VolumeViewer import open_volume_file\n',...
'  map1 = open_volume_file(map1_path)[0]\n',... 
'  map2 = open_volume_file(map2_path)[0]\n\n',...
'  if initial_map1_transform:\n',...
'    from Matrix import chimera_xform\n',...
'    xf = chimera_xform(initial_map1_transform)\n',...
'    map1.surface_model().openState.globalXform(xf)\n\n',...    
'  use_threshold = (map1_threshold != None)\n\n',...  
'  from FitMap.fitmap import map_points_and_weights, motion_to_maximum\n',...
'  points, point_weights = map_points_and_weights(map1, use_threshold)\n\n',...
'  move_tf, stats = motion_to_maximum(points, point_weights, map2, max_steps,\n',...
'                                     ijk_step_size_min, ijk_step_size_max,\n',...
'                                     optimize_translation, optimize_rotation)\n\n',...
'  import Matrix\n',...
'  if initial_map1_transform:\n',...
'    move_tf = Matrix.multiply_matrices(move_tf, initial_map1_transform)\n\n',...
'  f = open(xformName,''w'')\n',...
'  for i in range(4):\n',...
'    for j in range(3):\n',...
'      f.write(''{:f}''.format(move_tf[j][i]))\n',...
'      f.write("\\n")\n\n',...
'  f.close()\n',...
'  print move_tf\n',...
'  tfs = Matrix.transformation_description(move_tf)\n',...
'  print tfs\n\n',...
't = fit_map_in_map(argv[1],argv[2],argv[3])\n']);

fclose(fOUT);


end

