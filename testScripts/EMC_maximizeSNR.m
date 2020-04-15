function [ fParticle ] = EMC_maximizeSNR(ref,wgt,fscParams,iRef)




[ bhF, ref{1} ] = fourierTransformer(gpuArray(ref{1}{iRef}));
ref{2} = bhF.fwdFFT(gpuArray(ref{2}{iRef}));

pixelSize = 0.5/fscParams{4}(end);

wienerThreshold = zeros(2,1,'single','gpuArray');
for iGold = 1:2 
  [wgt{iGold}, wienerThreshold(iGold)] = BH_multi_cRef_wgtCritical(gpuArray(wgt{iGold}{iRef}));
       % TODO This is always even so this should be safe. 
   wgt{iGold} = wgt{iGold}([bhF.halfDimSize:size(wgt{iGold},1),1],:,:);
end

LSQ = real(bhF.invFFT(ref{1} ./ ( bhF.swapIndexINV(wgt{1}) + wienerThreshold(1)),2));


[~, MASKCORE ]                      = EMC_maskReference(LSQ , pixelSize, {'fsc',true});
MASKCORE = find(MASKCORE > 0.98);
cccLSQ = (LSQ(MASKCORE) - mean(LSQ(MASKCORE))) ./ std(LSQ(MASKCORE));

  
  
  radialGrid = BH_multi_gridCoordinates(bhF.inputSize,'Cartesian','GPU', ...
                                      {'none'},1,0,1,{'halfGrid'});
  radialGrid = radialGrid ./ pixelSize;  
  

  bFactor = 0;
  
  [ fsc3D, avgCTF ] = calc_anisoFSC(fscParams, radialGrid, wgt, bFactor, pixelSize, bhF);
  
  if any(bFactor)  
    [ bFactor, bandFilter ] = calc_bfact(fscParams, radialGrid,flgReference, bFactor, noForceMask,highPassFilter,pixelSize);
  else
    bFactor = {1};
  end
  
  for iGold = 1:2
    wgt{iGold} = bhF.swapIndexINV(wgt{iGold});
  end
  

  iBfact = 1;
  snrWeight = 0.5;
  search1 =  10 .^ [-4:0];
  ccc = 0.*search1;
  n=1;
  for fPfM = search1
    [ weightedImg ] = gather(apply_weights((fsc3D), ...
                            (avgCTF{2}), ...
                            (radialGrid < 0.5./pixelSize),...
                            (ref{2}), ...
                            wgt{2},...
                            fPfM, bFactor{iBfact}, ...
                            snrWeight,bhF));

    cccVal = cccLSQ .* (weightedImg(MASKCORE) - mean(weightedImg(MASKCORE))) ./ std(weightedImg(MASKCORE));
    ccc(n) = gather(sum(cccVal(:)) ./ numel(MASKCORE));
    n = n + 1; 
  end
  
  [m,c] = max(ccc)
  ccc
  search2 = search1(c-1):search1(c-1):search1(c+1);
  ccc = 0.*search2;
  n=1;
  for fPfM = search2
    [ weightedImg ] = gather(apply_weights((fsc3D), ...
                            (avgCTF{2}), ...
                            (radialGrid < 0.5./pixelSize),...
                            (ref{2}), ...
                            wgt{2},...
                            fPfM, bFactor{iBfact}, ...
                            snrWeight,bhF));

    cccVal = cccLSQ .* (weightedImg(MASKCORE) - mean(weightedImg(MASKCORE))) ./ std(weightedImg(MASKCORE));
    ccc(n) = gather(sum(cccVal(:)) ./ numel(MASKCORE));
    n = n + 1; 
  end
  
  ccc
  [m,c] = max(ccc)
  fParticle = search2(c)
  error('asdf')
end

function [ fsc3D, weights ] = calc_anisoFSC(fscParams, radialGrid, weights, bFactor, pixelSize, bhF)


  nCones = fscParams{7};
  if (nCones)
    firstCone = 1;
  else
    firstCone = 0;
  end

  coneList = fscParams{8};
  halfAngle = fscParams{9};
  samplingRate = fscParams{10};

  osX = fscParams{4};

  radialBinary = (radialGrid < 0.5./pixelSize);
    
  mtfX = 0:0.5/(length(osX)-1):0.5;
  
  
  coneMask = zeros(size(radialBinary) ,'single','gpuArray');
  fsc3D = zeros(size(radialBinary), 'single','gpuArray');
  
  bin = floor(size(weights{1},1)/1);
  inc = 0.5 / (bin*pixelSize);

  for iCone = firstCone:nCones
    fprintf('assembling the 3D FSC from cone %d/%d\n',iCone,nCones);

    iFSCfit = csape(fscParams{1}(:,1),fscParams{2}(:,iCone+1),'variational');

             
    % Don't lowpass in-case of FSC calculation
    if any(bFactor)         
      iForceMask = fscParams{6}{iCone+1}; 
    else
      iForceMask = 1;
    end
    
    iFSCclean = csape(osX,fnval(iFSCfit,osX).*iForceMask,'variational');

      

    if ( nCones )
      coneOrientation = BH_defineMatrix(coneList{iCone},'Bah','invVector');
      [ radius,~,height,~,~,~ ] = ...
                                      BH_multi_gridCoordinates( bhF.inputSize, ...
                                                                'Cylindrical', ...
                                                                'GPU', ...
                                                                {'single',...
                                                                coneOrientation,...
                                                                [0,0,0]', ...
                                                                'invVector',...
                                                                1,1},...
                                                                0, 0, 0, {'halfGrid'});
      iConeMask = (rad2deg(atan2(radius,abs(height))) < halfAngle{iCone});


      coneMask = coneMask + iConeMask;
    else
      iConeMask = radialBinary;
      coneMask = ones(size(radialBinary), 'single', 'gpuArray');
    end
  
 

    iConeMask = iConeMask & radialBinary;
    iFSC = iConeMask; 
    iFSC(iConeMask) = fnval(iFSCclean,radialGrid(iConeMask));
  
    fsc3D = fsc3D + iFSC;
   
   
    
    clear iFSC iConeMask radius height



  end




  fsc3D = fsc3D./coneMask;
  fsc3D(fsc3D < 10^-10) = 10^-10;


  fsc3D = bhF.swapIndexFWD(fsc3D);
  kernelDim = 5;
  padDim = ceil(kernelDim/2)+1;
  
  [ LIMITS ] = EMC_limits(size(fsc3D), size(fsc3D)+padDim, {});
  [ fsc3D ] = EMC_resize(fsc3D, LIMITS, {'taper',false;'value',nan});
  
  for iGold = 1:2
    % The weights start out centered so no need to fwdSwap
    weights{iGold} = EMC_resize(weights{iGold}, LIMITS, {'taper',false});
  end
  
  % Set the padded region to something close by. It may be better to think
  % more about this
  m = isnan(fsc3D(:));
  fsc3D(m) = 0;
  dilationKernel = zeros([1,kernelDim],'single','gpuArray');
  dilationKernel(1) = 1.0;
  tmp = EMC_convn(fsc3D, dilationKernel);
  fsc3D(m) = fsc3D(m) + tmp(m); 
  
  for iGold = 1:2
    tmp = EMC_convn(weights{iGold}, dilationKernel);
    weights{iGold}(m) = weights{iGold}(m) + tmp(m); 
  end
  clear tmp m dilationKernel


  for i = 1.5:-0.5:0.5    
    KERNEL = EMC_gaussianKernel([1,kernelDim], i, 'gpu', {});
    fsc3D = EMC_convn(single(fsc3D), KERNEL);
    for iGold = 1:2
      weights{iGold} = EMC_convn(weights{iGold},KERNEL);                                  
    end
  end
  [ LIMITS ] = EMC_limits(size(fsc3D), size(fsc3D)-padDim, {});
  [ fsc3D ] = EMC_resize(fsc3D, LIMITS, {'taper',false});
  fsc3D = bhF.swapIndexINV(fsc3D);


  for iGold = 1:2
    [ weights{iGold} ] = EMC_resize(weights{iGold}, LIMITS, {'taper',false});
    weights{iGold} = bhF.swapIndexINV(weights{iGold});
  end  

% 
%   clear radialBinary coneMask 
end

function [ weightedImg ] = apply_weights(anisoFSC, avgCTF, nyquistLimit, img, ...
                                         wgt, fPfM, bFactor, snrWeight,bhF)
                                       


% bFactor contains any sharpening, and/or forced cutoffs
weight = bFactor.* ...
  (wgt + fPfM .*(snrWeight.*2.*anisoFSC./(1-anisoFSC)).^-1 .* avgCTF ).^-1 .* nyquistLimit;



weightedImg = real(bhF.invFFT(img .* weight,2));
                                        
                         
end


