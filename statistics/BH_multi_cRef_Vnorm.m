function [weightedImgs] = BH_multi_cRef_Vnorm( ...                                          
                                          fscParams, aliParams, mskParams,...
                                          imgs, weights, ...
                                          flgCombine, flgReference, ...
                                          pixelSize, bFactor,varargin)


% Don't force to zero - particularly for comparison with ground truth, should
% only be used with a very small (i.e. non-zero but no real amplification)
% value.
if any(bFactor < 0)
  bFactor = abs(bFactor)
  noForceMask = 1
else 
  noForceMask = 0
end

if nargin == 10
  gpuDevice(varargin{1})
end

if nargin == 11
  highPassFilter = varargin{2};
else
  highPassFilter = [0,0];
end                                        

[ padVal ] = BH_multi_padVal(size(imgs{1}), size(weights{1}));       
        


% Only minor differences in half sets are expected, so make a mask from combined
% values
nnz(imgs{1}(:));
sum(isnan(imgs{1}(:)))
if (flgCombine < 0)
  flgCombine = 1;
else
  % Resample the odd half which was aligned using odd as ref in fsc calc,
  % but keep an un-aligned copy to apply the changes too.
  img1_orig = imgs{1};
  img1_orig = img1_orig - mean(img1_orig(:));
  img1_orig = img1_orig ./ rms(img1_orig(:));
  imgs{1} = BH_resample3d(gather(img1_orig), ...
                                    aliParams(1,:), ...
                                    aliParams(2,1:3), ...
                                    {'Bah',1,'spline'}, 'cpu', 'forward');
end
                                

% Could replace the spline with Fourier resampling.                            


for iWgt = 1:2
% 
  imgs{iWgt} = imgs{iWgt} - mean(imgs{iWgt}(:));
  imgs{iWgt} = gpuArray(imgs{iWgt} ./rms(imgs{iWgt}(:)));

  [weights{iWgt}, ~] = BH_multi_cRef_wgtCritical(gpuArray(weights{iWgt}));

end





% The mask used in the FSC calc
% % % % % particleMask = BH_mask3d(mskParams{1},2.*mskParams{2},mskParams{3},mskParams{4});

particleMask = BH_mask3d(mskParams{1},mskParams{2},mskParams{3},mskParams{4});

padMask = BH_multi_padVal(size(particleMask), size(imgs{1}));

particleMask = BH_padZeros3d(particleMask,padMask(1,:),padMask(2,:),'GPU','single');
particleMask = BH_multi_randomizeTaper(particleMask);

if (mskParams{5})
  % Soft shape mask used for FSC calculation
%   [ mShape2 ]= BH_mask3d(imgs{1} + imgs{2}, 1.*pixelSize,'','');
% FIXME these should come from the call 
shape_mask_lowpass = 14;
shape_mask_threshold = 2.4;

     
%       padIMG = real(ifftn(padIMG./(padWGT+wienerThreshold)));
    mShape2 = BH_padZeros3d((imgs{1} + imgs{2}),'fwd',padVal,'GPU','single');
    mShape2 = real(ifftn(fftn(mShape2)./(ifftshift(weights{1}+weights{2})+100)));
  [ mShape2 ] = EMC_maskReference(gpuArray(mShape2), pixelSize, {'fsc', true; 'lowpass', shape_mask_lowpass; 'threshold', shape_mask_threshold});  
    mShape2 = BH_padZeros3d(mShape2,'inv',padVal,'GPU','single');
  particleMask = mShape2 .* particleMask;
end







radialGrid = BH_multi_gridCoordinates(size(weights{1}),'Cartesian','GPU', ...
                                      {'none'},1,0,1);
radialGrid = radialGrid ./ pixelSize;  

[ anisoFSC, avgCTF ] = calc_anisoFSC(fscParams, radialGrid,weights, bFactor, pixelSize);



if any(bFactor)  
  [ bFactor, bandFilter ] = calc_bfact(fscParams, radialGrid,flgReference, bFactor, noForceMask,highPassFilter,pixelSize);
else
  bFactor = {1};
end



% With the FSC already compensated for the solvent content, iFpFm = 1
fPfM= 1;


weightedImgs = cell(length(bFactor));
if ~(flgCombine)
  % restore the original orientation of the eve image
  imgs{1} = img1_orig; clear img1_orig
end

for iBfact = 1:length(bFactor)
  if (flgCombine)
    snrWeight = 1;

    [ weightedImgs{iBfact} ] = gather(apply_weights((anisoFSC), ...
                                            (avgCTF{1}+avgCTF{2})./2, ...
                                            (radialGrid < 0.5./pixelSize),...
                                            (imgs{1}+imgs{2}), ...
                                            ifftshift(weights{1}+weights{2}),...
                                            particleMask ,...
                                            padVal, fPfM, bFactor{iBfact}, snrWeight));
  else
    snrWeight = 0.5;
    for iGold = 1:2
          [ weightedImgs{iGold} ] = gather(apply_weights((anisoFSC), ...
                                            (avgCTF{iGold}), ...
                                            (radialGrid < 0.5./pixelSize),...
                                            (imgs{iGold}), ...
                                            ifftshift(weights{iGold}),...
                                            particleMask ,...
                                            padVal, fPfM, bFactor{iBfact}, snrWeight));
    end
  end
  
  
end


end

function [ fsc3D, weights ] = calc_anisoFSC(fscParams, radialGrid, weights, bFactor, pixelSize)


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
    iCone

    iFSCfit = csape(fscParams{1}(:,1),fscParams{2}(:,iCone+1),'variational');

             
    % Don't lowpass in-case of FSC calculation
    if any(bFactor)         
      iForceMask = fscParams{6}{iCone+1}; 
    else
      iForceMask = 1
    end
    
    iFSCclean = csape(osX,fnval(iFSCfit,osX).*iForceMask,'variational');

      

    if ( nCones )
      coneOrientation = BH_defineMatrix(coneList{iCone},'Bah','invVector');
      [ radius,~,height,~,~,~ ] = ...
                                      BH_multi_gridCoordinates( size(radialBinary), ...
                                                                'Cylindrical', ...
                                                                'GPU', ...
                                                                {'single',...
                                                                coneOrientation,...
                                                                [0,0,0]', ...
                                                                'invVector',...
                                                                1,1},...
                                                                0, 0, 0 );
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


  

  fsc3D = fftshift(fsc3D./coneMask);
  fsc3D(fsc3D < 10^-10) = 10^-10;
  


    for i = 1.5:-0.5:0.5
      g = gpuArray(BH_multi_gaussian3d(5.*[1,1,1],i));
      fsc3D = convn(fsc3D,g,'same');
      for iGold = 1:2
        weights{iGold} = convn(weights{iGold},g,'same');                                  
      end
    end
    
  for iGold = 1:2
    weights{iGold} = ifftshift(weights{iGold});
  end
  fsc3D = ifftshift(fsc3D);

  clear radialBinary coneMask 
end

% Should be adapted to work with CONES
function [ avgCTF ] = calc_avgCTF( radialGrid, weights,pixelSize)
  
  radialGrid = fftshift(radialGrid);
  g = gpuArray(BH_multi_gaussian3d(5.*[1,1,1],.75));
  avgCTF = cell(2,1);
  for iGold = 1:2
    bin = floor(size(weights{iGold},1)/4);
    inc = 0.5 / (bin*pixelSize);
    
    avgWgt = zeros(size(weights{iGold}),'single','gpuArray');
  
    for q = 1:bin
       iMask = gpuArray((q-1)*inc <= radialGrid & radialGrid < (q)*inc);
       avgWgt(iMask) =  sum(weights{iGold}(iMask))./(sum(iMask(:)));
    end
    
    avgCTF{iGold} = ifftshift(convn(avgWgt, g, 'same'));
    

    clear tmpIMG
  end
  clear avgWgt iMask radialGrid
end

function [ weightedImg ] = apply_weights(anisoFSC, avgCTF, nyquistLimit, img, ...
                                         wgt, mShape2, ...
                                         padVal, fPfM, bFactor, snrWeight)
                                       


% bFactor contains any sharpening, and/or forced cutoffs
weight = bFactor.* ...
  (wgt + fPfM .*(snrWeight.*2.*anisoFSC./(1-anisoFSC)).^-1 .* avgCTF ).^-1 .* nyquistLimit;

% Prevent very strong low-pass so res-map calc can be run.
if mShape2(1) == -9999
  mShape2 = 1;
  max(weight(:))
  weight = (wgt + 1).^-1;
end



clear bFactor wgt anisoFSC avgCTF nyquistLimit
%img = img ./ max(img(:));


weightedImg = BH_bandLimitCenterNormalize(img-mean(img(:)),weight, ...
                                          '',padVal,'doubleTaper');
                                        
clear img                                        
weightedImg = real(ifftn(weightedImg));

size(weightedImg)
padVal
size(mShape2)
weightedImg = weightedImg(padVal(1,1)+1:end-padVal(2,1),...
                          padVal(1,2)+1:end-padVal(2,2),...
                          padVal(1,3)+1:end-padVal(2,3)) .* mShape2;
                        
clear mShape2
                                         
end

function [ bFactorFilter, bandFilter ] = calc_bfact(fscParams, radialGrid, flgReference, bFactor,noForceMask,highPassFilter,pixelSize)
% First check to see if any cones were calculated in addition to spherical
% shells

nBfactors = length(bFactor);
coneOverlap = zeros(size(radialGrid),'uint8','gpuArray');

if any(fscParams{2}(:,2:end))
  nCones = fscParams{7};
  % Don't include the spherical shells
  coneOffset = 1;
else
  nCones = 0;
  coneOffset = 0;
end

coneList = fscParams{8};
halfAngle= fscParams{9};

bFactorFilter = cell(nBfactors,1);
bandFilter = 0;%zeros(size(radialGrid),'single');
for iBfact = 1:nBfactors
  bFactorFilter{iBfact} = zeros(size(radialGrid),'single');
end
  % Should just add pixel Size to fitFSC.Raw1
  
  osX = fscParams{4};

  mtfX = 0:0.5/(length(osX)-1):0.5;

flgPrintUsage = 1;
for iFilter = 1+coneOffset:1+nCones
 iFilter;

  if (flgReference)
    forceMask = fscParams{5}{iFilter};
  else
    forceMask = fscParams{6}{iFilter};
  end
  if iFilter > 1

    coneOrientation = BH_defineMatrix(coneList{iFilter-1},'Bah','invVector');
    [ radius,~,height,~,~,~ ] = ...
                                BH_multi_gridCoordinates( size(radialGrid), ...
                                                          'Cylindrical', ...
                                                          'GPU', ...
                                                          {'single', ...
                                                          coneOrientation,...
                                                          [0,0,0]', ...
                                                          'invVector',...
                                                          1,1},...
                                                          0, 0, 0 );

    iConeMask = (rad2deg(atan2(radius,abs(height))) < halfAngle{iFilter-1});
    clear radius height
    coneOverlap = coneOverlap + uint8(iConeMask);
  else
    iConeMask = 1;

  end
  
if any(highPassFilter)
  bandPassFilter = (BH_bandpass3d(size(radialGrid),highPassFilter(1),highPassFilter(2),pixelSize,'GPU',pixelSize));
end

 for iBfact = 1:nBfactors
%     bFit  =  fit(osX, exp(bFactor(iBfact).*osX.^2) .* ...
%                  ((exp(-10.*mtfX'.^1.25)+0.06)./1.06).^-1 .* ...
%                  forceMask, 'cubicSpline');

    % detector = -10 is very close to the mtf for a falcon DDD
    
    switch fscParams{3}{3}
      case 0
        adHocMTF = 1;
        if flgPrintUsage 
          fprintf('\n\nUsing MTF 0\n\n');
          flgPrintUsage = 0;
        end
      case 1
        detector = -100;
        capVal = .05;
        adHocMTF =((exp(detector.*osX.^1.25)+capVal)./capVal).^-1;        

        if flgPrintUsage 
          fprintf('\n\nUsing MTF new\n\n');
          flgPrintUsage = 0;
        end     
      case 2
        detector = -20;
        capVal = .13;
        adHocMTF = ((exp(detector.*osX.^1.25)+capVal)./capVal).^-1;
        
        if flgPrintUsage 
          fprintf('\n\nUsing MTF orig\n\n');
          flgPrintUsage = 0;
        end
      otherwise
        detector = -1.*round(fscParams{3}{3});
        capVal = fscParams{3}{3}-round(fscParams{3}{3});
        adHocMTF = ((exp(detector.*osX.^1.25)+capVal)./capVal).^-1;
        
        if flgPrintUsage 
          fprintf('\n\nUsing MTF exp %d %3.3f\n\n',detector,capVal);
          flgPrintUsage = 0;
        end       
        
    end
   
    if (noForceMask)
      bFit  =  fit(osX, exp(bFactor(iBfact)./4.*osX.^2) .* ...
                   adHocMTF,'cubicSpline'); 
    else
      bFit  =  fit(osX, exp(bFactor(iBfact)./4.*osX.^2) .* ...
                   adHocMTF .* forceMask, 'cubicSpline'); 
    end
    
    
    

    bFactorFilter{iBfact} = bFactorFilter{iBfact} + ...
                      iConeMask .* reshape(bFit(radialGrid),size(radialGrid));
  
    if any(highPassFilter)
      bFactorFilter{iBfact} = bFactorFilter{iBfact} .* bandPassFilter;
    end
 
   % if iBfact ==1 
   %   bLimit = fit(osX, forceMask, 'cubicSpline'); 
   %   bandFilter = bandFilter + iConeMask .* reshape(bLimit(radialGrid),size(radialGrid));
   % end

 end
end

%if any(highPassFilt)
%  bandFilter = bandFilter .* BH_bandpass3d(size(bandFilter),highPassFilter(1),highPassFilter(2),pixelSize,'GPU',pixelSize);
%end

if (nCones)
  coneOverlap = single(coneOverlap);
%  figure, imshow3D(coneOverlap);
  divZeroMask = (coneOverlap~=0);
  [ gaussKernel ] = gpuArray(BH_multi_gaussian3d(7, 1.5 ));
  for iBfact = 1:nBfactors
    bFactorFilter{iBfact}(divZeroMask)  = bFactorFilter{iBfact}(divZeroMask) ./coneOverlap(divZeroMask);
    bFactorFilter{iBfact} = ifftshift(gather(convn(fftshift(bFactorFilter{iBfact}), gaussKernel, 'same')));
   % if iBfact == 1
   %   bandFilter(divZeroMask)  =bandFilter(divZeroMask) ./coneOverlap(divZeroMask);
   %   bandFilter = ifftshift(gather(convn(fftshift(bandFilter), gaussKernel, 'same')));
   % end
  end
 
end
end




