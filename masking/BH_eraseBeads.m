function [ Stack ] = BH_eraseBeads( Stack, beadRadius, fileName, scalePixelsBy, mapBackIter, TLT)
% Erase beads from a non-CTF corrected stack

% I'm a bit conflicted on the timing of the removal. On one hand it seems
% to make more sense to remove the beads after ctf multiplication, this way
% the gold beads are re-localized to the disc we are erasing. This would
% mean that non-CTF corrected template matching would be affected though,
% and also that the ctf estimate will be more strongly biased by the beads,
% which might just be on one surface.


[d1,d2,d3] = size(Stack);
modelName = sprintf('fixedStacks/%s_ali%d.erase',fileName,mapBackIter + 1);


if (mapBackIter)
  % The bead model needs to be updated.
  tiltxf = sprintf('mapBack%d/%s_ali%d_ctf.tltxf',mapBackIter,fileName,mapBackIter);
  old_model = sprintf('fixedStacks/%s_ali%d.erase',fileName, mapBackIter);
  
  % Make sure the erase model is set up for the fixed stack. This model
  % will be updated if tomoCPR is run.
  if isfile(old_model)
    [fail] = system(sprintf('imodtrans -2 %s %s fixedStacks/%s_ali%d.erase', ...
      tiltxf, old_model , fileName, mapBackIter + 1));
    if (fail)
      fprintf('model %s exists : %d\n', old_model,exist(old_model,'file'));
      fprintf('xf %s exists : %d\n', tiltxf,exist(tiltxf,'file'));
      
      error('imodtrans failed to update the bead erase model fixedStacks/%s_ali%d.erase',fileName, mapBackIter)
    end
    
  else
    fprintf('WARNING: skipping bead erasing, b/c no file %s is found',old_model);
    return
  end
  
  
else
  
  % Make sure the erase model is set up for the fixed stack. This model
  % will be updated if tomoCPR is run.
  if ~isfile(sprintf('fixedStacks/%s.erase',fileName))
    fprintf('WARNING: skipping bead erasing, b/c no file fixedStacks/%s.erase is found',fileName);
    return
  end
  
  [fail] = system(sprintf('imodtrans -i fixedStacks/%s.fixed fixedStacks/%s.erase %s',fileName,fileName,modelName));
  if  (fail)
    error('imodtrans failed to set the original bead erase model to the fixed stack header');
  end
  
end

[fail] =  system(sprintf('model2point %s %s_txt',modelName,modelName));

if (fail)
  error('model2point failed to convert the bead erase model to text')
end




beadModel = load(sprintf('%s_txt',modelName));
beadModel(:,1:2) = beadModel(:,1:2) ./ scalePixelsBy;

if isa(Stack(1), 'gpuArray')
  gKernel = BH_mask3d('cylinder',beadRadius.*[2,2]+12,beadRadius.*[1,1],[0,0],'2d');
  useGPU = 1;
  deviceFlag='GPU';
else
  gKernel = gather(BH_mask3d('cylinder',beadRadius.*[2,2]+12,beadRadius.*[1,1],[0,0],'2d'));
  useGPU = 0;
  deviceFlag='cpu';
end

for iPrj = 1:d3
  
  iProjection = Stack(:,:,iPrj);
  iPrj_inFullStack = TLT(iPrj, 23);
  
  oX = round(beadModel( abs(beadModel(:,3) - iPrj_inFullStack +1) < 10^-1, 1));
  oY = round(beadModel( abs(beadModel(:,3) - iPrj_inFullStack +1) < 10^-1, 2));
  
  if (useGPU)
    iMask = zeros(size(iProjection),'single','gpuArray');
  else
    iMask = zeros(size(iProjection),'single');
  end
  
  for iBead = 1:length(oX)
    
    startX = max(1,oX(iBead)-beadRadius-6);
    endX   = min(oX(iBead)+beadRadius+5, size(iProjection,1));
    startY = max(1,oY(iBead)-beadRadius-6);
    endY   = min(oY(iBead)+beadRadius+5, size(iProjection,2));
    
    if (endX - startX <=1) || (endY - startY <= 1)
      fprintf('Skipping bead from %d-%d X, %d-%d, Y iPrj %d\n', ...
        startX, endX, startY, endY, iPrj);
    else
      % Avoid overlap
      currentTile = iMask(startX:endX,startY:endY);
      padVal = BH_multi_padVal(size(currentTile),size(gKernel));
      currentTile = BH_padZeros3d(currentTile,'fwd',padVal,deviceFlag,'single');
      
      beadMask = currentTile > gKernel;
      
      iMask(startX:endX,startY:endY) = BH_padZeros3d((gKernel.*(~beadMask)+beadMask.*currentTile),'inv',padVal,deviceFlag,'single');
    end
    
  end
  
  iMaskInv = 1- iMask;
  % AVG = real(ifftn(fftn(iProjection).*gKernel));
  
  % maxP = max(AVG(:));
  % minP = min(AVG(:));
  % hInc = (maxP-minP)/1000;
  % hVal = minP:hInc:maxP-hInc;
  %
  % h = hist(AVG(:),hVal);
  %
  % dVal = diff(h);
  %
  % cutOFF = find(dVal > 0.33*max(dVal(:)),1,'first');
  %
  % intCutOFF = hVal(cutOFF);
  %
  % iMask = AVG < intCutOFF;
  %
  dataMean = mean(iProjection(iMask<0.01));
  dataRMS  = rms(iProjection(iMask<0.01)-dataMean);
  
  %noiseImage = fftn(iProjection.*iMaskInv);
  %noiseImage = real(ifftn(abs(noiseImage).*exp(1i.* (-pi + 2.*pi.*rand(size(noiseImage),'single','gpuArray')))));
  %Stack(:,:,iPrj) = iProjection.*iMaskInv + iMask.*noiseImage;
  
  Stack(:,:,iPrj) = iProjection.*iMaskInv + iMask.*(randn([d1,d2]).*dataRMS+dataMean);
  
  
  
end

end

