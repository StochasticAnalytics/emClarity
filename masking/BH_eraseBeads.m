function [ Stack ] = BH_eraseBeads( Stack, beadRadius, beadModel)
% Erase beads from a non-CTF corrected stack


[d1,d2,d3] = size(Stack);

% gKernel = BH_multi_gaussian2d(beadRadius.*[6,6],(1/sigma)*beadRadius,0);

if isa(Stack(1), 'gpuArray')
  gKernel = BH_mask3d('cylinder',beadRadius.*[2,2]+12,beadRadius.*[1,1],[0,0],'2d');
  useGPU = 1;
else
  gKernel = gather(BH_mask3d('cylinder',beadRadius.*[2,2]+12,beadRadius.*[1,1],[0,0],'2d'));
  useGPU = 0;
end

for iPrj = 1:d3
  
iProjection = Stack(:,:,iPrj);


oX = round(beadModel( abs(beadModel(:,3) - iPrj +1) < 10^-1, 1));
oY = round(beadModel( abs(beadModel(:,3) - iPrj +1) < 10^-1, 2));

if (useGPU)
  iMask = zeros(size(iProjection),'single','gpuArray');
else
  iMask = zeros(size(iProjection),'single');
end

  for iBead = 1:length(oX)
    try
      % Avoid overlap
      currentTile = iMask(oX(iBead)-beadRadius-6:oX(iBead)+beadRadius+5, ...
            oY(iBead)-beadRadius-6:oY(iBead)+beadRadius+5);

      beadMask = currentTile > gKernel;

      iMask(oX(iBead)-beadRadius-6:oX(iBead)+beadRadius+5, ...
            oY(iBead)-beadRadius-6:oY(iBead)+beadRadius+5) = ...
               (gKernel.*(~beadMask)+beadMask.*currentTile);
    catch
      fprintf('Skipping edge on bead %d prj %d\n',iBead, iPrj);
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


Stack(:,:,iPrj) = iProjection.*iMaskInv + iMask.*(randn([d1,d2]).*dataRMS+dataMean);


  
end

end

