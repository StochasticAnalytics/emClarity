function [ Stack ] = BH_eraseBeads( Stack, beadRadius, beadModel)
% Erase beads from a non-CTF corrected stack


[d1,d2,d3] = size(Stack);


% gKernel = BH_multi_gaussian2d(beadRadius.*[6,6],(1/sigma)*beadRadius,0);

 gKernel = BH_mask3d('cylinder',beadRadius.*[2,2]+12,beadRadius.*[1,1],[0,0],'2d');

for iPrj = 1:d3
  
iProjection = gpuArray(Stack(:,:,iPrj));


oX = round(beadModel( abs(beadModel(:,3) - iPrj +1) < 10^-1, 1));
oY = round(beadModel( abs(beadModel(:,3) - iPrj +1) < 10^-1, 2));

  iMask = zeros(size(iProjection),'single','gpuArray');

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

dataMean = mean(iProjection(iMask<0.01));
dataRMS  = rms(iProjection(iMask<0.01)-dataMean);

if isa(Stack, 'gpuArray')
  Stack(:,:,iPrj) = (iProjection.*iMaskInv + iMask.*(randn([d1,d2],'single','gpuArray').*dataRMS+dataMean));
else
  Stack(:,:,iPrj) = gather(iProjection.*iMaskInv + iMask.*(randn([d1,d2],'single','gpuArray').*dataRMS+dataMean));
end


  
end

end

