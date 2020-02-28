function [  ] = BH_multi_cRef_VnormApply( matFile, CYCLE, outputPrefix, ...
                                          imgBaseName, ...
                                          idxVect, bFactor, varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

load(matFile);
idxVect = str2num(idxVect);
img1 = sprintf('%s_ODD_NoWgt.mrc',imgBaseName);
img2 = sprintf('%s_EVE_NoWgt.mrc',imgBaseName);
wgt1 = sprintf('%s_ODD_Wgt.mrc',imgBaseName);
wgt2 = sprintf('%s_EVE_Wgt.mrc',imgBaseName);

if nargin > 6
  mFactor = str2num(varargin{1});
else
  mFactor = 0;
end

bFactor = str2num(bFactor);
flgFSCeach = 1;
CYCLE = str2num(CYCLE)
if (CYCLE < 0)
  CYCLE = abs(CYCLE);
  flgFSCeach = 0;
  % Just just the fsc from class 1 for all
end

cycleNumber = sprintf('cycle%0.3u',CYCLE);

h_img = getHeader(MRCImage(img1,0));
h_wgt = getHeader(MRCImage(wgt1,0));

h_img.nZ
img1 = BH_unStackMontage4d(idxVect,img1,h_img.nX/h_img.nZ.*[1,1],'');
img2 = BH_unStackMontage4d(idxVect,img2,h_img.nX/h_img.nZ.*[1,1],'');
wgt1 = BH_unStackMontage4d(idxVect,wgt1,h_wgt.nX/h_wgt.nZ.*[1,1],'');
wgt2 = BH_unStackMontage4d(idxVect,wgt2,h_wgt.nX/h_wgt.nZ.*[1,1],'');


size(img1{1})

for i = idxVect
  
  if (flgFSCeach)
    fscArgs = i;
  else
    fscArgs = 1;
  end

  try
    fscParams = subTomoMeta.(cycleNumber).fitFSC.(sprintf('REF%d',fscArgs));
    aliParams = subTomoMeta.(cycleNumber).fitFSC.(sprintf('ResampleREF%d',fscArgs));
    mskParams = subTomoMeta.(cycleNumber).fitFSC.(sprintf('MaskREF%d',fscArgs));
  catch
    fprintf('\nDid not find REF, looking at Raw assuming this is a class avg\n.');
    try
      fscParams = subTomoMeta.(cycleNumber).fitFSC.(sprintf('Raw%d',fscArgs));
      aliParams = subTomoMeta.(cycleNumber).fitFSC.(sprintf('ResampleRaw%d',fscArgs));
      mskParams = subTomoMeta.(cycleNumber).fitFSC.(sprintf('MaskRaw%d',fscArgs));
    catch
      fprintf('\nJust using fsc info from the Raw and only class 1\n');
      fscParams = subTomoMeta.(cycleNumber).fitFSC.(sprintf('Raw%d',1));
      aliParams = subTomoMeta.(cycleNumber).fitFSC.(sprintf('ResampleRaw%d',1));
      mskParams = subTomoMeta.(cycleNumber).fitFSC.(sprintf('MaskRaw%d',1));
    end
  end

  pixelSize = 0.5/fscParams{4}(end);
  if mFactor
    % Override default mtf filter
    fscParams{3}{3} = mFactor;
  end
  

  [w] = BH_multi_cRef_Vnorm(fscParams, aliParams, mskParams, ...
                            {img1{i},img2{i}},{wgt1{i},wgt2{i}}, ...
                            1,0, pixelSize,bFactor);
                          
  for iBfact = 1:length(bFactor)
    SAVE_IMG(MRCImage(gather(w{iBfact})), ...
             sprintf('%s-n%0.3u-bf-%d.mrc',outputPrefix,i,bFactor(iBfact)),pixelSize);
  end
end

end
