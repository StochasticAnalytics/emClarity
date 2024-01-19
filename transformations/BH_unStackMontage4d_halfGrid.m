function [ STACK_OUT ] = BH_unStackMontage4d_halfGrid(  IMAGES , NAME, totalNumber)
%Unstack a 4d stack of 3d images saved as one 3d volume.
%
%
%   Input variables:
%
%   IMAGES = vector with list of images to extract
%
%   NAME = string where montage is.
%
%   LOCATIONS = cell with each i a 6 vector (xlow;xhigh...zhigh) index for
%               extraction.
%
%   Output variables:
%
%   OUTPUT_PREFIX = string specifing PREFIX-mont.mrc
%
%   MONTAGE = cell of mrc images
%
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Goals & limitations
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   TODO:
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



montage = getVolume(MRCImage(NAME));

nVolumes = numel(IMAGES);

STACK_OUT = cell(nVolumes,1);

[ d1 ] = size(montage,1);

if (rem(d1,totalNumber))
  error('The first dimension of the montage is not evenly divisible by the expected number of images.');
else
  d1 = d1 / totalNumber;
end

nIMG = 1;
for iPos = IMAGES
  
  STACK_OUT{nIMG} = montage(1+(iPos-1)*d1:d1 +(iPos-1)*d1,:,:);
  nIMG = nIMG + 1;
  
end

clear montage

end % end of unStackMontage4d function

