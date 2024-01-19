function [ STACK_OUT ] = BH_unStackMontage4d(  IMAGES , NAME, LOCATIONS, sizeWINDOW )
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
ind = LOCATIONS;



% keep original montage size even if a subset of images is selected
if isa(ind, 'cell')
  nVolumes = length(ind);
  sizeRef = ind{1}(2:2:6)';
  
  if any(sizeRef - sizeWINDOW < 0)
    fprintf('SIZE_REF %d %d %d\nSIZE_WINDOW %d %d %d', sizeRef, sizeWINDOW);
    error(['The stored reference is too small to be masked appopriately,' , ...
      'choose a smaller mask radius or find a better solution.\n'])
  else
    refTrim =  BH_multi_padVal( sizeWINDOW, sizeRef );
  end
elseif isnumeric(ind) && numel(ind) == 2
  refTrim = [0,0,0;0,0,0];
  nVolumes = ind(1)*ind(2);
  if rem(size(montage,1),ind(1)) || rem(size(montage,2),ind(2))
    error('Montage dimensions are not divisible by the provided row, col #');
  else
    nX = ind(1); nY = ind(2);
    iDim = size(montage,1)/nX;
    jDim = size(montage,2)/nY;
  end
  ind = cell(nVolumes,1);
  iVol = 1;
  for iY = 1:nY
    for iX = 1:nX
      ind{iVol} = [1 + (iX-1).* iDim,(iX) .* iDim, ...
        1 + (iY-1).* iDim,(iY) .* jDim, ...
        1,size(montage,3) ];
      
      iVol = iVol + 1;
    end
  end
else
  error('Locations must be a cell or 2 vector = MxN');
end

STACK_OUT = cell(nVolumes,1);

for iPos = 1:nVolumes
  
  if ismember(iPos, IMAGES)
    
    
    iIMG = iPos;
    
    STACK_OUT{iPos} = montage(ind{iIMG}(1)+refTrim(1,1) : ind{iIMG}(2)-refTrim(2,1), ...
      ind{iIMG}(3)+refTrim(1,2) : ind{iIMG}(4)-refTrim(2,2), ...
      ind{iIMG}(5)+refTrim(1,3) : ind{iIMG}(6)-refTrim(2,3));
    
    %figure, imshow3D(gather(STACK_OUT{iPos})), pause(5); close(gcf)
  else
    
    STACK_OUT{iPos} = [];
  end
  
  
  
end






end % end of unStackMontage4d function

