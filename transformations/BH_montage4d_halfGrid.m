function [ MONTAGE, IMG_LOC ] = BH_montage4d(  IMAGES , OUTPUT_PREFIX, varargin)
%Create a 4d stack of 3d images saved as one 3d volume.
%
%
%   Input variables:
%
%   IMAGES =  cell with either image names or image volumes
%
%   Output variables:
%
%   MONTAGE = mrc image, or if reading in from disk, an empty string.
%
%   OUTPUT_PREFIX = string specifing PREFIX-mont.mrc
%
%   IMG_LOC = locations for extraction of the images.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Goals & limitations
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   TODO:
%     - check on image sizes
%     - option for padding, for non-square montages,
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~iscell(IMAGES)
  error('Input must be a cell with image names, or images')
end

if nargin > 2
  doHalfGrid = 1;
end

nVolumes = length(IMAGES);



if ischar(IMAGES{1})
  img1 = getVolume(MRCImage(IMAGES{1}));
else
  img1 = IMAGES{1};
end

[dx,dy,dz] = size(IMAGES{1});


MONTAGE = zeros(dx.*nVolumes,dy,dz,'single');

% Track image locations. Return and store in the subTomoMeta file.
imageLocations = cell(nVolumes,1);
for iIMG = 1:nVolumes
  
  if isnumeric(IMAGES{iIMG})
    MONTAGE(1 + (iIMG-1)*dx:dx + (iIMG-1)*dx,:,:) = gather(IMAGES{iIMG});
    IMAGES{iIMG} = [];
  else
    
    MONTAGE(1 + (iIMG-1)*dx:dx + (iIMG-1)*dx,:,:) =  getVolume(MRCImage(IMAGES{iIMG}));
    
  end
  
  
  
  imageLocations{iIMG} = [1 + (iIMG-1)*dx ;dx + (iIMG-1)*dx; ...
    1;dy;1;dz];
end

if ischar(IMAGES{1})
  save(MRCImage(MONTAGE), sprintf('%s-mont.mrc', OUTPUT_PREFIX))
  MONTAGE = '';
end

IMG_LOC = imageLocations;
end % end of montage4d function

