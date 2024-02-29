function [ MONTAGE, IMG_LOC ] = BH_montage4d(  IMAGES , OUTPUT_PREFIX)
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

nVolumes = length(IMAGES);

nDim = ceil(sqrt(nVolumes));


if ischar(IMAGES{1})
  img1 = OPEN_IMG('single', IMAGES{1});
else
  img1 = IMAGES{1};
  
end

imgSize = size(img1);
padSize = 0; % don't change, this is used to access class info, might could store it

clear img1
dx = imgSize(1); dy = imgSize(2); dz = imgSize(3);

ox = 1 + padSize; oy = 1 + padSize; oz = 1;

MONTAGE(nDim.*(dx+padSize)+padSize, nDim.*(dy+padSize)+padSize, dz) = single(0);

% Track image locations. Return and store in the subTomoMeta file.
imageLocations = cell(nVolumes,1);
for iIMG = 1:nVolumes
  
  
  if isnumeric(IMAGES{iIMG})
    img = IMAGES{iIMG};
    
  else
    img = OPEN_IMG('single', IMAGES{iIMG});
    
  end
  
  
  if size(img) ~= imgSize
    error('All subimages must have the same size.)')
  end
  
  if iIMG <= nVolumes
    if mod(iIMG, nDim)
      imageLocations{iIMG} = [ox;ox + dx-1;oy;oy + dy-1;oz; oz + dz-1];
      MONTAGE(ox: ox + dx-1, oy: oy + dy-1, oz: oz + dz-1) = img;
      ox = ox + dx + padSize;
    else
      imageLocations{iIMG} = [ox;ox + dx-1;oy;oy + dy-1;oz; oz + dz-1];
      MONTAGE(ox: ox + dx-1, oy: oy + dy-1, oz: oz + dz-1) = img;
      ox = 1 + padSize ;
      oy = oy + dy + padSize;
    end
    
  end
  
end



if ischar(IMAGES{1})
  save(MRCImage(MONTAGE), sprintf('%s-mont.mrc', OUTPUT_PREFIX))
  MONTAGE = ''
end

IMG_LOC = imageLocations;
end % end of montage4d function

