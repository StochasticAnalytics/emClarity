function [ outputVol ] = BH_reScale3d( inputVol, nameOUT, MAG, METHOD, varargin )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if (isnumeric(MAG))
  mag = MAG; 
else
  mag = str2num(MAG);
end

if nargin > 4
  xyzShift = varargin{1};
else
  xyzShift = [0,0,0];
end

pixelSize = 1.0;
if isa(inputVol, 'cell')
  nVols = length(inputVol);
  outPutArray = false;
  writeOut = false;
elseif isnumeric(inputVol)
  nVols = 1;
  outPutArray = true;
  inputVol = {inputVol};
  writeOut = false;
elseif ischar(inputVol)
  [imgPath, imgName, imgExt] = fileparts(inputVol);
  if isempty(imgPath)
    imgPath = '.';
  end
    % Read in the image
  nVols = 1;
  mrcImage = MRCImage(inputVol,0);
  header = getHeader(mrcImage);
  pixelSizeX = header.cellDimensionX / header.nX;
  pixelSizeY = header.cellDimensionY / header.nY;
  if pixelSizeX ~= pixelSizeY
    fprintf('\npixel size in X (%2.2f), Y (%2.2f) inconsistent, leaving unset\n',pixelSizeX,pixelSizeY);
  else
    pixelSize = pixelSizeX;
  end
  inputVol = {getVolume(mrcImage)};
  writeOut = true;
  outPutArray = false;

end
outputVol = cell(nVols,1);

% Assuming all images are the same size
sizeVol = size(inputVol{1});

if length(sizeVol) ~= 3
  error('rescale only set to work on 3d volumes');
end
% % % if mag > 1
% % %   sizeOut = round(sizeVol./mag);
% % % else
% % %   sizeOut = sizeVol;
% % % end
sizeOut = round(sizeVol.*mag);

% Bandlimit first and clear mask to save memory 
freqCutoff  = mag.*0.475;
bandPass = BH_bandpass3d(sizeVol,0,0,1/freqCutoff,METHOD,1);
  
for iVol = 1:nVols
    inputVol{iVol} = real(ifftn(fftn(inputVol{iVol}).*bandPass));                     
end

clear bandPass




[~,~,~,x,y,z] = BH_multi_gridCoordinates(sizeVol,'Cartesian',METHOD,...
                                        {'single',[1,0,0;0,1,0;0,0,1],...
                                        xyzShift','forward',1,mag},0,1,0);

[X,Y,Z,~,~,~] = BH_multi_gridCoordinates(sizeOut,'Cartesian',METHOD,...
                                        {'single',[1,0,0;0,1,0;0,0,1],...
                                        xyzShift','forward',1,mag},0,1,0);


  
for iVol = 1:nVols
  if strcmp(METHOD,'GPU')
    outputVol{iVol} = interpn(x,y,z, inputVol{iVol},X,Y,Z ,'linear',0);
  else
    outputVol{iVol} = interpn(x,y,z, inputVol{iVol},X,Y,Z ,'spline',0);       
    outputVol{iVol}(isnan(outputVol{iVol})) = 0;
  end
  clear inputVol{iVol}                                 
end

clear x y z X Y Z

if ( outPutArray )
  outputVol = outputVol{1};
elseif( writeOut )
  SAVE_IMG(MRCImage(gather(single(outputVol{1}))), nameOUT,pixelSize./mag);
end

end
