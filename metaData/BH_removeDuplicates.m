function [  ] = BH_removeDuplicates( PARAMETER_FILE, CYCLE )
%Remove duplicates where some particles have drifted to the same postion.
%   Detailed explanation goes here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Goals & Limitations:
%
%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   TODO
%
%    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% radius +/- for duplicate tolerance. Even though you may not expect any
% particles w/in a radius equal to that of the particle, this would require
% a much more expensive calculation. It is likely that any particle that
% comes within an "event horizon" will end up at a nearly identical origin 
% volume; a smaller radius is probably okay

CYCLE = EMC_str2double(CYCLE); 
cycleNumber = sprintf('cycle%0.3u', CYCLE);
emc = BH_parseParameterFile(PARAMETER_FILE);

dupSampling = ceil(10e-10 / emc.('PIXEL_SIZE'));

pixelSize = emc.('PIXEL_SIZE').*dupSampling.*10^10;
if emc.('SuperResolution')
  pixelSize = pixelSize * 2;
end
latticeRadius = emc.('particleRadius');

dupRadius = max(1,floor(0.2*min(latticeRadius)/pixelSize));
dupTolerance = (2.*dupRadius)+1;

dupMask(dupTolerance, dupTolerance, dupTolerance) = gpuArray(single(0));
dupMask = dupMask + 1;

if (nargin ~= 2)
  error('args = PARAMETER_FILE, CYCLE')
end


% Backup the current geometry
  system(sprintf('cp %s.mat preDupRemoval_%s.mat',emc.('subTomoMeta'),emc.('subTomoMeta')));
  load(sprintf('%s.mat', emc.('subTomoMeta')), 'subTomoMeta');
  geometry = subTomoMeta.(cycleNumber).RawAlign;
  masterTM = subTomoMeta; clear subTomoMeta


% Get the number of tomograms to process.
tomoList = fieldnames(geometry);
nTomograms = length(tomoList);
nRemoved = 0;
nTotal = 0;
for iTomo = 1:nTomograms
  
  % Load in the geometry for the tomogram, and get number of subTomos.
  positionList = geometry.(tomoList{iTomo});
  includeList = find(positionList(:,26) ~= -9999);
 
  nTotal = nTotal + length(includeList);
  
  tomoNumber = masterTM.mapBackGeometry.tomoName.(tomoList{iTomo}).tomoNumber;
  tiltName = masterTM.mapBackGeometry.tomoName.(tomoList{iTomo}).tiltName;
  tomoName = sprintf('%s_%d',tiltName,tomoNumber);
  
  recGeom = masterTM.reconGeometry.(tomoName);
%   iHeader = getHeader(MRCImage(tomoName));
  
  clear postionMatrix positionIDX
  % leave IDX in main memory because it is just for reference.
%   sx = floor(iHeader.nX ./ dupSampling);
%   sy = floor(iHeader.nY ./ dupSampling);
%   sz = floor(iHeader.nZ ./ dupSampling);
   sx = floor(recGeom(1,1)./dupSampling);
   sy = floor(recGeom(1,2)./dupSampling);
   sz = floor(recGeom(1,3)./dupSampling);
  positionMatrix = zeros([sx,sy,sz],'single','gpuArray');
  positionIDX = zeros([sx,sy,sz],'single');
  
  % Make a volume with ones in the position of the centers of the tomos.
  for iSubTomo = includeList'
   
    subTomoOrigin = round(positionList(iSubTomo,11:13)./dupSampling);
    if any(subTomoOrigin < 1 + dupRadius) || any([sx,sy,sz] < subTomoOrigin + dupRadius)
      positionList(iSubTomo,26) = -9999;
    else
    positionMatrix(subTomoOrigin(1),subTomoOrigin(2),subTomoOrigin(3)) = 1;
    positionIDX(subTomoOrigin(1),subTomoOrigin(2),subTomoOrigin(3)) = ...
                                                    positionList(iSubTomo, 4);
    end
       
  end % loop building position matrix
  
  % Convolve positionmatrix with duplicate mask. The numbers in the
  % resulting matrix will correspond to the number of duplicates within the
  % specified radius.
  
  overlapMatrix = convn(positionMatrix, dupMask, 'same');
  
  % Positions inbetween particle origins will also be non-zero, so restrict
  % search to be particle origins that are within radius.
  duplicateList = find( (positionMatrix) & (overlapMatrix > 1) );
  length(duplicateList)
  for iDup = duplicateList'
    
    % get the positions within radius, get corresponding particle ids, find
    % highest CCC, set class id to -9999 for others, also remove those ids
    % from duplicate list by setting value in duplicate list to -1, which
    % is ignored.
    
   
    [i,j,k] = ind2sub([sx,sy,sz], iDup);
    % Check that the duplicate hasn't already been evaluated
    try
    if positionList((positionList(:,4) == positionIDX(iDup)),3) ~= -9999
      dupWindow = positionIDX( i - dupRadius : i + dupRadius, ...
                               j - dupRadius : j + dupRadius, ...
                               k - dupRadius : k + dupRadius ) ;
      % From window, select only real particle ids                      
      idxList = dupWindow(dupWindow ~= 0);
      % Logical translating particle ids to postions in geometry file
      posList = ismember(positionList(:,4), idxList);
      % Replace ones in logical with CCC from previous raw Alignment
      cccList = max(positionList(:,1:26:26*emc.nPeaks),[],2).*posList;
      
      [~ , maxCCCcoord] = max(cccList);
      % set highest CCC to zero, so the remaining are all inferior
      % duplicates, set these to ignore class, and also keep a record at
      % column 3
      posList(maxCCCcoord) = 0;
      nRemoved = nRemoved + sum(posList);
      positionList(posList,26:26:26*emc.nPeaks) = -9999;
      positionList(posList,3:26:26*emc.nPeaks)  = -9999;
      
    end
    catch
      iDup
    end
    
    
    
  end % loop over duplicates

  
  
  
  
  % For now just use a convolution, later change this to fourier based
  % which will be much faster.
  
  
  
  
  geometry.(tomoList{iTomo}) = positionList;
end % end loop over tomorams

fprintf('%d of %d particles removed\n', nRemoved, nTotal);

subTomoMeta = masterTM;
subTomoMeta.(cycleNumber).RawAlign = geometry;
save(emc.('subTomoMeta'), 'subTomoMeta');

end
 
