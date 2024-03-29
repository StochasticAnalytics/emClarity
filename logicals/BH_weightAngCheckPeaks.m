function [ peakWgt, sortedList ] = BH_weightAngCheckPeaks(positionList, nPeaks, ...
                                                          score_sigma, iSubTomo, tomoName,...
                                                          track_stats)
                                                          
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% A larger tolerance would probably catch orientations that should
% converege to the same value earlier on and therby save computation. For
% testing, set a bit more conservatively.
angleTolerance = 5;
compressByFactor = 2;  

% Check that the class id is not set to ignore
includedPeaks = reshape(find(positionList(26:26:26*nPeaks) ~= -9999),[],1);
excludedPeaks = reshape(find(positionList(26:26:26*nPeaks) == -9999),[],1);
peakWgt = zeros(nPeaks,1);
sortedList = zeros(size(positionList));

if numel(includedPeaks) == 1
  % No need to wgt or look for angular sep
  incIndex = [1+26*(includedPeaks(1)-1):26+26*(includedPeaks(1)-1)];
  sortedList(1:26) = positionList(incIndex);
  sortedList(27:end) = positionList(~ismember(1:length(positionList),incIndex));
  peakWgt(1) = 1;
  peakWgt(2:end) = -9999;
  fprintf('Weighting %d volumes for st %d from tomo %s\n',1,iSubTomo,tomoName);
  return
else
  peakMax = [ includedPeaks, [positionList(1+26*(includedPeaks-1))]' ];
  peakMax = sortrows(peakMax, -2);
  
  pIDX = 1;
  for iPeak = 1:length(includedPeaks)
    incIndex = 1+26*(peakMax(iPeak,1)-1):26+26*(peakMax(iPeak,1)-1);
    sortedList(1+26*(pIDX-1):26+26*(pIDX-1)) = positionList(incIndex);
    pIDX = pIDX +1;
  end
  for iPeak = 1:length(excludedPeaks)
    % The order of the excluded peaks with respect to each other doesn't matter.
    incIndex = 1+26*(excludedPeaks(iPeak)-1):26+26*(excludedPeaks(iPeak)-1);
    sortedList(1+26*(pIDX-1):26+26*(pIDX-1)) = positionList(incIndex);
    pIDX = pIDX +1;
  end
  
  % re-determine sorted included positions
  includedPeaks = reshape(find(sortedList(26:26:26*nPeaks) ~= -9999),[],1);
  excludedPeaks = reshape(find(sortedList(26:26:26*nPeaks) == -9999),[],1);
  peakWgt(includedPeaks) = 1; 
  peakWgt(excludedPeaks) = -9999;
  
  nIncluded = numel(includedPeaks);
end

% List to loop over to check for unique angles
combinations = [nIncluded:-1:2,1];
nCHk = nchoosek(combinations,2);
 
for iMax = nIncluded:-1:2
  % Rotation matrix for the peak in question

  rNxMx = reshape(sortedList(17+26*(iMax-1):25+26*(iMax-1)),3,3);
  for iComb = find(nCHk(:,1) == iMax)
    rComb = reshape(sortedList(17+26*(nCHk(iComb,2)-1):25+26*(nCHk(iComb,2)-1)),3,3);
%     distVect = (180/pi)*sqrt(sum((rNxMx-rComb).^2,1));
    % From http://www.boris-belousov.net/2016/12/01/quat-dist/#using-rotation-matrices
    R = rNxMx * transpose(rComb);
    angDist = acosd((trace(R)-1)./2);


    if all(angDist < angleTolerance)
      % this is a duplicate rotation, set its class to -9999 to ignore.
      sortedList(26+26*(iMax-1)) = -9999;
      peakWgt(iMax) = -9999;
      break
    end
      
  end
end
  
% % Sort with peak # in col one descending on the CCC
% rankedScores = sortrows([includedPeaks, positionList(1 + 26.*(includedPeaks-1))'],-2);

% The columns of a rotation matrix tell you where the basis vectors are
% mapped to by the rotation, and since these are orthonormal, the distance
% between any two pairs of basis vectors is the angle between them. So
% sqrt(sum((R1-R2).^2,1)*(180/pi) gives the distance in degrees for each of
% the three basis vectors under the two given rotations.

% FIXME change this to a proper angular difference.

minScore = 1e-1;

peakScore = sortedList(1:26:26*nPeaks)' ./ score_sigma;
peakKeep = ( peakWgt~=-9999 & abs(peakScore) > minScore);
peakScore(~peakKeep) = -9999;
peakWgt(~peakKeep) = -9999;
if (track_stats)
  compressByFactor = sqrt(max(peakScore(peakKeep)));
end
peakWgt(peakKeep) = (peakScore(peakKeep) ./ max(peakScore(peakKeep))).^compressByFactor;
peakWgt(peakKeep) = peakWgt(peakKeep) ./ sum(peakWgt(peakKeep));
peakWgt = real(peakWgt);
peakWgt(peakWgt < minScore & peakKeep) = -9999;
if any(isnan(peakWgt))
  error('Found a nan in the re-weighted peaks')
end

if mod(iSubTomo,25)
  
  fprintf('Weighting %d volumes with cf %3.3f for st %d from tomo %s, ',nPeaks,compressByFactor,iSubTomo,tomoName);
  fprintf('raw score [');
  for i = 1:nPeaks
    fprintf('%1.3f ', peakScore(i));
  end
  fprintf('] --> [');
  for i = 1:nPeaks
    fprintf('%1.3f ', peakWgt(i));
  end
  fprintf(']\n');

  
end

end



