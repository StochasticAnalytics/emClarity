function [ peakWgt, sortedList ] = BH_weightAngCheckPeaks(positionList, nPeaks, ...
                                                          symmetry, iSubTomo, tomoName)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

compressByFactor = 2;  
% A larger tolerance would probably catch orientations that should
% converege to the same value earlier on and therby save computation. For
% testing, set a bit more conservatively.
angleTolerance = 5;
                     
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
    distVect = (180/pi)*sqrt(sum((rNxMx-rComb).^2,1));
    if all(distVect < angleTolerance)
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


peakKeep = ( peakWgt~=-9999 );
peakScore = sortedList(1:26:26*nPeaks);
peakScore(~peakKeep) = -9999;
peakWgt(peakKeep) = (peakScore(peakKeep) ./ max(peakScore(peakKeep))).^compressByFactor;
peakWgt(peakKeep) = peakWgt(peakKeep) ./ sum(peakWgt(peakKeep));

fprintf('Weighting %d volumes for st %d from tomo %s, ',nPeaks,iSubTomo,tomoName);
switch nPeaks
  case 2
    fprintf('raw score [%2.2f %2.2f] --> [%2.2f %2.2f]\n',peakScore, peakWgt);
  case 3
    fprintf('raw score [%2.2f %2.2f %2.2f] --> [%2.2f %2.2f %2.2f]\n',peakScore, peakWgt);
  case 4
    fprintf('raw score [%2.2f %2.2f %2.2f %2.2f] --> [%2.2f %2.2f %2.2f %2.2f]\n',peakScore, peakWgt);
  case 5
    fprintf('raw score [%2.2f %2.2f %2.2f %2.2f %2.2f] --> [%2.2f %2.2f %2.2f %2.2f %2.2f]\n',peakScore, peakWgt);
  otherwise
    error('More than five peaks is a bit too spunky.')
end
          
end

