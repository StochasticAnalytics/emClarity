function [ GEOMETRY_UPDATED, nTOTAL, nSUBSET ] = ...
  BH_randomSubset( GEOMETRY, PCAorSNR, nPARTICLES, halfSET )
%Count the particles to be used in PCA, optionally making a random subset.
%
%
%   Input variables:
%
%   GEOMETRY = geometry structure.
%
%   nPARTICLES = -1, count all non-ignored particles (class -9999)
%                float, randomly select this many particles for the
%                decomposition, denote by updating the flag in column 8 to be 1
%                or 0.
%
%   Output variables:
%
%   GEOMETRY_UPDATED = updated geometry
%
%   nTOTAL = total number to be used for classification
%
%   nSUBSET = total number to be used for decomposition
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Goals & limitations:
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   TODO:
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(PCAorSNR, 'pca')
  % randomized subsetfor determining eigenvectors
  flgSNR = 0;
  flgPCA = 1;
elseif strcmp(PCAorSNR, 'snr')
  % randomly split into 10 sets for snr(nSubtomo) estimation
  flgSNR = 1;
  flgPCA = 0;
else
  error('\nPCAorSNR must be pca or %s\n','snr');
end

if isstruct(GEOMETRY)
  try
    % Get the number of tomograms to process.
    tomoList = fieldnames(GEOMETRY);
    geometry = GEOMETRY; clear GEOMETRY
  catch
    error('Could not access the fieldnames in the struct geometry.')
  end
else
  error('GEOMETRY is not a structure')
end

if ~isnumeric(nPARTICLES)
  error('nPARTICLES must be -1 or a float, not %s', nPARTICLES)
end

nTomograms = length(tomoList);
nIncluded = 0;
includedIDX = [];
% Loop through counting particles available.
for iTomo = 1:nTomograms
  % Read in the geometry for each tomogram, update column eight and count the
  % number included.
  positionList = geometry.(tomoList{iTomo});
  % Count the max number available, use is member for STD alignment
  includedPositions  = ( positionList(:,26) ~= -9999 & ismember(positionList(:,7),halfSET) );
  nIncluded = nIncluded + sum(includedPositions);
  % Save a record of unique particle idx for randomization
  includedIDX = [includedIDX; positionList(includedPositions,4)];
  % Update the position lists if all available are to be used.
  if ( nPARTICLES == -1  && flgPCA)
    positionList(includedPositions, 8) = 1;
    positionList(~includedPositions, 8) = 0;
    geometry.(tomoList{iTomo}) = positionList;
    nSUBSET = nIncluded;
  end
  
end

nTOTAL = nIncluded

if ( nPARTICLES > 0 || flgSNR )
  
  if ( nIncluded >= nPARTICLES )
    fprintf('%d  particles are being used out of %d available\n',...
      nPARTICLES, nIncluded);
  else
    error('%d particles requested, but only %d are available', ...
      nPARTICLES,nIncluded)
  end
  
  rng('shuffle');
  
  if ( flgSNR )
    if nTOTAL < 1500
      error('For now assuming atlease 1500 total particles available.');
    end
    
    
    snrIDX = zeros(numel(includedIDX),2,'uint16');
    snrIDX(:,1) = uint16(includedIDX);
    % 4 replicates at 4 16 36 64 100 144 particles each
    nSubset = [3:3:15].^2
    nSUBSET = nTOTAL;
    for iSubset = 1:length(nSubset)
      for iReplicate = 1:5
        classIDX = (iSubset-1).*5 + iReplicate;
        if iSubset > 1
          nClass = nSubset(iSubset) - nSubset(iSubset-1);
        else
          nClass = nSubset(1);
        end
        [pID, pIDX] = datasample(includedIDX,nClass,'Replace',false);
        % Assign subset id to those selected
        snrIDX(ismember(snrIDX(:,1),pID),2) = classIDX;
        remainingSet = true(size(includedIDX));
        remainingSet(pIDX) = 0;
        includedIDX = includedIDX(remainingSet);
      end
    end
    
    
    
    for iTomo = 1:nTomograms
      
      positionList = geometry.(tomoList{iTomo});
      iSnrSet = ismember(snrIDX(:,1),positionList(:,4));
      
      iSnrSet = snrIDX(iSnrSet,:);
      
      for iParticle = 1:size(iSnrSet,1)
        iParticleIDX = find(positionList(:,4) == iSnrSet(iParticle,1),1,'first');
        positionList(iParticleIDX,10) = iSnrSet(iParticle,2);
      end
      
      geometry.(tomoList{iTomo}) = positionList;
      
    end
    
  else
    
    randomSubset = datasample(includedIDX, nPARTICLES, 'Replace', false);
    nSUBSET = numel(randomSubset);
    randomSubset = reshape(randomSubset, nSUBSET, 1);
    
    for iTomo = 1:nTomograms
      
      positionList = geometry.(tomoList{iTomo});
      positionList(:,8) = single(ismember(positionList(:,4), randomSubset));
      geometry.(tomoList{iTomo}) = positionList;
      
    end
  end
  
end

GEOMETRY_UPDATED = geometry;

end % end of pcaRandomize

