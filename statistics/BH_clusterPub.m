function [ ] = BH_clusterPub(PARAMETER_FILE, CYCLE)
%Kmeans based classification
%   
%
%   Input Variables:
%
%   GEOMETRY =
%
%   COEFF_MAT = matfile with previous pca decomposition
%
%   nCLUSTERS = vector with number of clusters to try.
%
%   COEFFS = cell with 1x2 vectors giving ranges of coeffs to try
%            e.g. {[2,40], [7,40]}
%
%   kDIST = distance measure to use. I have observed some improved seperation
%           for my data using 'cosine' rather than the default.
%
%           'sqeuclidean', 'cityblock', 'cosine', 'correlation'
%
%   kREP = number of replicates for each
%
%   Output Variables:
%
%   None - writes out an updated geometry for each combination of nClUSTERS and
%          COEFFS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Goals & limitations:
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   TODO
%
%   Change parpool to 48 prior to testing on archer. Also take a look into the
%   available GPU accelerated K means. Check the memory used for coeff and
%   whether or not this is limiting. It should not be.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (nargin ~= 2)
  error('PARAMETER_FILE, CYCLE')
end

startTime =  clock;
CYCLE = str2num(CYCLE);
cycleNumber = sprintf('cycle%0.3u', CYCLE);

pBH = BH_parseParameterFile(PARAMETER_FILE);

flgClassify = pBH.('flgClassify');
%%% For general release, I've disabled class average alignment and
%%% multi-reference alignment, so set the default to OFF. If either of
%%% these features are re-introduced, this will need to be reverted.
if ( flgClassify ); flgClassify = -1 ; end
if flgClassify < 0
  flgGold = 0;
else
  flgGold = 1;
end

try
  nPeaks = pBH.('nPeaks');
catch
  nPeaks = 1;
end

featureVector = cell(2,1);
if flgGold
  featureVector{1,1} = pBH.('Pca_coeffs_odd');
  featureVector{2,1} = pBH.('Pca_coeffs_eve');
else
  featureVector{1,1} = pBH.('Pca_coeffs')
  featureVector{1,1}
end


clusterVector= pBH.('Pca_clusters');

try
  kDIST        = pBH.('Pca_distMeasure');
catch
  kDIST = 'sqeuclidean';
end
try
  kREP         = pBH.('Pca_nReplicates');
catch
  kREP = 128;
end

nCores       = BH_multi_parallelWorkers(pBH.('nCpuCores'));

% try
%   relativeScale = pBH.('Pca_relativeScale')
% catch
%   relativeScale= ones(size(clusterVector,1),1);
% end

try
  flgFlattenEigs = pBH.('Pca_flattenEigs');
catch
  flgFlattenEigs=0;
end

try
  load(sprintf('%s.mat', pBH.('subTomoMeta')), 'subTomoMeta');
  geometry_clean = subTomoMeta.(cycleNumber).Avg_geometry;
  geometry       = subTomoMeta.(cycleNumber).Avg_geometry; 
  masterTM = subTomoMeta; clear subTomoMeta
catch 
  error('failed to load geometry')
end

for iGold = 1:1+flgGold
  if (flgGold)
    if iGold == 1
      halfSet = 'ODD';
      randSet = 1;
    else
      halfSet = 'EVE';
      randSet = 2;
    end
  else
    halfSet = 'STD';
    randSet = [1,2];
  end 
  
  coeffMatrix  = sprintf('%s_%s_%s_pcaFull.mat',cycleNumber,pBH.('subTomoMeta'),halfSet);
  outputPrefix = sprintf('%s_%s', cycleNumber, pBH.('subTomoMeta'));
  % Get the number of tomograms to process.
  tomoList = fieldnames(geometry_clean);
  nTomograms = length(tomoList);


  flgRefineKmeans = false;
  kAlgorithm = 'kMeans';
 % kAlgorithm = 'neuralNetwork'

  kDist = sprintf('%s', kDIST);

  switch kDist
    case 'sqeuclidean'
      kDistMeasure = 'sqeuclidean' 
    case 'cityblock'
      kDistMeasure = 'cityblock'
    case 'cosine'
      kDistMeasure = 'cosine'
 
    case 'correlation'
      kDistMeasure = 'correlation'
    case 'ward'
      kDistMeasure = 'ward'
      kAlgorithm = 'HAC'
    case 'neural'
      kDistMeasure = 'neural'
      kAlgorithm = 'neuralNetwork'
    otherwise
      kDistMeasure = 'sqeuclidean'
      fprintf(['\nDefaulting to sqeuclidean b/c %s was not recognized'] ...
              , kDist);
  end

  kReplicates =  kREP;

 %kDistMeasure = 'euclidean'

  try
    oldPca = load(coeffMatrix);
    coeffsUNTRIMMED = oldPca.coeffs
    idxList = oldPca.idxList;
    if nPeaks > 1
      peakList = oldPca.idxList;
    else
      peakList = [];
    end
    
    numParticles = oldPca.nTOTAL;
    clear oldPca;
  catch 
    error('trouble loading the previous pcs mat file.')
  end
  
  try
    EMC_parpool(nCores);
  catch
    delete(gcp('nocreate'));
    pause(3)
    EMC_parpool(nCores);
  end

  %%% experimental part of pcaMS
  nScaleSpace = size(featureVector{iGold},1)
  featureVector{1}
  nFeatures = zeros(1,nScaleSpace)
  
%   if length(relativeScale) ~= nScaleSpace
%     error('relativeScale has %d elements for %d scaleSpaces', ...
%           length(relativeScale), nScaleSpace);
%   end
%   
  if isa(coeffsUNTRIMMED, 'cell')
    [nI,nJ] = size(coeffsUNTRIMMED{1});
    for iScale = 1:nScaleSpace
      nFeatures(iScale) = nnz(featureVector{iGold}(iScale,:));
    end
    coeffMat = zeros(sum(nFeatures), nJ, 'single');
    nAdded = 0
    for iScale = 1:nScaleSpace
      fV = featureVector{iGold}(iScale,:);
      fV = sort(fV(fV~=0))
     
      coeffMat(1+nAdded:nAdded+nFeatures(iScale),:) = ...
                            double(coeffsUNTRIMMED{iScale}(ismember(1:nI,fV),:));
       
    % normalizing the variance of the rows gives equal weight to each eigenvector which
    % is not reasonable as they are by their nature scaled by the amount of
    % variance explained across a given dimension.
% % %       coeffMat(1+nAdded:nAdded+nFeatures(iScale),:) = ...
% % %         coeffMat(1+nAdded:nAdded+nFeatures(iScale),:) ./ ...
% % %         repmat(rms(coeffMat(1+nAdded:nAdded+nFeatures(iScale),:),2),1,nJ).*iScale;
    
      if (flgFlattenEigs)
        coeffMat(1+nAdded:nAdded+nFeatures(iScale),:) = ...
          coeffMat(1+nAdded:nAdded+nFeatures(iScale),:) ./ ...
          repmat(rms(coeffMat(1+nAdded:nAdded+nFeatures(iScale),:),2),1,nJ);
      end
    % Instead, maintain option to weight the features from different scale
    % spaces relative to each other.
%       coeffMat(1+nAdded:nAdded+nFeatures(iScale),:) = ...
%         coeffMat(1+nAdded:nAdded+nFeatures(iScale),:) .* relativeScale(iScale);    
      
      nAdded = nAdded + nFeatures(iScale)
    end
    %coeffsUNTRIMMED = coeffMat; clear coeffMat
  end
  
  for iCluster = 1:length(clusterVector)
    nClusters = clusterVector(iCluster);
    

      if strcmpi(kAlgorithm, 'kMeans')
        [class, classCenters, sumd, D] = kmeans(coeffMat', nClusters, ...
                                    'replicates', kReplicates, ...
                                    'Distance', kDistMeasure, ...
                                    'MaxIter', 50000, ... % Default was 100
                                    'Options', statset('UseParallel', 1) );
                                  
      elseif strcmpi(kAlgorithm, 'kMedoids')
        [class, classCenters, sumd, D] = kmedoids(coeffMat', nClusters, ...
                                    'replicates', kReplicates, ...
                                    'Distance', kDistMeasure, ...
                                     'Options', statset('UseParallel', 1, ...
                                                        'MaxIter', 50000) );
                                                      
      elseif strcmpi(kAlgorithm, 'HAC')
        [class] = clusterdata(coeffMat', ...
                              'maxclust', nClusters, ...
                              'linkage', 'ward', ...
                              'distance', 'euclidean', ... 
                              'savememory', 'off');
                            sumd = 0;

      elseif strcmpi(kAlgorithm, 'neuralNetwork')
        net = selforgmap([1 nClusters]);
        [net, tr] = train(net, coeffMat);
        y = net(coeffMat)
        class = vec2ind(y)

      else
        error('kAlgorithm must be kMeans, or kMedoids, not %s', kAlgorithm);
      end
      
     
      fprintf('\n\nSum of dist to all centroids for each replicate.\n\n')
      fprintf('%g\n',sumd)
      totSum1 = sum(sumd);
      totStd1 = std(sumd);
      fprintf('Total kmeans dist = %g\n', totSum1)
      fprintf('Total kmeans std  = %g\n', totStd1)

      if (flgRefineKmeans)
        % Using the postions found, refine the original estimates

        kMin = min(classCenters,[],1);
        kMax = max(classCenters,[],1);
        kRange = kMax - kMin;

        % Maximum percentages of the range to search around

        k01 = 0.0001 .* kRange;
        k05 = 0.001 .* kRange;
        k10 = 0.01 .* kRange;
        k25 = 0.1 .* kRange;

        seedMatrix = rand([size(classCenters),512],'single');

        seedMatrix(:,:,1:128)  = repmat(k01,size(classCenters,1),1,128) .* ... 
                                                            seedMatrix(:,:,1:128);

        seedMatrix(:,:,129:256) = repmat(k05,size(classCenters,1),1,128) .* ... 
                                                          seedMatrix(:,:,129:256);

        seedMatrix(:,:,257:384) = repmat(k10,size(classCenters,1),1,128) .* ... 
                                                          seedMatrix(:,:,257:384);

        seedMatrix(:,:,385:512) = repmat(k25,size(classCenters,1),1,128) .* ... 
                                                          seedMatrix(:,:,385:512);


%         [class, classCenters, sumd] = kmeans(coeffs(features, :)', nClusters, ...
%                                     'Start', seedMatrix, ...
%                                     'Distance', kDistMeasure, ...
%                                     'MaxIter', 20000, ... % Default was 100
%                                     'Options', statset('UseParallel', 1) );
        if strcmpi(kAlgorithm, 'kMeans')
          [class, classCenters, sumd,D] = kmeans(coeffMat', nClusters, ...
                                      'replicates', kReplicates, ...
                                      'Distance', kDistMeasure, ...
                                      'MaxIter', 50000, ... % Default was 100
                                      'Options', statset('UseParallel', 1) );

        elseif strcmpi(kAlgorithm, 'kMedoids')
          [class, classCenters, sumd,D] = kmedoids(coeffMat', nClusters, ...
                                      'replicates', kReplicates, ...
                                      'Distance', kDistMeasure, ...
                                       'Options', statset('UseParallel', 1, ...
                                                          'MaxIter', 50000) );
        else
          error('kAlgorithm must be kMeans, or kMedoids, not %s', kAlgorithm);
        end
        
        fprintf('\n\nSum of dist to all centroids for each replicate.\n\n')
        fprintf('%g\n',sumd);
        totSum2 = sum(sumd);
        totStd2 = std(sumd);
        fprintf('Total kmeans dist = %g\n', totSum2);
        fprintf('Total kmeans std  = %g\n', totStd2);
        fprintf('percent change in mean %g\n', ...
                                (totSum2 - totSum1)./max(totSum1,totSum2) .* 100);
        fprintf('percent change in std  %g\n', ...
                                (totStd2 - totStd1)./max(totStd1,totStd2) .* 100);


      end

      % This leaves each cluster untouched, but changes the cluster label
      % such that the cluster lablelled 1 is also the most populated
      % cluster.
      classCount = zeros(nClusters,1);
      for i = 1:nClusters
        % counts of class numbers
        classCount(i) = sum(class(:) == i);
      end
      % list of classIDX with highest count first
      [~, ndx] = sort(classCount, 'descend');
      newClass = class;
      for i = 1:nClusters
        newClass(class == ndx(i)) = i;
      end

      fileOUT = fopen(sprintf('%s_%s_ClassIDX.txt',pBH.('subTomoMeta'),cycleNumber), 'a');
      fprintf(fileOUT, '\n\n%s, %s, %s\n','position','idx','count'); 
      for iClass = 1:nClusters
        fprintf(fileOUT, '%d, %d\n',iClass,sum(newClass == iClass));
      end
      fclose(fileOUT);

      class = newClass;
      clear newClass ndx;


      if length(idxList) ~= length(class)
        error('idxList ~= class')
      end

      % This isn't great, and maybe my brain is just tired.
      save(sprintf('clusterTrouble_%d.mat',iCluster), 'idxList', 'class','classCenters','D');

      for iTomo = 1:nTomograms
        positionList = geometry_clean.(tomoList{iTomo});
        includedClass = ( positionList(:,26) ~= -9999 & ismember(positionList(:,7),randSet));
        positionList(:,8) = includedClass;
        particleIDX = positionList( includedClass, 4);
  
         % Returns the lowest index where this is true
        [~, lIndClass] = ismember(particleIDX, idxList);
        [~, lIndPart]  = ismember(particleIDX, positionList(:,4));

        % for trouble shooting
        try
      
        if (nPeaks > 1) 
          for thisIDX = 1:length(lIndClass)
            for iPeak = 0:nPeaks-1
              positionList(lIndPart(thisIDX), 26 + 26*iPeak) = class(lIndClass(thisIDX)+iPeak);
%               fprintf('iTomo %d iSubtomo %d iPeak %d Class %d\n',iTomo,lIndPart(thisIDX),iPeak+1,class(lIndClass(thisIDX)+iPeak));
            end
          end
          geometry.(tomoList{iTomo}) = positionList;
          fprintf('Size iTomo %d %d\n',size(positionList));
        else
          positionList(lIndPart, 26) = class(lIndClass);
          geometry.(tomoList{iTomo}) = positionList;
          
        end
        catch
          save('ClusterLine391Err.mat')
          error('Caught error, saving workspace for evaluation.\n')
        end
      end
      fout = sprintf('%s_%d_%d_nClass_%d_%s', outputPrefix, featureVector{iGold}(1,1), ...
                                               featureVector{iGold}(1,end), nClusters, halfSet);

      % Save a copy of the geometry in the subTomoMeta and also save the name for
      % easy reference in a text file.
      masterTM.(cycleNumber).('ClusterResults').(fout) = geometry;


   
  end % loop over cluster size

  subTomoMeta = masterTM;

  save(pBH.('subTomoMeta'), 'subTomoMeta');

  %save(sprintf('%s_pca.mat',OUTPUT_PREFIX), 'nTOTAL','U', 'S', 'V', 'coeffs')
  fprintf('Total execution time on set %s: %f seconds\n', halfSet,etime(clock, startTime));
  delete(gcp('nocreate'));
end % end of Gold loop
end % end of cluster function

