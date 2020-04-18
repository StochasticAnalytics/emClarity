function [ GEOMETRY_SPLIT ] = BH_fscSplit( GEOMETRY, splitOnTomos, nPeaks )
%Randomly divide a data set into two halves for fsc-gold analysis.
%  
%
%   Input variables:
%
%   GEOMETRY = A structure with tomogram names as the field names, and geometry
%              information in a 20 column array.
%              Additionally, a field called 'source_path' has a value with the
%              absolute path to the location of the tomograms.
%
%              The input is a string 'Geometry_templatematching.mat' for
%              example, and it is expected that the structure is saved as the
%              variable named geometry.
%
%   OUTPUT_NAME = Something like 'Geometry_fscGold-cycle-000.mat'
%
%   Output variables:
%
%   None - mat file with 1/2 in column 7 written to disk.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Goals & limitations:
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   TODO:
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isstruct(GEOMETRY)
  subTomoMeta = GEOMETRY;
  geometry = subTomoMeta.('cycle000').geometry;
  returnStruct = true;
else
  try
    load(GEOMETRY, 'subTomoMeta');
    geometry = subTomoMeta.('cycle000').geometry;
    returnStruct = false;
    GEOMETRY_SPLIT = '';
  catch error
    error('failed to load geometry')
  end
end

% Get the number of tomograms to process.
tomoList = fieldnames(geometry);
nTomograms = length(tomoList);

floorCeil = 0;
for iTomo = 1:nTomograms
  
  
  % Load in the geometry for the tomogram, and get number of subTomos.
  positionList = geometry.(tomoList{iTomo});


  nSubTomos = size(positionList,1);
  % Generate random half sets for fsc
  if (floorCeil)
    halfValue = floor(nSubTomos./2);
    floorCeil = 0;
  else
    halfValue = ceil(nSubTomos./2);
    floorCeil = 1;
  end
  
  
% % %   % Add randomized positions for ML intialization
% % %   if (nPeaks > 1)
% % %     positionList = repmat(positionList,1,nPeaks);
% % %     for iPeak = 1:nPeaks-1
% % %       
% % %       for iRand = 1:size(positionList,1)
% % %         positionList(iRand,[17:25] + 26*iPeak) = ...
% % %            reshape(BH_defineMatrix('rand','Bah','inv')*reshape(positionList(iRand,[17:25]),3,3),1,9);
% % %       end
% % %       
% % %     end  
% % %   end
  
  if (splitOnTomos)
    positionList(:,7:26:26*nPeaks) = 1 + floorCeil;
  else
    halfData = datasample([1:nSubTomos], halfValue, 'Replace', false);
    % Set all values column 7 = 1
    positionList(:,7:26:26*nPeaks) = 1;
    % Replace half data column 7 with 2
    positionList(halfData,7:26:26*nPeaks) = 2;
    %numel(find(positionList(:,7)==1))
    % Update geometry
  end
  

  geometry.(tomoList{iTomo}) = positionList;
  
  
end

subTomoMeta.('cycle000').geometry=geometry;

if (returnStruct)
  GEOMETRY_SPLIT = subTomoMeta;
else
  save(GEOMETRY, 'subTomoMeta');
end

end % end of fscSplit

