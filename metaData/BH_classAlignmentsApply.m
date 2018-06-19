function [UPDATED_GEOMETRY] = ...
                            BH_classAlignmentsApply( INPUT_GEOMETRY ,...
                                      BEST_ANGLES, SAMPLING, CLASS_VECTOR, iGold, halfSet)
%Apply class alignments to the full set of subTomograms
%   
%
%   Input variables
%
%   ALIGN_GEOM = string with mat file that has the geometry output from 
%                alignment.
%
%   
%   Output variables = none, the clustering geometry updated - eventually this
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
%   TODO
%           1)add check that new origin is not too close to the edge.
%           2)get rid of 'find' and switch from index to logical indexing
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




inputGeometry = INPUT_GEOMETRY;
alignmentGeometry = BEST_ANGLES;
classVector = CLASS_VECTOR(1,:);


% Get the number of tomograms to proc

tomoList = fieldnames(inputGeometry);
nTomograms = length(tomoList);

for iTomo = 1:nTomograms
  
  positionList = inputGeometry.(tomoList{iTomo});
  newAlignment = alignmentGeometry;
    
  if strcmpi(halfSet,'ODD') || strcmpi(halfSet,'EVE')
    includeList =  find(ismember(positionList(:,26), classVector) & ...
                   ismember(positionList(:,7),  iGold));
% % %     removeList  = ~ismember(positionList(:,26), classVector) & ...
% % %                   ismember(positionList(:,7),  iGold);
% % %     positionList(removeList,26) = -9999;
  else
    includeList =  find(ismember(positionList(:,26), classVector));
% % %     removeList  = ~ismember(positionList(:,26), classVector);
% % %     positionList(removeList,26) = -9999;
  end

  for iParticle = includeList'
    % assuming all classes are sequential, only discarded between cycles.
    class =  positionList(iParticle,26);
    pIndex = find(newAlignment(:,2) == class); %%% 1-->2

    if isempty(pIndex)
      fprintf('pIndex is empty for %s pos %d class %d\n', tomoList{iTomo}, iParticle,class);
    else
     classRot = BH_defineMatrix(newAlignment(pIndex,3:5),'Bah','inv');

     oldRot = reshape(positionList(iParticle,17:25),3,3);

     newRot = gather(reshape( oldRot * classRot, 1,9));

     positionList(iParticle,17:25) = newRot;
   
                                                                   
                                       

     inputGeometry.(tomoList{iTomo}) = positionList; 

    end

  end

UPDATED_GEOMETRY = inputGeometry;

end

