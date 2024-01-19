function [UPDATED_GEOMETRY] = ...
  BH_refAlignmentsApply( INPUT_GEOMETRY,BEST_ANGLES,...
  SAMPLING, ...
  REF_VECTOR,REF_GROUP,flgGold)
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
refVector = REF_VECTOR(1,:);
refGroup  = REF_GROUP;


% Get the number of tomograms to proc

tomoList = fieldnames(inputGeometry);
nTomograms = length(tomoList);

for iTomo = 1:nTomograms
  
  positionList = inputGeometry.(tomoList{iTomo});
  newAlignment = alignmentGeometry;
  
  
  includeList =  find(ismember(positionList(:,26), refVector) & ...
    ismember(positionList(:,7),  flgGold));
  removeList  = ~ismember(positionList(:,26), refVector) & ...
    ismember(positionList(:,7),  flgGold);
  positionList(removeList,26) = -9999;
  
  
  for iParticle = includeList'
    % assuming all classes are sequential, only discarded between cycles.
    class =  positionList(iParticle,26);
    pIndex = find(newAlignment(:,2) == class); %%% 1-->2
    
    if ~isempty(pIndex)
      
      classRot = BH_defineMatrix(newAlignment(pIndex,3:5),'Bah','inv');
      
      oldRot = reshape(positionList(iParticle,17:25),3,3);
      newRot = gather(reshape( oldRot * classRot, 1,9));
      positionList(iParticle,17:25) = newRot;
      
      shifts = gather((oldRot*classRot )*(1.*newAlignment(pIndex,8:10).*SAMPLING)')';
      %         shifts = gather(positionList(iParticle,14:16) - newAlignment(pIndex,8:10).*SAMPLING);
      % replace the classIDX with the groupIDX
      positionList(iParticle,11:13) = gather(positionList(iParticle,11:13)) + ...
        shifts;
      find(refVector == class);
      positionList(iParticle,26) = refGroup(find(refVector == class));
      
      
      
      inputGeometry.(tomoList{iTomo}) = positionList;
    end
  end
  
end

UPDATED_GEOMETRY = inputGeometry;

end

