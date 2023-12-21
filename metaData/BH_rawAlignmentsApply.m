
function [UPDATED_GEOMETRY] = ...
                            BH_rawAlignmentsApply( INPUT_GEOMETRY , ...
                            BEST_ANGLES, SAMPLING, nPeaks, rotConvention, updateWeights, updateClassByBestReferenceScore)
%Apply class alignments to the full set of subTomograms
%   
%
%   Input variables
%
%   ALIGN_GEOM = string with mat file that has the geometry output from 
%                alignment.
%
%
%   BEST_ANGLES =
%
%       1 - iSubTomo
%       2 - Best reference
%       3 - CCC
%       4 - wedgeWeight
%       5:7 - Translational Shifts
%       8:16 - Rotation matrix for the highest CCC
%
%
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
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




inputGeometry = INPUT_GEOMETRY;
alignmentGeometry = BEST_ANGLES;



% Get the number of tomograms to proc

tomoList = fieldnames(inputGeometry);
nTomograms = length(tomoList);

for iTomo = 1:nTomograms
 
    positionList = inputGeometry.(tomoList{iTomo});
    newAlignment = alignmentGeometry.(tomoList{iTomo});

    if ~isempty(newAlignment)
      
 for iPeak = 1:nPeaks   
    includeList = find(positionList(:,26 + 26*(iPeak-1)) ~= -9999)
   
    for iParticle = includeList'
     % assuming all classes are sequential, only discarded between cycles.
     particleIDX =  positionList(iParticle,4 + 26*(iPeak-1));
      
     pIndex = find(newAlignment(:,2) == particleIDX);
    
     
     newAngles = newAlignment(pIndex,[3:5] + 10*(iPeak-1));
    if isempty(pIndex)
          positionList(iParticle,26 + 26*(iPeak-1)) = -9999;
          fprintf('particle instance %d in catch clause rawAlignmentApply set to ignore\n',particleIDX);
          continue
    else
 
     classRot = BH_defineMatrix(newAngles,rotConvention,'inv');
    end
% % % % %      catch
    
      

% % % % %      end
     
  
     oldRot = reshape(positionList(iParticle,[17:25] + 26*(iPeak-1)),3,3);
     newRot = gather(reshape(oldRot * classRot, 1,9));
     
     shifts = gather(newAlignment(pIndex,[8:10] + 10*(iPeak-1)).*SAMPLING);
     positionList(iParticle,[17:25] + 26*(iPeak-1)) = newRot;
     positionList(iParticle,[11:13] + 26*(iPeak-1)) = ...
                                gather(positionList(iParticle,[11:13]+ 26*(iPeak-1))) + ...
                                                                         shifts;
     

     positionList(iParticle,1+26*(iPeak-1)) = gather(newAlignment(pIndex,[6] + 10*(iPeak-1)));
     % I don't think I'm using this column anywhere else but double check
     % classification maybe. FIXME.
     if (updateWeights)
      positionList(iParticle,2+26*(iPeak-1)) = gather(newAlignment(pIndex,[7] + 10*(iPeak-1)));
     end

     if (updateClassByBestReferenceScore)
      positionList(iParticle,26+26*(iPeak-1)) = gather(newAlignment(pIndex,[1] + 10*(iPeak-1)));
     end

    end
 end % loop over peaks
    end
    inputGeometry.(tomoList{iTomo}) = positionList;
end

UPDATED_GEOMETRY = inputGeometry;

end

