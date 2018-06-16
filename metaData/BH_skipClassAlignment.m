function [  ] = BH_skipClassAlignment( PARAMETER_FILE, CYCLE, STAGEofALIGNMENT, flgGold )
%Skip the alignment of class averages, while updating the appropriate
%fields in subTomoMeta
%   In Particular, useful for the first cycle (0) so that averges may be
%   made from classification, and then used to refine the template matching
%   results. This is important for the monolayer, because template matching
%   takes into account the lattice (aka a. priori information) but is (by
%   design) fairly indeterministic when it comes to the angular search.
%   This may lead to gross errors in the translational search for class
%   members that have been placed in the wrong class due to missing wedge
%   bias, which is only weakly accounted for in template matching.s

if (nargin ~= 4)
  error('args = PARAMETER_FILE, CYCLE, STAGEofALIGNMENT, flgGold')
end

flgGold = str2num(flgGold);
if (flgGold)
  sfx = {'ODD','EVE'}
else
  sfx = {'STD','STD'}
end
% Read in the geometry from the selected classification, save as if it were
% the result of class alignment. For now this doesn't mean anything as
% class identity changes cycle to cycle, and is not considered in between.
% This may not be the case in the future, so the "correct" geometry is
% explicity copied here.
cycleNumber = sprintf('cycle%0.3u', str2num(CYCLE));

pBH = BH_parseParameterFile(PARAMETER_FILE);

flgClassify = pBH.('flgClassify');
try
  flgMultiRefAlignment = pBH.('flgMultiRefAlignment');
catch
  flgMultiRefAlignment = 0;
end

load(sprintf('%s.mat', pBH.('subTomoMeta')), 'subTomoMeta');
outputPrefix = sprintf('%s_%s', cycleNumber, pBH.('subTomoMeta'));
%className    = pBH.('Cls_className');
%refName = pBH.('Ref_className');

if strcmpi(STAGEofALIGNMENT, 'ClassAlignment')
 
  if (flgMultiRefAlignment && ~flgClassify)
    subTomoMeta.(cycleNumber).('ClassAlignment') = ...
                                subTomoMeta.(cycleNumber).('Avg_geometry');

	

  elseif (flgMultiRefAlignment && flgClassify)
     subTomoMeta.(cycleNumber).('ClassAlignment') = ...
                                subTomoMeta.(cycleNumber).('ClusterRefGeom');

  else

    try 
     subTomoMeta.(cycleNumber).('ClassAlignment') = ...
                                subTomoMeta.(cycleNumber).('ClusterClsGeom');
    catch
     subTomoMeta.(cycleNumber).('ClassAlignment') = ...
                                subTomoMeta.(cycleNumber).('ClusterRefGeom');
    end

 %   if (flgGold)
 %     features1     = pBH.('Pca_coeffs_odd');
 %     features2     = pBH.('Pca_coeffs_eve');
 %     GEOM1 = sprintf('%s_%d_%d_nClass_%d_ODD',outputPrefix,features1(1,1), ...
 %                                             features1(1,end), className);
 %     GEOM2 = sprintf('%s_%d_%d_nClass_%d_EVE',outputPrefix,features2(1,1), ...
 %                                             features2(1,end), className);  
 %                                           
 %     % if class averages weren't extracted, just references, check and use these
 %    % instead. build in tool "isfield" doesn't work for dynamic names, so use
 %     % eval, which will work if the field exists, otherwise throw an error.
 %     
 %    try
 %    %checkClassAvg = isfield(subTomoMeta, sprintf( '(%s).(%s).(%s)', cycleNumber, 'ClusterResults',GEOM1));
 %     eval(sprintf( 'subTomoMeta.(''%s'').(''%s'').(''%s'')', cycleNumber, 'ClusterResults',GEOM1));
 %    catch
 % 
 %    GEOM1 = sprintf('%s_%d_%d_nClass_%d_ODD',outputPrefix,features1(1,1), ...
 %                                             features1(1,end), refName);
 %     GEOM2 = sprintf('%s_%d_%d_nClass_%d_EVE',outputPrefix,features2(1,1), ...
 %                                            features2(1,end), refName); %%
%
%     end
%    
%  
%  
%    [ GEOM_OUT ] = BH_mergeClassGeometry( subTomoMeta.(cycleNumber).('ClusterResults').(GEOM1), ...
%                                          subTomoMeta.(cycleNumber).('ClusterResults').(GEOM2));%
%
%    else
%      features1     = pBH.('Pca_coeffs_odd');
%      GEOM1 = sprintf('%s_%d_%d_nClass_%d_STD',outputPrefix,features1(1,1), ...
%                                              features1(1,end), className)
%      GEOM_OUT = subTomoMeta.(cycleNumber).('ClusterResults').(GEOM1);
%   end
% 
%    subTomoMeta.(cycleNumber).('ClassAlignment') = GEOM_OUT;
  end
elseif strcmpi(STAGEofALIGNMENT, 'RawAlignment')
  subTomoMeta.(cycleNumber).('RawAlign') = ...
                                subTomoMeta.(cycleNumber).('ClassAlignment');
else
  error(['STAGEofALIGNMENT to skip may be RawAlignment or ClassAlignment,'],...
        ['the former requires the latter to exist.\n']);
end

save(pBH.('subTomoMeta'), 'subTomoMeta');   



end

