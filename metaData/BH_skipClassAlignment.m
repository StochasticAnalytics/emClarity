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

flgGold = EMC_str2double(flgGold);
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
cycleNumber = sprintf('cycle%0.3u', EMC_str2double(CYCLE));

emc = BH_parseParameterFile(PARAMETER_FILE);

flgClassify = emc.('flgClassify');
try
  flgMultiRefAlignment = emc.('flgMultiRefAlignment');
catch
  flgMultiRefAlignment = 0;
end

load(sprintf('%s.mat', emc.('subTomoMeta')), 'subTomoMeta');
outputPrefix = sprintf('%s_%s', cycleNumber, emc.('subTomoMeta'));


if strcmpi(STAGEofALIGNMENT, 'RawAlignment')
 
  if (flgMultiRefAlignment && ~flgClassify)
    subTomoMeta.(cycleNumber).('RawAlign') = ...
                                subTomoMeta.(cycleNumber).('Avg_geometry');

  elseif (flgMultiRefAlignment && flgClassify)
     subTomoMeta.(cycleNumber).('RawAlign') = ...
                                subTomoMeta.(cycleNumber).('ClusterClsGeom');
  else

    try 
     subTomoMeta.(cycleNumber).('RawAlign') = ...
                                subTomoMeta.(cycleNumber).('ClusterClsGeom');
    catch
     subTomoMeta.(cycleNumber).('RawAlign') = ...
                                subTomoMeta.(cycleNumber).('ClusterRefGeom');
    end

  end
else
  error(['STAGEofALIGNMENT to skip may be RawAlignment'],...
        ['the former requires the latter to exist.\n']);
end

save(emc.('subTomoMeta'), 'subTomoMeta');   


end

