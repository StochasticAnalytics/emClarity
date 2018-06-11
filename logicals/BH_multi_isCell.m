function [ idxOut ] = BH_multi_isCell( cellInput )
%Return a vector of non empty cell contents
%   Detailed explanation goes here

idxOut = [];
for iCell = 1:length(cellInput)
  if ~isempty(cellInput{iCell})
    idxOut = [idxOut; iCell];
  end
end

end

