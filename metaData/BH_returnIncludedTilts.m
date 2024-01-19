function [tiltList,nTilts] = BH_returnIncludedTilts(mapBackGeom)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

STACK_LIST_tmp = fieldnames(mapBackGeom);
STACK_LIST_tmp = STACK_LIST_tmp(~ismember(STACK_LIST_tmp,{'tomoName','viewGroups'}));
tiltList = cell(length(STACK_LIST_tmp),1);

nTilts = 0;
for iStack = 1:length(STACK_LIST_tmp)
  if mapBackGeom.(STACK_LIST_tmp{iStack}).nTomos
    nTilts = nTilts + 1;
    tiltList{nTilts} = STACK_LIST_tmp{iStack};
  end
end

tiltList = tiltList(1:nTilts);
end

