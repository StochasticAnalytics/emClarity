function [ GEOM_OUT ] = BH_mergeClassGeometry( GEOM1 , GEOM2)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


try
  % Get the number of tomograms to process.
  tomoList1 = fieldnames(GEOM1);
  geometry1 = GEOM1; clear GEOM1
  tomoList2 = fieldnames(GEOM2);
  geometry2 = GEOM2; clear GEOM2
catch error
  error('Could not access the fieldnames in the struct geometry.')
end

nTomograms = length(tomoList1);

% Loop through counting particles available. 
for iTomo = 1:nTomograms
  % Read in the geometry for each tomogram, update column eight and count the
  % number included.
   positionList1 = geometry1.(tomoList1{iTomo});
   positionList2 = geometry2.(tomoList2{iTomo});

   positionList1((positionList2(:,7)==2),:) = positionList2((positionList2(:,7)==2),:);
   
   geometry1.(tomoList1{iTomo}) = positionList1;
   
   
   
end

GEOM_OUT = geometry1;
end

