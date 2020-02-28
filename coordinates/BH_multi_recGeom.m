function [ recGeom, tiltName, nTomos] = BH_multi_recGeom( reconCoordName )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% could use import data in later versions of matlab (it wasn't working in
% the compiled binaries with < 15a for some reason.)
recFile = fopen(reconCoordName,'r');
tiltName = textscan(recFile,'%s',1) ;
tiltName = tiltName{1}{1};
nTomos = textscan(recFile,'%d',1); nTomos = nTomos{1};
recCoords = textscan(recFile,'%f');
fclose(recFile);

%%% Some sanity checks 
% First line should be the name of the tilt-series the tomo is
% reconstructed from.
[~,tiltNameFromTomo,~] = fileparts(reconCoordName);
tiltStr = strsplit(tiltNameFromTomo,'_');
tiltNameFromTomo = strjoin(tiltStr(1:end-1),'_');
if ~strcmp(tiltNameFromTomo, tiltName)
  error('the tomo base name (%s) does not match the tiltName in the coords file (%s)',tiltNameFromTomo,tiltName);
end

recGeom = zeros(nTomos,6);
for iSt = 1:nTomos
  recGeom(iSt,:) = recCoords{1}(1 + (iSt-1)*6: 6 + (iSt-1)*6);              
end
           
% Note that the x/z shifts (col 5,6) are shifts given to IMOD, which are the opposite of the location of the origin (relative to the center)
% To make it more confusing, since the reconstruction is done in a ref frame rotated about X, the Z shift is flipped so it matches the origin in Z
end

