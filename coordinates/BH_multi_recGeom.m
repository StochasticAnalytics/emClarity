function [ recGeom, tiltName, nTomos, tilt_geometry ] = BH_multi_recGeom( reconCoordName, mapBackIter )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% File format:
% 1. Name of the tilt-series the tomo is reconstructed from
% 2. Number of tomograms
% 3. For each tomogram, numbered sequentially:
%   1. NX
%   2. NY start
%   3. NY end
%   4. NZ
%   5. X shift
%   6. Z shift
recFile = importdata(reconCoordName);
tiltName = recFile.textdata{1};
nTomos = recFile.data(1);
recCoords = recFile.data(2:end);

tilt_geometry_name = sprintf('fixedStacks/ctf/%s_ali%d_ctf.tlt', tiltName, mapBackIter+1);
try
  tilt_geometry  = load(tilt_geometry_name);
catch
  error('Could not load the tilt geometry file: %s', tilt_geometry_name);
end

% Check that we have a multiple of 6 entries. This should have been caught in recSCript2.sh
if mod(numel(recCoords),6) ~= 0
  error('The number of entries in the reconCoord file is not a multiple of 6');
end


%%% Some sanity checks
% First line should be the name of the tilt-series the tomo is reconstructed from.
[~,tiltNameFromTomo,~] = fileparts(reconCoordName);
tiltStr = strsplit(tiltNameFromTomo,'_');
tiltNameFromTomo = strjoin(tiltStr(1:end-1),'_');
if ~strcmp(tiltNameFromTomo, tiltName)
  error('the tomo base name (%s) does not match the tiltName in the coords file (%s)',tiltNameFromTomo,tiltName);
end
% The tilt_geometry should have rows that are a multiple of 26 (> 26 means more than one orientation per peak)
if mod(size(tilt_geometry,2),26) ~= 0 && mod(size(tilt_geometry,2),23) ~= 0
  tilt_geometry
  error('The tilt geometry file does not have a multiple of 26 entries/row');
end

% This includes all possible tomos from a given tilt-series when defined. 
% Tomos may be ignored when cleaning template matching results, or later if set to be ignored
% in geometryAnalysis or if there are zero sub-tomos left.
recGeom = cell(nTomos,1);
for iTomo = 1:nTomos
  read_in_Coords = recCoords(1 + (iTomo-1)*6: 6 + (iTomo-1)*6);
  tomoName = sprintf('%s_%d',tiltName, iTomo);

  % Check to make sure no out of bounds conditions were created in X Y
  % when the user created the model or point file
  % TODO: only checking Y b/c that results in a crash. Checking X would probably make sense too.
  if read_in_Coords(2) < -75
    error(sprintf('Out of bounds condition for %s yMin at %f, please change recon.txt recon.coords', tomoName, read_in_Coords(2)))
  elseif read_in_Coords(2) < 1
    % If not too extreme, just clamp it to 1
    read_in_Coords(2) = 1;
  end
  yMax = tilt_geometry(1,21);
  if read_in_Coords(3) > yMax + 75
    error(sprintf('Out of bounds condition for %s ymin at %f ymax at %f, please change recon.txt recon.coords', tomoName, read_in_Coords(3), yMax))
  elseif read_in_Coords(3) > yMax
    read_in_Coords(3) = yMax;
  end


  tomoCoords = struct();
  tomoCoords.('y_i') = (read_in_Coords(2));
  tomoCoords.('y_f') = (read_in_Coords(3));
  tomoCoords.('NX') = (read_in_Coords(1));
  tomoCoords.('NY') = (read_in_Coords(3) - read_in_Coords( 2) + 1);
  tomoCoords.('NZ') = (read_in_Coords(4));
  tomoCoords.('dX_specimen_to_tomo') = -1*read_in_Coords(5);
  tomoCoords.('dY_specimen_to_tomo') = ...
  (emc_get_origin_index(tomoCoords.('NY')) ...
    + read_in_Coords(2)) ... % origin of the tomogram in the full tilt projection
    - emc_get_origin_index(tilt_geometry(1,21));
  tomoCoords.('dZ_specimen_to_tomo') = read_in_Coords( 6);

  % Check that NX, NY, NZ are all positive
  if tomoCoords.('NX') <= 0
    error('NX is not positive for %s', tomoName);
  end
  if tomoCoords.('NY') <= 0
    error('NY is not positive for %s', tomoName);
  end
  if tomoCoords.('NZ') <= 0
    error('NZ is not positive for %s', tomoName);
  end
  

  recGeom{iTomo} = tomoCoords;  
end



% Note that the x/z shifts (col 5,6) are shifts given to IMOD, which are the opposite of the location of the origin (relative to the center)
% To make it more confusing, since the reconstruction is done in a ref frame rotated about X, the Z shift is flipped so it matches the origin in Z
end

