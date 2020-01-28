function [gX, gY, gZ, vX, vY, vZ] = EMC_maskGrids(TYPE, SIZE, METHOD, OPTION)
%
% [gX, gY, gZ, vX, vY, vZ] = EMC_maskGrids(TYPE, SIZE, METHOD, OPTION)
% Compute ndgrids of different size and different coordinate systems.
%
% Input:
%   TYPE (str):             	'cartesian', 'spherical', 'cylindrical' or 'radial'.
%
%   SIZE (vector):            	Size (in pixel) of the grids to compute; [x, y, z] or [x, y].
%                               NOTE: [1, N] or [N, 1] is not allowed.
%
%   METHOD (str):             	Device to use; 'gpu' or 'cpu'.
%
%   OPTION (cell|struct):       Optional parameters.
%                               If cell: {field,value ; ...}, note the ';' between parameters.
%                               NOTE: Can be empty.
%                               NOTE: Unknown fields will raise an error.
%
%     -> 'shift' (vector):     	[x, y, z] or [x, y] translations to apply (should correspond to SIZE).
%                               Shifts are not allowed with half=true or origin=-1.
%                               NOTE: NaNs or Inf are not accepted.
%                               default = no shifts
%
%     -> 'origin' (int):      	Origin convention
%                               default = 1
%
%     -> 'normalize' (bool):	Normalize the vectors between -0.5 and 0.5.
%                               default = false
%
%     -> 'isotrope' (bool):    	Stretch the vector values to the smallest dimensions.
%                               default = false
%
%     -> 'half' (bool):        	Compute 'half' of the vX vector.
%                               default = false
%
%     -> 'precision' (str):   	Precision of the vectors; 'single' or 'double'.
%                               default = 'single'
%
% Output:
%   gX, gY, gZ (num arrays):    X, Y and Z mask grids.
%                               If 2d, gZ = NaN.
%                               If TYPE = 'radial', gX is the radial grid and gX = gY = NaN.
%
%   vX, vY, vZ (num vectors):   Coordinate vectors. See EMC_coordVectors for more details.
%
% Example:
%   [X, Y, ~, x, y, ~] =  EMC_maskGrids('radial', [100,90], 'cpu', {'origin',-1});
%   [X, Y, Z, x, y, z] =  EMC_maskGrids('spherical', [100,90,100], 'cpu', {});
%
% Other EMC-files required:
%   EMC_coordVectors
%
% See also EMC_coordVectors, EMC_is3d
%

% Created:  18Jan2020, R2019a
% Version:  v.1.0   New SIZE convention (see EMC_is3d).
%                   Now follow MATLAB convention for scalars and vectors (TF, 23Jan2020)
%           v.1.1   Rename (EMC_multi_gridMasks to EMC_maskGrids) and unittest (TF, 24Jan2020).
%

%% checkIN
[vX, vY, vZ] = EMC_coordVectors(SIZE, METHOD, OPTION);  % all optional inputs are accepted.

% None of the vectors should be NaN <=> SIZE should describe a 2d or 3d array.
if any(SIZE(1:2) == 1)
    error('EMC_maskGrids:SIZE', 'SIZE should describe a 2d or 3d array, got %s', mat2str(SIZE))
elseif ~(ischar(TYPE) || ~isstring(TYPE))
    error('EMC_maskGrids:TYPE', 'TYPE should be a string|char, got %s', class(TYPE))
end

if numel(SIZE) == 3 && SIZE(3) ~= 1
    if strcmpi(TYPE, 'cartesian')
        [gX, gY, gZ] = ndgrid(vX, vY, vZ);
    elseif strcmpi(TYPE, 'radial') || strcmpi(TYPE, 'radius')
        gX = sqrt(vX'.^2 + vY.^2 + reshape(vZ, 1, 1, []).^2);
        gY = nan; gZ = nan;
    elseif strcmpi(TYPE, 'spherical')
        [X, Y, Z] = ndgrid(vX, vY, vZ);
        gX = sqrt(vX'.^2 + vY.^2 + reshape(vZ, 1, 1, []).^2);
        gY = atan2(Y, X);
        gY(gY < 0) = gY(gY < 0) + 2.*pi;  % set from [-pi,pi] to [0,2pi]
        gZ = acos(Z ./ gX);  % [0,pi]
    elseif strcmpi(TYPE, 'cylindrical')
        [X, Y, gZ] = ndgrid(vX, vY, vZ);
        gX = sqrt(vX'.^2 + vY.^2 + reshape(vZ, 1, 1, []).^2);
        gY = atan2(Y, X);
        gY(gY < 0) = gY(gY < 0) + 2.*pi;  % set from [-pi,pi] to [0,2pi]
    else
        if ~(isstring(TYPE) || ischar(TYPE))
            error('EMC_maskGrids:TYPE', ...
                  "SYSTEM should be  'cartesian', 'spherical', 'cylindrical' or 'radial', got %s", TYPE)
        else
            error('EMC_maskGrids:TYPE', ...
                  "SYSTEM should be  'cartesian', 'spherical', 'cylindrical' or 'radial', got %s", class(TYPE))
        end
    end
else % 2d
    if strcmpi(TYPE, 'cartesian')
        [gX, gY] = ndgrid(vX, vY); gZ = nan;
    elseif strcmpi(TYPE, 'radial')
        gX = sqrt(vX'.^2 + vY.^2);
        gY = nan; gZ = nan;
    elseif strcmpi(TYPE, 'spherical') || strcmpi(TYPE, 'cylindrical')
        [X, Y] = ndgrid(vX, vY);
        gX = sqrt(vX'.^2 + vY.^2);
        gY = atan2(Y, X);
        gY(gY < 0) = gY(gY < 0) + 2.*pi;  % set from [-pi,pi] to [0,2pi]
        gZ = nan;
    else
       	if ~(isstring(TYPE) || ischar(TYPE))
            error('EMC_maskGrids:TYPE', ...
                  "SYSTEM should be  'cartesian', 'spherical', 'cylindrical' or 'radial', got %s", TYPE)
        else
            error('EMC_maskGrids:TYPE', ...
                  "SYSTEM should be  'cartesian', 'spherical', 'cylindrical' or 'radial', got %s", class(TYPE))
        end
    end
end

end  % EMC_maskGrids
