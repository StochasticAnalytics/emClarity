function [gX, gY, gZ, vX, vY, vZ] = EMC_coordGrids(SYSTEM, SIZE, METHOD, OPTION)
%
% [gX, gY, gZ, vX, vY, vZ] = EMC_coordGrids(SYSTEM, SIZE, METHOD, OPTION)
% Compute ndgrids of different size and different coordinate systems.
%
% Input:
%   SYSTEM (str):            	'cartesian', 'spherical', 'cylindrical' or 'radial'.
%
%   SIZE (vector):            	Size (in pixel) of the grids to compute; [x, y, z] or [x, y].
%                               NOTE: [1, N], [N, 1] or [1,1] are not allowed.
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
%                               If SYSTEM = 'radial', gX is the radial grid and gX = gY = NaN.
%
%   vX, vY, vZ (num vectors):   Coordinate vectors. See EMC_coordVectors for more details.
%
% Example:
%   [X, Y, ~, x, y, ~] =  EMC_coordGrids('radial', [100,90], 'cpu', {'origin',-1});
%   [X, Y, Z, x, y, z] =  EMC_coordGrids('spherical', [100,90,100], 'cpu', {});
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
%           v.1.1.1 Rename (masking/EMC_maskGrids to coord/EMC_coordGrids), switch to new
%                   error identifier convention and unittest (TF, 30Jan2020).
%

%% checkIN
[vX, vY, vZ] = EMC_coordVectors(SIZE, METHOD, OPTION);  % all optional inputs are accepted.

% None of the vectors should be NaN <=> SIZE should describe a 2d or 3d array.
if any(SIZE(1:2) == 1)
    error('EMC:SIZE', 'SIZE should describe a 2d or 3d array, got %s', mat2str(SIZE))
elseif ~(ischar(SYSTEM) || isstring(SYSTEM))
    error('EMC:SYSTEM', 'SYSTEM should be a string|char, got %s', class(SYSTEM))
end

if numel(SIZE) == 3 && SIZE(3) ~= 1
    if strcmpi(SYSTEM, 'cartesian')
        [gX, gY, gZ] = ndgrid(vX, vY, vZ);
    elseif strcmpi(SYSTEM, 'radial') || strcmpi(SYSTEM, 'radius')
        gX = sqrt(vX'.^2 + vY.^2 + reshape(vZ, 1, 1, []).^2);
        gY = nan; gZ = nan;
    elseif strcmpi(SYSTEM, 'spherical')
        [X, Y, Z] = ndgrid(vX, vY, vZ);
        gX = sqrt(vX'.^2 + vY.^2 + reshape(vZ, 1, 1, []).^2);
        gY = atan2(Y, X);
        gY(gY < 0) = gY(gY < 0) + 2.*pi;  % set from [-pi,pi] to [0,2pi]
        gZ = acos(Z ./ gX);  % [0,pi]
    elseif strcmpi(SYSTEM, 'cylindrical')
        [X, Y, gZ] = ndgrid(vX, vY, vZ);
        gX = sqrt(vX'.^2 + vY.^2 + reshape(vZ, 1, 1, []).^2);
        gY = atan2(Y, X);
        gY(gY < 0) = gY(gY < 0) + 2.*pi;  % set from [-pi,pi] to [0,2pi]
    else
        error('EMC:SYSTEM', "SYSTEM should be  'cartesian', 'spherical', 'cylindrical' or 'radial'")
    end
else % 2d
    if strcmpi(SYSTEM, 'cartesian')
        [gX, gY] = ndgrid(vX, vY); gZ = nan;
    elseif strcmpi(SYSTEM, 'radial')
        gX = sqrt(vX'.^2 + vY.^2);
        gY = nan; gZ = nan;
    elseif strcmpi(SYSTEM, 'spherical') || strcmpi(SYSTEM, 'cylindrical')
        [X, Y] = ndgrid(vX, vY);
        gX = sqrt(vX'.^2 + vY.^2);
        gY = atan2(Y, X);
        gY(gY < 0) = gY(gY < 0) + 2.*pi;  % set from [-pi,pi] to [0,2pi]
        gZ = nan;
    else
       	error('EMC:SYSTEM', "SYSTEM should be  'cartesian', 'spherical', 'cylindrical' or 'radial'")
    end
end

end  % EMC_coordGrids
