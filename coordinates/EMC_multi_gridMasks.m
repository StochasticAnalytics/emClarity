function [gX, gY, gZ, vX, vY, vZ] = EMC_multi_gridMasks(TYPE, SIZE, METHOD, OPTION)
% [gX, gY, gZ, vX, vY, vZ] = EMC_multi_gridMasks(TYPE, SIZE, METHOD, OPTIONAL)
%
% Compute ndgrids of different size and different coordinate systems.
%
% TYPE (str):                   'cartesian', 'spherical', 'cylindrical' or 'radial'.
%
% SIZE (vector):                Size of the grid to compute (x, y, z) or (x, y).
%                               Sould be a 2d/3d row vector of integers.
%
% METHOD (str):                 Device to use; 'gpu' or 'cpu'.
%
% OPTION (cell | struct):       Optional parameters.
%                               If cell: {field,value ; ...}, note the ';' between parameters.
%                               NOTE: Can be empty.
%                               NOTE: Unknown fields will raise an error.
%
%   -> 'shift' (vector):        (x, y, z) or (x, y) translations to apply; should correspond to SIZE.
%                               default = no shifts
%
%   -> 'origin' (int):          Origin convention
%                               -1: zero frequency first (fft output)
%                               0: real origin (if even nb of pixel, the center is in between 2 pixels)
%                               1: right origin (extra pixel to the left; ceil((N+1)/2))
%                               2: left origin (extra pixel to the right)
%                               default = 1
%
%   -> 'normalize' (bool):      Normalize the gridVectors between 0 and 1.
%                               default = false
%
%   -> 'half' (bool):           Compute half of the grids (half of the x axis).
%                               default = false
%
%   -> 'isotrope' (bool):       Stretch the grid values to the smallest dimensions.
%                               default = false
%
%   -> 'precision' (str):       Precision of the grids and vectors; 'single' or 'double'
%                               default = 'single'
%
%--------
% RETURN:                       gX, gY, gZ are the gridMasks.
%                               vX, vY, vZ are the vectorCoordinates.
%
%--------
% EXAMPLE: [X, Y, Z, x, y, z] =  EMC_multi_gridMasks('radial', [100,90], 'cpu', {'origin',-1});
%          [X, Y, Z, x, y, z] =  EMC_multi_gridMasks('spherical', [100,90,100], 'cpu', {});

%% checkIN
[vX, vY, vZ] = EMC_multi_vectorCoordinates(SIZE, METHOD, OPTION);  % all optional inputs are accepted.

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
        gY(gY < 0) = gY(gY < 0) + 2.*pi;  % set from [-pi,pi] --> [0,2pi]
        gZ = acos(Z ./ gX);  % [0,pi]
    elseif strcmpi(TYPE, 'cylindrical')
        [X, Y, gZ] = ndgrid(vX, vY, vZ);
        gX = sqrt(X.^2 + Y.^2);
        gY = atan2(Y, X);
        gY(gY < 0) = gY(gY < 0) + 2.*pi;  % set from [-pi,pi] --> [0,2pi]
    else
        error("SYSTEM should be  'cartesian', 'spherical', 'cylindrical' or 'radial', got %s", TYPE)
    end
else % 2d
    if strcmpi(TYPE, 'cartesian')
        [gX, gY] = ndgrid(vX, vY); gZ = nan;
    elseif strcmpi(TYPE, 'radial')
        gX = sqrt(vX'.^2 + vY.^2);
        gY = nan; gZ = nan;
    elseif strcmpi(TYPE, 'spherical') || strcmpi(TYPE, 'cylindrical')
        [X, Y] = ndgrid(vX, vY);
        gX = sqrt(X.^2 + Y.^2);
        gY = atan2(Y, X);
        gZ = nan;
    else
        error("SYSTEM should be  'cartesian', 'spherical', 'cylindrical' or 'radial', got %s", TYPE)
    end
end
end



