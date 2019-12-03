function [gX, gY, gZ, vX, vY, vZ] = EMC_multi_gridMasks(TYPE, SIZE, METHOD, SHIFT, ORIGIN, ...
                                                        flgNormalize, flgHalf, flgIsotrope)
% [gX, gY, gZ, vX, vY, vZ] = EMC_multi_gridMasks(TYPE, SIZE, METHOD, SHIFT, ORIGIN, ...
%                                                flgNormalize, flgHalf, flgIsotrope)
% Compute mesh grids of different size and different coordinate systems.
% Mostly used to compute masks.
%
% TYPE (str):           'Cartesian', 'Spherical', 'Cylindrical' or 'Radial'.
%
% SIZE (vector):        Size of the grid to compute (x, y, z) or (x, y).
%                       Sould be a 2d/3d row vector of integers.
%
% METHOD (str):         Device to use; 'gpu' or 'cpu'.
%
% SHIFT (vector):       (x, y, z) or (x, y) translations to apply; should correspond to SIZE.
%
% ORIGIN (int):         Origin convention
%                       -1: zero frequency first (fft output)
%                       0: real origin (if even nb of pixel, the center is in between 2 pixels)
%                       1: right origin (extra pixel to the left; ceil((N+1)/2))
%                       2: left origin (extra pixel to the right)
%                       emClarity uses the right origin convention (ORIGIN=1).
%
% flgNormalize (bool): 	Normalize the gridVectors between 0 and 1.
% flgHalf (bool):       Compute half of the grids (half of the x axis).
% flgIsotrope (bool):   Stretch the grid values to the smallest dimensions.
%                       NOTE: Useful with TYPE='Radial', for filtering.
%
% EXAMPLE: [X, Y, Z, x, y, z] =  EMC_multi_gridMasks('Radial', [100,90], 'cpu', [0,0], 1, false, false, 1);

%% 
[vX, vY, vZ] = EMC_multi_vectorCoordinates(SIZE, METHOD, SHIFT, ORIGIN, flgNormalize, flgHalf);
if numel(SIZE) == 3
    flg3d = true;
else
    flg3d = false;
end

% This could be more efficient by doing it when creating the vectors directly, with linspace.
if flgIsotrope
    radius = floor((SIZE + 1) / 2);
    radius_min = min(radius);
    vX = vX .* (radius_min / radius(1));
    vY = vY .* (radius_min / radius(2));
    if (flg3d)
        vZ = vZ .* (radius_min / radius(3));
    else
        vZ = nan;
    end
end

if (flg3d)
    switch TYPE
        case 'Cartesian'
            [gX, gY, gZ] = ndgrid(vX, vY, vZ);
        case 'Radial'
            gX = vX'.^2 + vY.^2 + reshape(vZ, 1, 1, []).^2;
            gY = nan; gZ = nan;
        case 'Spherical'
            [X, Y, Z] = ndgrid(vX, vY, vZ);
            gX = sqrt(vX'.^2 + vY.^2 + reshape(vZ, 1, 1, []).^2);
            gY = atan2(Y, X);
            gY(gY < 0) = gY(gY < 0) + 2.*pi;  % set from [-pi,pi] --> [0,2pi]
            gZ = acos(Z ./ gX);  % [0,pi]
        case 'Cylindrical'
            [X, Y, gZ] = ndgrid(vX, vY, vZ);
            gX = sqrt(X.^2 + Y.^2);
            gY = atan2(Y, X);
            gY(gY < 0) = gY(gY < 0) + 2.*pi;  % set from [-pi,pi] --> [0,2pi]
        otherwise
            error("SYSTEM should be  'Cartesian', 'Spherical', 'Cylindrical' or 'Radial', got %s", SYSTEM)
    end
else % 2d
    switch TYPE
        case 'Cartesian'
            [gX, gY] = ndgrid(vX, vY); gZ = nan;
        case 'Radial'
            gX = vX'.^2 + vY.^2;
            gY = nan; gZ = nan;
        case {'Spherical', 'Cylindrical'}
            [X, Y] = ndgrid(vX, vY);
            gX = sqrt(X.^2 + Y.^2);
            gY = atan2(Y, X);
            gZ = nan;
        otherwise
            error("SYSTEM should be  'Cartesian', 'Spherical', 'Cylindrical' or 'Radial', got %s", SYSTEM)
    end
end
end



