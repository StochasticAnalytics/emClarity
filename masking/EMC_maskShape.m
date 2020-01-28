function [MASK] = EMC_maskShape(SHAPE, SIZE, RADIUS, METHOD, OPTION)
%
% MASK = EMC_maskShape(SHAPE, SIZE, RADIUS, METHOD, OPTION)
% Create a real space mask of a given shape.
%
% Input:
%   SHAPE (str | numeric):          'sphere', 'cylinder' or 'rectangle'.
%
%   SIZE (int vector):              Size (in pixel) of the mask to compute; [x, y, z] or [x, y].
%                                   NOTE: [1, N] or [N, 1] is not allowed.
%
%   RADIUS (int vector):            Radius (in pixel) of the shape to compute; [x, y, z] or [x, y].
%                                   Should correspond to SIZE, and be greater or equal to 2.
%
%   METHOD (str):                 	Device to use; 'gpu' or 'cpu'.
%
%   OPTION (cell | struct):       	Optional parameters.
%                                   If cell: {field, value; ...}, note the ';' between parameters.
%                                   NOTE: Can be empty.
%                                   NOTE: Unknown fields will raise an error.
%
%     -> 'shift' (vector):          [x, y, z] or [x, y] translations to apply (should correspond to SIZE).
%                                   NOTE: NaNs or Inf are not accepted.
%                                   default = no shifts
%
%  	  -> 'origin' (int):            Origin convention - Center of the shape.
%                                   0, 1 or 2; see EMC_multi_gridVectors for more details.
%                                   NOTE: origin=-1 is not allowed as this function is only for real
%                                         space masks.
%                                   NOTE: origin=0 is currently not allowed for SHAPE = 'rectangle' and
%                                         'cylinder'.
%                                   default = 1
%
%  	  -> 'taper' (bool|float|matrix):
%                                   Gaussian taper to apply to the boolean mask.
%                                   If bool: apply or not the default taper.
%                                   If float: the size of the taper, in percentage (between 0 and 1) of
%                                   pixels, relative to the smallest axis. The minimum size of the taper
%                                   is automatically set to 8 pixels.
%                                   If matrix: specify your own kernel (odd dimensions) to convolve
%                                   with the boolean mask.
%                                   NOTE: The variance of the gaussian is automatically adjusted to keep
%                                         the gaussian at ~0.5 at the middle of the taper.
%                                   default = 0.04
%
%  	  -> 'sym' (int):           	Restrict the mask to the first asymmetric unit.
%                                   Should correspond to the central symmetry (positive int).
%                                   If 1, return the full mask.
%                                   default = 1
%
%	  -> 'precision' (str):       	Precision of the output MASK; 'single' or 'double'.
%                                   default = 'single'
%
% Output:
%   MASK (num array):           	2d/3d mask.
%
% Note:
%   - If taper, the function forces the edges of the mask to be blurred.
%     The pixels at the edges of the array will be 0.
%   - I [TF] have tried doing the convolution in Fourier space, but in most situations it is
%     considerably slower...
%
% Example:
%   [MASK] = EMC_maskShape('cylinder', [128,128], [30,30], 'cpu', {})
%   [MASK] = EMC_maskShape('cylinder', [128,128], [30,30], 'cpu', {'shift',[0,10]})
%

%% 
[SIZE, OPTION, flg, ndim] = checkIN(SIZE, RADIUS, METHOD, OPTION);

cutoffLow = 0.001;  % everything below this value is set to 0.

if flg.taper && ~flg.ownTaper
    % Compute the size of the kernel in pixel.
    smallestDimension = min(SIZE);
    kernelSize = round((smallestDimension * OPTION.taper - 1) / 2) * 2 + 1;  % closest odd int
    if kernelSize < 9; kernelSize = 9; end  % at least 9 pixels, which is equivalent to a 6 pixel taper.

    % Keep the gaussian at ~0.5 at the middle of the taper.
    middle = ceil(kernelSize / 2) / 2;
    sigma = sqrt(-1 * middle^2 / (2 * log(0.5)));
    kernel = EMC_gaussianKernel(zeros(1,ndim) + kernelSize, sigma, {'precision', OPTION.precision; 'method', METHOD});
end

%% Compute the desired binary shape
if strcmpi(SHAPE, 'sphere')
    % If any blurring is meant to be applied to this mask, the radius of the ellipsoid needs to be
    % increased by half of the kernel size to preserve the desired radius.
    RADIUS = RADIUS + ceil(kernelSize/2);
    
    % Define the center of the sphere/circle using the specified origin and shifts.
    [vX, vY, vZ] = EMC_coordVectors(SIZE, METHOD, {'origin', OPTION.origin; ...
                                                   'shift', OPTION.shift; ...
                                                   'precision', OPTION.precision});

    % Compute a radial cartesian grid and apply the equation of the sphere/ellipsoid.
    % The surface of the ellipsoid is at 1. At this point, the binary ellipsoid mask is computed.
    if (flg.is3d)
        MASK = (vX'./RADIUS(1)).^2 + (vY./RADIUS(2)).^2 + reshape(vZ./RADIUS(3),1,1,[]).^2 <= 1;
    else
        MASK = (vX'./RADIUS(1)).^2 + (vY./RADIUS(2)).^2 <= 1;
    end
    
    if strcmpi(OPTION.precision, 'single')
        MASK = single(MASK);
    else
        MASK = double(MASK);
    end
    
elseif strcmpi(SHAPE, 'cylinder')    
    % Adjust the radius for the same reason as explain with sphere/ellipsoids.
    RADIUS(1:2) = RADIUS(1:2) + ceil(kernelSize/2);

    % To compute a (3d) cylinder, this algorithm first compute a 2d sphere/ellipsoid and then
    % broadcast it along the Z (depth) axis.
    if (flg.is3d)
        % First compute the 2d ellipsoid with the x and y shifts.
        [vX, vY, ~] = EMC_coordVectors(SIZE(1:2), METHOD, {'origin', OPTION.origin; ...
                                                           'precision', OPTION.precision;
                                                           'shift', OPTION.shift(1:2)});

        % Broadcaste in Z (the ellipsoid is invariant in Z <=> cylinder).
        RADIUS(3) = RADIUS(3) + floor(kernelSize/2);  % adjust Z for blurring.
        if (flg.gpu)
            vZ = zeros([1, 1, RADIUS(3) .*2 + 1], OPTION.precision, 'gpuArray');
        else
            vZ = zeros([1, 1, RADIUS(3) .*2 + 1], OPTION.precision);
        end
      	MASK = (vX' ./ RADIUS(1)).^2 + (vY ./ RADIUS(2)).^2 + vZ <= 1;  % broadcast

        if strcmpi(OPTION.precision, 'single')
            MASK = single(MASK);
        else
            MASK = double(MASK);
        end
    
        % Resize the cylinder to the desired SIZE in Z, taking into account the z shift.
        limits = EMC_limits(size(MASK), SIZE, {'origin', OPTION.origin; 'shift', [0,0,OPTION.shift(3)]});
        MASK = EMC_resize(MASK, limits, {'taper', false});

    else  % 2d Cylinder; this block is equivalent to SHAPE='sphere'.
        [vX, vY] = EMC_coordVectors(SIZE, METHOD, {'origin', OPTION.origin; ...
                                                   'precision', OPTION.precision; ...
                                                   'shift', OPTION.shift});
        MASK = (vX' ./ RADIUS(1)).^2 + (vY ./ RADIUS(2)).^2 <= 1;
        if strcmpi(OPTION.precision, 'single')
            MASK = single(MASK);
        else
            MASK = double(MASK);
        end
    end

elseif strcmpi(SHAPE, 'rectangle')
    % Adjust the radius for the convolution.
    RADIUS = RADIUS + floor(kernelSize/2);

    if (flg.gpu)
        MASK = ones(RADIUS .*2 + 1, OPTION.precision, 'gpuArray');
    else
        MASK = ones(RADIUS .*2 + 1, OPTION.precision);
    end
    limits = EMC_limits(size(MASK), SIZE, {'origin', OPTION.origin; 'shift', OPTION.shift});
    MASK = EMC_resize(MASK, limits, {'taper', false});
else
    error("SHAPE should be 'sphere', 'cylinder' or 'rectangle', got %s", SHAPE)
end

%% Restrict the mask to the first symmetry pair.
if (flg.sym)
    [~, angles, ~] = EMC_maskGrids('cylindrical', SIZE, METHOD, {'shift', OPTION.shift; ...
                                                                 'origin', OPTION.origin; ...
                                                                 'precision', OPTION.precision});
    sectorMax = 2 * pi / OPTION.sym * 1.025;  % small overlap
    angles = (angles > (2*pi-sectorMax/2) | angles < sectorMax/2);
    MASK = MASK .* angles;
end

%% Finally, apply the blurring kernel to the binary mask.
% Use EMC_resize to tape the current binary MASK edges with zeros. Convolved with the blurring
% kernel, it creates a nice roll off to zero at the edges of the MASK.
if (flg.taper)
    % Force taper at the edges.
    MASK = EMC_resize(MASK, zeros(1,ndim*2), {'force_taper', true; 'taper', zeros(1, ceil(kernelSize/2))});
    MASK = convn(MASK, kernel, 'same');  % this is the most expensive part
    MASK = MASK ./ max(MASK(:));  % rounding errors; max=1
    MASK(MASK <= cutoffLow) = 0;  % predictable end of the pass
end

end  % EMC_shapeMask


function [SIZE, OPTION, flg, ndim] = checkIN(SIZE, RADIUS, METHOD, OPTION)

[flg.is3d, SIZE, ndim] = EMC_is3d(SIZE);

if ~isnumeric(RADIUS) || ~isvector(RADIUS) || any(rem(SIZE,1)) || any(SIZE < 1)
    error('EMC_maskShape:RADIUS', 'RADIUS should be a vector of %d integers greater than 1', ndim)
end

if strcmpi(METHOD, 'gpu')
    flg.gpu = true;
elseif strcmpi(METHOD, 'cpu')
    flg.gpu = false;
else
    error('EMC_maskShape:METHOD', "METHOD should be 'gpu' or 'cpu', got %s", METHOD)
end

OPTION = EMC_getOption(OPTION, {'shift', 'origin', 'taper', 'sym', 'precision'}, false);

% shift
if isfield(OPTION, 'shift')
    if ~isnumeric(OPTION.shift) || ~isvector(OPTION.shift)
        error('EMC_maskShape:shift', ...
              'shift should be a vector of float|int, got %s', class(OPTION.shift))
    elseif any(isnan(OPTION.shift)) || any(isinf(OPTION.shift))
        error('EMC_maskShape:shift', ...
              'shift should not contain NaNs or Inf, got %s', mat2str(OPTION.shift, 2))
   	elseif numel(OPTION.shift) ~= ndim
        error('EMC_maskShape:shift', ...
              'For a %dd SIZE, shift should be a vector of %d float|int, got %s', ...
              ndim, ndim, mat2str(OPTION.shift, 2))
    end
else
    OPTION.shift = zeros(1, ndim);  % default
end


% origin
if isfield(OPTION, 'origin')
    if ~isscalar(OPTION.origin) || ~(OPTION.origin == 0 || OPTION.origin == 1 || OPTION.origin == 2)
        % EMC_resize (used with 'rectangle' and 3d 'cylinders') will raise an error if origin=0
        error('EMC_maskShape:origin', "origin should be 0, 1, or 2, got %.04f", OPTION.origin)
    end
else
    OPTION.origin = 1;  % default
end

% taper
flg.ownTaper = false;
if isfield(OPTION, 'taper')
    if islogical(OPTION.taper)
        if OPTION.taper
            flg.taper = true;
            OPTION.taper = 0.04;  % default
        else
            flg.taper = false;
            OPTION.taper = 0;
        end
    elseif isfloat(OPTION.taper) || isinteger(OPTION.taper)
        if ~(0 <= OPTION.taper <= 1)
            error('taper should be between 0 and 1, got %.04f', OPTION.taper)
        end
        if OPTION.taper == 0 || isnan(OPTION.taper)
            flg.taper = false;
        else
            flg.taper = true;
        end
    elseif isnumeric(OPTION.taper)
        taperSize = size(OPTION.taper);
        if ~all(taperSize(1) == taperSize)
            error('kernel should have the same number of elements in each dimension.')
        elseif ~mod(taperSize(1), 2)
            error('kernel size should be odd')
        end
        flg.taper = true;
        flg.ownTaper = true;
    else
        error('taper should be a boolean or a positive float, got %s', class(OPTION.taper))
    end
else
     OPTION.taper = 0.04;  % default
     flg.taper = true;
end

% sym
if isfield(OPTION, 'sym')
    if isnumeric(OPTION.sym) && isscalar(OPTION.sym) && ~rem(OPTION.sym, 1) && OPTION.sym > 0
        if OPTION.sym == 1
            flg.sym = false;
        else
            flg.sym = true;
        end
    else
        error('sym should be a positive integer')
    end
else
    flg.sym = false;  % default
end

% precision
if isfield(OPTION, 'precision')
    if ~(strcmpi(OPTION.precision, 'single') || strcmpi(OPTION.precision, 'double'))
        error("precision should be 'single' or 'double', got %s", OPTION.precision)
    end
else
    OPTION.precision = 'single';
end

end  % checkIN
