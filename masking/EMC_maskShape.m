function MASK = EMC_maskShape(SHAPE, SIZE, RADIUS, METHOD, OPTION)
%
% MASK = EMC_maskShape(SHAPE, SIZE, RADIUS, METHOD, OPTION)
% Create a real space mask of a given shape.
%
% Input:
%   SHAPE (str):                    'sphere', 'cylinder' or 'rectangle'.
%                                   NOTE: If 'rectangle' or 3d 'cylinder' (z only), the function
%                                         has to operate in pixel space. As such, in these cases,
%                                         the RADIUS and eventual shifts must only be integers,
%                                         and origin=0 is not allowed.
%
%   SIZE (int vector):              Size (in pixel) of the mask to compute; [x, y, z] or [x, y].
%                                   NOTE: [1, N] or [N, 1] is not allowed.
%
%   RADIUS (int vector):            Radius of the shape to compute; [x, y, z] or [x, y].
%                                   Should correspond to SIZE and be greater than 1.
%
%   METHOD (str):                   Device to use; 'gpu' or 'cpu'.
%
%   OPTION (cell | struct):         Optional parameters.
%                                   If cell: {field, value; ...}, note the ';' between parameters.
%                                   NOTE: Can be empty.
%                                   NOTE: Unknown fields will raise an error.
%
%     -> 'shift' (vector):          [x, y, z] or [x, y] shifts to apply in pixel.
%                                   Should correspond to SIZE.
%                                   NOTE: NaNs or Inf are not accepted.
%                                   NOTE: positive shifts translate the shape to the right.
%                                         negative shifts translate the shape to the left.
%                                   default = no shifts
%
%     -> 'origin' (int):            Origin convention - Center of the shape.
%                                   0, 1 or 2; see EMC_coordVectors for more details.
%                                   NOTE: origin=-1 is not allowed as this function is only for real
%                                         space masks.
%                                   default = 1
%
%     -> 'kernel' (bool|float|vector):
%                                   (gaussian) Kernel to convolve to the computed boolean mask.
%                                   If bool: apply or not the default kernel.
%                                   If float: the size of the kernel, in percentage (between 0 and 1) of
%                                             pixels, relative to the smallest axis. The minimum size of
%                                             the kernel is set to 9 pixels.
%                                   If vector: row vector (1xn, n is odd) used as separable kernel to
%                                              convolve with the boolean mask.
%                                   NOTE: For the default kernels, the variance of the gaussian is
%                                         automatically adjusted to keep the gaussian at ~0.5 at
%                                         the middle of the kernel.
%                                   NOTE: The kernel is imposing a minimum SIZE, which is equal to
%                                         ceil(kernelSize/2)*2+1. With the default kernels, the
%                                         minimum SIZE is 11.
%                                   default = 0.04 (4% of min(SIZE))
%
%     -> 'sym' (int):               Restrict the mask to the first asymmetric unit.
%                                   Should correspond to the central symmetry (positive int).
%                                   If 1, return the full mask.
%                                   default = 1
%
%     -> 'precision' (str):         Precision of the output MASK; 'single' or 'double'.
%                                   default = 'single'
%
% Output:
%   MASK (num array):               2d/3d mask with desired shape.
%
% Note:
%   - If a kernel is applied, the function forces the edges of the mask to go to zeros.
%     In that case, the pixels at the edges of the mask are set to 0.
%
% Example:
%   [MASK] = EMC_maskShape('cylinder', [128,128], [30,30], 'gpu', {})
%   [MASK] = EMC_maskShape('cylinder', [128,128,128], [30,30,30], 'gpu', {'shift',[0,10]})
%
% Other EMC-files required:
%   EMC_setPrecision, EMC_setMethod, EMC_is3d, EMC_getOption, EMC_gaussianKernel,
%   EMC_limits, EMC_resize, EMC_coordVectors, EMC_coordGrids, EMC_convn
%

% Created:  18Jan2020, R2019a
% Version:  v.1.0.0 switch to separable kernel (TF, 31Jan2020).
%           v.1.1.0 unittest (TF, 2Feb2020).
%           v.1.1.1 even kernels are now accepted. If you use your own vector kernel,
%                   it will be casted to the METHOD and precision before convolution (4Feb2020).
%           v.1.1.2 SIZE, RADIUS and shifts must be row vectors (columns are illegal now) (4Feb2020).
%

%% 
[SIZE, OPTION, flg] = checkIN(SIZE, RADIUS, METHOD, OPTION);

cutoffLow = 0.001;  % everything below this value is set to 0.

if flg.kernel
    if ~flg.ownKernel
        % Compute the size of the kernel in pixel.
        kernelSize = round((min(SIZE) * OPTION.kernel - 1) / 2) * 2 + 1;  % closest odd int
        if kernelSize < 9; kernelSize = 9; end  % at least 9 pixels

        % Keep the gaussian at ~0.5 at the middle of the roll off.
        middle = ceil(kernelSize / 2) / 2;
        sigma = sqrt(-1 * middle^2 / (2 * log(0.5)));
        OPTION.kernel = EMC_gaussianKernel([1,kernelSize], sigma, {'precision', OPTION.precision; ...
                                                                   'method', METHOD});
    else
        kernelSize = length(OPTION.kernel);
    end

    % To make the shape go to zeros at the edges of the mask, the function tape the mask
    % with zeros. This create a minimum size on the mask.
    if any(SIZE < ceil(kernelSize/2) * 2 + 1)
        error('EMC:kernel', 'with a kernel size of %d, the minimum SIZE is %d, got %s', ...
              kernelSize, ceil(kernelSize/2) * 2 + 1, mat2str(SIZE))
    end
else
    kernelSize = 0;
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
    
    % Cast from logical to float.
    MASK = EMC_setPrecision(MASK, OPTION.precision);
    
elseif strcmpi(SHAPE, 'cylinder')    
    % Adjust the radius for the same reason as explain with sphere/ellipsoids.
    RADIUS(1:2) = RADIUS(1:2) + ceil(kernelSize/2);

    % To compute a (3d) cylinder, this algorithm first compute a 2d sphere/ellipsoid and then
    % broadcast it along the Z (depth) axis.
    if flg.is3d
        % First compute the 2d ellipsoid with the x and y shifts.
        [vX, vY] = EMC_coordVectors(SIZE(1:2), METHOD, {'origin', OPTION.origin; ...
                                                        'precision', OPTION.precision;
                                                        'shift', OPTION.shift(1:2)});

        % Broadcaste in Z (the ellipsoid is invariant in Z <=> cylinder).
        RADIUS(3) = RADIUS(3) + floor(kernelSize/2);  % adjust Z for blurring.
        if flg.gpu
            vZ = zeros([1, 1, RADIUS(3) .*2 + 1], OPTION.precision, 'gpuArray');
        else
            vZ = zeros([1, 1, RADIUS(3) .*2 + 1], OPTION.precision);
        end
        MASK = (vX'./RADIUS(1)).^2 + (vY./RADIUS(2)).^2 + vZ <= 1;  % broadcast

        % Cast from logical to float.
        MASK = EMC_setPrecision(MASK, OPTION.precision);
    
        % Resize the cylinder to the desired SIZE in Z, taking into account the z shift.
        limits = EMC_limits(size(MASK), SIZE, {'origin', OPTION.origin; 'shift', [0,0,OPTION.shift(3)]});
        MASK = EMC_resize(MASK, limits, {'taper', false});

    else  % 2d Cylinder; this block is equivalent to SHAPE='sphere'.
        [vX, vY] = EMC_coordVectors(SIZE, METHOD, {'origin', OPTION.origin; ...
                                                   'precision', OPTION.precision; ...
                                                   'shift', OPTION.shift});
        MASK = (vX'./RADIUS(1)).^2 + (vY./RADIUS(2)).^2 <= 1;

        % Cast from logical to float.
        MASK = EMC_setPrecision(MASK, OPTION.precision);
    end

elseif strcmpi(SHAPE, 'rectangle')
    % Adjust the radius for the convolution.
    RADIUS = RADIUS + floor(kernelSize/2);

    if flg.gpu
        MASK = ones(RADIUS .*2 + 1, OPTION.precision, 'gpuArray');
    else
        MASK = ones(RADIUS .*2 + 1, OPTION.precision);
    end
    limits = EMC_limits(size(MASK), SIZE, {'origin', OPTION.origin; 'shift', OPTION.shift});
    MASK = EMC_resize(MASK, limits, {'taper', false});

else
    if ~(ischar(SHAPE) || isstring(SHAPE))
        error('EMC:SHAPE', "SHAPE should be string or char array, got %s", class(SHAPE))
    else
        error('EMC:SHAPE', "SHAPE should be 'sphere', 'cylinder' or 'rectangle', got %s", SHAPE)
    end
end

%% Restrict the mask to the first symmetry pair.
if flg.sym
    [~, angles, ~] = EMC_coordGrids('cylindrical', SIZE, METHOD, {'shift', OPTION.shift; ...
                                                                  'origin', OPTION.origin; ...
                                                                  'precision', OPTION.precision});
    sectorMax = 2*pi / OPTION.sym * 1.025;  % small overlap
    angles = (angles > (2*pi-sectorMax/2) | angles < sectorMax/2);
    MASK = MASK .* angles;
end

%% Finally, apply the blurring kernel to the binary mask.
% Use EMC_resize to tape the current 'binary' MASK edges with zeros. Convolved with the blurring
% kernel, it creates a nice roll off to zero at the edges of the MASK.
if flg.kernel
    % Force taper at the edges.
    MASK = EMC_resize(MASK, nan, {'force_taper', true; 'taper', zeros(1, ceil(kernelSize/2))});
    MASK = EMC_convn(MASK, OPTION.kernel); % this is the most expensive part
    MASK = MASK ./ max(MASK, [], 'all');  % rounding errors; max=1
    MASK(MASK <= cutoffLow) = 0;  % predictable end of the pass
end

end  % EMC_maskShape


function [SIZE, OPTION, flg] = checkIN(SIZE, RADIUS, METHOD, OPTION)

[flg.is3d, SIZE, ndim] = EMC_is3d(SIZE);
if SIZE(1) == 1 || SIZE(2) == 1
    error('EMC:SIZE', 'SIZE should be the size of a 2d or 3d array, got size %s', mat2str(SIZE))
end

if ~isnumeric(RADIUS) || ~isvector(RADIUS) || ~all(RADIUS > 1) || any(rem(RADIUS,1)) || ...
   ~isequal(size(SIZE), size(RADIUS)) || any(isinf(RADIUS))
    error('EMC:RADIUS', 'RADIUS should be a vector of %d integers greater than 1', ndim)
end

if strcmpi(METHOD, 'gpu')
    flg.gpu = true;
elseif strcmpi(METHOD, 'cpu')
    flg.gpu = false;
else
    error('EMC:METHOD', "METHOD should be 'gpu' or 'cpu'")
end

OPTION = EMC_getOption(OPTION, {'shift', 'origin', 'kernel', 'sym', 'precision'}, false);

% shift
if isfield(OPTION, 'shift')
    if ~isnumeric(OPTION.shift) || ~isrow(OPTION.shift)
        error('EMC:shift', 'OPTION.shift should be a vector of float|int, got %s', class(OPTION.shift))
    elseif any(isnan(OPTION.shift)) || any(isinf(OPTION.shift))
        error('EMC:shift', 'OPTION.shift should not contain NaNs or Inf, got %s', mat2str(OPTION.shift, 2))
   	elseif numel(OPTION.shift) ~= ndim
        error('EMC:shift', 'For a %dd SIZE, OPTION.shift should be a vector of %d float|int, got %s', ...
              ndim, ndim, mat2str(OPTION.shift, 2))
    end
else
    OPTION.shift = zeros(1, ndim);  % default
end

% origin
if isfield(OPTION, 'origin')
    if ~isscalar(OPTION.origin) || ~isnumeric(OPTION.origin) || ...
       ~(OPTION.origin == 0 || OPTION.origin == 1 || OPTION.origin == 2)
        % EMC_resize (used with 'rectangle' and 3d 'cylinders') will raise an error if origin=0
        error('EMC:origin', "OPTION.origin should be 0, 1, or 2, got %.04f", OPTION.origin)
    end
else
    OPTION.origin = 1;  % default
end

% precision
if isfield(OPTION, 'precision')
    if ~(ischar(OPTION.precision) || isstring(OPTION.precision))
        error('EMC:precision', 'OPTION.precision should be a string or char, got %s', class(OPTION.precision))
    elseif ~(strcmpi(OPTION.precision, 'single') || strcmpi(OPTION.precision, 'double'))
        error('EMC:precision', "OPTION.precision should be 'single' or 'double', got %s", OPTION.precision)
    end
else
    OPTION.precision = 'single';  % default
end

% kernel
flg.ownKernel = false;
if isfield(OPTION, 'kernel')
    if islogical(OPTION.kernel) && isscalar(OPTION.kernel)  % bool
        if OPTION.kernel
            flg.kernel = true;
            OPTION.kernel = 0.04;  % default
        else
            flg.kernel = false;
        end
    elseif isscalar(OPTION.kernel) && isnumeric(OPTION.kernel)  % float-nan
        if OPTION.kernel < 0 || OPTION.kernel > 1
            error('EMC:kernel', 'OPTION.kernel should be between 0 and 1, got %.04f', OPTION.kernel)
        end
        if OPTION.kernel == 0
            flg.kernel = false;
        else
            flg.kernel = true;
        end
    elseif isrow(OPTION.kernel) && isnumeric(OPTION.kernel)  % row vector
        OPTION.kernel = EMC_setMethod(EMC_setPrecision(OPTION.kernel, OPTION.precision), METHOD);
        flg.kernel = true;
        flg.ownKernel = true;
    else
        error('EMC:kernel', ['OPTION.kernel should be a boolean, a positive float between 0 and 1,', ...
              'or a row numeric vector, got %s', class(OPTION.taper)])
    end
else
     OPTION.kernel = 0.04;  % default
     flg.kernel = true;
end

% sym
if isfield(OPTION, 'sym')
    if isnumeric(OPTION.sym) && isscalar(OPTION.sym) && OPTION.sym > 0 && ...
       ~isinf(OPTION.sym) && ~rem(OPTION.sym, 1)
        if OPTION.sym == 1
            flg.sym = false;
        else
            flg.sym = true;
        end
    else
        error('EMC:sym', 'OPTION.sym should be a positive integer')
    end
else
    flg.sym = false;  % default
end

end  % checkIN
