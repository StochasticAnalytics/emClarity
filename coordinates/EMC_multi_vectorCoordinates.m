function [vX, vY, vZ] = EMC_multi_vectorCoordinates(SIZE, METHOD, OPTION)
% [vX, vY, vZ] = EMC_multi_vectorCoordinates(SIZE, METHOD, OPTION)
%
% Compute gridVectors.
%
% SIZE (vector):                Size of the vectors to compute; [x, y, z] or [x, y]
%                              	Values correspond to a number of pixels, as such, it should be integers.
%
% METHOD (str):                 Device to use; 'gpu' or 'cpu'.
%
% OPTION (cell | struct):       Optional parameters.
%                               If cell: {field,value ; ...}, note the ';' between parameters.
%                               NOTE: Can be empty.
%                               NOTE: Unknown fields will raise an error.
%
% -> 'shift' (vector):          [x, y, z] or [x, y] translations to apply (should correspond to SIZE)
%                               default = no shifts
%
% -> 'origin' (int):            Origin convention
%                               -1: zero frequency first (fft output)
%                               0: real origin (if even nb of pixel, the center is in between 2 pixels)
%                               1: right origin (extra pixel to the left; ceil((N+1)/2))
%                               2: left origin (extra pixel to the right)
%                               emClarity uses the right origin convention (ORIGIN=1).
%                               default = 1
%
% -> 'normalize' (bool):        Normalize the vectors between -0.5 and 0.5.
%                               default = false
%
% -> 'isotrope' (bool):         Stretch the vector values to the smallest dimensions.
%                               default = false
%
% -> 'half' (bool):             Compute half of the X vectors (useful for rfft)
%                               default = false
%
% -> 'precision' (str):         Precision of the vectors; 'single' or 'double'.
%                               default = 'single'
%
%--------
% EXAMPLE: [x,y,z] = EMC_multi_vectorCoordinates([10,9,8], 'gpu', {});
% EXAMPLE: [x, y]  = EMC_multi_vectorCoordinates([64,64], 'cpu', {'origin',-1 ; 'normalize',true});

%% checkIN
[SIZE, flg3d, ndim] = EMC_is3d(SIZE);
validateattributes(SIZE, {'numeric'}, {'vector', 'numel', ndim, 'nonnegative', 'integer'}, '', 'SIZE');

if ~(strcmpi('gpu', METHOD) || strcmpi('cpu', METHOD))
    error("method should be 'cpu' or 'gpu', got %s", METHOD)
end

% Extract optional parameters
OPTION = EMC_extract_option(OPTION, ...
                            {'origin', 'shift', 'normalize', 'isotrope', 'half', 'precision'}, ...
                            false);

if isfield(OPTION, 'origin')
    if ~(OPTION.origin == -1 || OPTION.origin == 0 || OPTION.origin == 1 || OPTION.origin == 2)
        error("origin should be 0, 1, 2, or -1, got %d", OPTION.origin)
    end
else
    OPTION.origin = 1;  % default
end

if isfield(OPTION, 'shift')
    validateattributes(OPTION.shift, {'numeric'}, {'vector', 'numel', ndim, 'finite', 'nonnan'}, ...
                       '', 'shift');
else
    OPTION.shift = zeros(1, ndim);  % default
end

if isfield(OPTION, 'normalize')
    if ~islogical(OPTION.normalize)
        error('normalize should be a boolean, got %s', class(OPTION.normalize))
    end
else
    OPTION.normalize = false;  % default
end

if isfield(OPTION, 'half')
    if ~islogical(OPTION.half)
        error('half should be a boolean, got %s', class(OPTION.half))
    end
else
    OPTION.half = false;  % default
end

if isfield(OPTION, 'isotrope')
    if ~islogical(OPTION.isotrope)
        error('isotrope should be a boolean, got %s', class(OPTION.isotrope))
    end
else
    OPTION.isotrope = false;  % default
end

if isfield(OPTION, 'precision')
    if ~(strcmpi(OPTION.precision, 'single') || strcmpi(OPTION.precision, 'double'))
        error("precision should be 'single' or 'double', got %s", OPTION.precision)
    end
else
    OPTION.precision = 'single';  % default
end

%% Create vectors with defined origin and shifts.
% For efficiency, shifts are applied to the boundaries and not to the vectors directly.
% On the other hand, normalization is done on the vectors, as rounding errors on the
% boundaries might lead to significative errors on the vectors.

% By default, the vectors are set to compute the real origin (origin = 0).
% To adjust for origin = 1|2, compute an offset to add to the vectors limits.
limits = zeros(2, ndim, OPTION.precision);
if OPTION.origin > 0 || OPTION.origin == -1
    if OPTION.origin == 2
        direction = -1;
    else  % origin = -1 or 1
        direction = 1;
    end
    for dim = 1:ndim
        if ~mod(SIZE(dim), 2)  % even dimensions: shift half a pixel
            limits(1, dim) = -SIZE(dim)/2 + 0.5 - direction * 0.5;
            limits(2, dim) =  SIZE(dim)/2 - 0.5 - direction * 0.5;
        else  % odd dimensions
            limits(1, dim) = -SIZE(dim)/2 + 0.5;
            limits(2, dim) =  SIZE(dim)/2 - 0.5;
        end
    end
else  % OPTIONAL.origin == 0
    limits(1,:) = -SIZE./2 + 0.5;
    limits(2,:) =  SIZE./2 - 0.5;
end

% centered
if OPTION.origin >= 0
    if (OPTION.half)
        if any(OPTION.shift)
            error('shifts are not allowed with half = true, got %s', mat2str(OPTION.shift, 2))
        end
        if OPTION.origin == 0 && ~mod(SIZE(1), 2)  % real center and even pixels
            if strcmpi(OPTION.precision, 'single')
                vX = 0.5:single((SIZE(1)/2));
            else
                vX = 0.5:(SIZE(1)/2);
            end
        else
            if strcmpi(OPTION.precision, 'single')
                vX = 0:single(floor(SIZE(1)/2));
            else
                vX = 0:floor(SIZE(1)/2);
            end
        end
    else
        vX = (limits(1, 1) - OPTION.shift(1)):(limits(2, 1) - OPTION.shift(1));
    end
    vY = (limits(1, 2) - OPTION.shift(2)):(limits(2, 2) - OPTION.shift(2));
    if (flg3d)
        vZ = (limits(1, 3) - OPTION.shift(3)):(limits(2, 3) - OPTION.shift(3));
    else
        vZ = nan;
    end

% not centered
else
    if any(OPTION.shift)
        error('shifts are not allowed with origin = -1, got %s', mat2str(OPTION.shift, 2))
    end
    if (OPTION.half)
        if strcmpi(OPTION.precision, 'single')
            vX = 0:single(floor(SIZE(1)/2));
        else
            vX = 0:floor(SIZE(1)/2);
        end
    else
        vX = [0:limits(2,1), limits(1,1):-1];
    end
    vY = [0:limits(2,2), limits(1,2):-1];
    if (flg3d); vZ = [0:limits(2,3), limits(1,3):-1]; else; vZ = nan; end
end

if strcmpi(METHOD, 'gpu')
    vX = gpuArray(vX);
    vY = gpuArray(vY);
    if (flg3d); vZ = gpuArray(vZ); end
elseif ~strcmpi(METHOD, 'cpu')
    error("METHOD must be 'gpu' or 'cpu', got %s", METHOD);
end

if (OPTION.isotrope)
    radius = min(abs(limits));
    radius_min = min(radius);
    vX = vX .* (radius_min / radius(1));
    vY = vY .* (radius_min / radius(2));
    if (flg3d); vZ = vZ .* (radius_min / radius(3)); else; vZ = nan; end
    if (OPTION.normalize)
        size_min = min(SIZE);
        vX = vX ./ size_min;
        vY = vY ./ size_min;
        if (flg3d); vZ = vZ ./ size_min; end
    end
elseif (OPTION.normalize)
    vX = vX ./ SIZE(1);
    vY = vY ./ SIZE(2);
    if (flg3d); vZ = vZ ./ SIZE(3); end
end

end  % EMC_multi_vectorCoordinates
