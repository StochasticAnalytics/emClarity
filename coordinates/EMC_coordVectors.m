function [vX, vY, vZ] = EMC_coordVectors(SIZE, METHOD, OPTION)
%
% [vX, vY, vZ] = EMC_coordVectors(SIZE, METHOD, OPTION)
% Compute coordinate vectors.
%
% Input:
%   SIZE (int vector):          Size (in pixel) of the vectors to compute; [x, y, z] or [x, y].
%                               NOTE: [1, N] is not equivalent to [N, 1].
%                                     See 'Output' and 'Note' for more details.
%
%   METHOD (str):              	Device to use; 'gpu' or 'cpu'.
%
%   OPTION (cell | struct):   	Optional parameters.
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
%                               -1: zero frequency first (fft output)
%                               0: real origin (if even nb of pixel, the center is in between 2 pixels)
%                               1: right origin (extra pixel to the left; ceil((N+1)/2))
%                               2: left origin (extra pixel to the right)
%                               emClarity uses the 'right' origin convention (ORIGIN=1).
%                               default = 1
%
%     -> 'normalize' (bool):   	Normalize the vectors between -0.5 and 0.5.
%                               default = false
%
%     -> 'isotrope' (bool):    	Stretch the vector values to the smallest dimensions.
%                               default = false
%
%     -> 'half' (bool):         Compute 'half' of the vX vector. Originally made for rfft (the FT
%                               of a real function is Hermitian).
%                               default = false
%
%     -> 'precision' (str):   	Precision of the vectors; 'single' or 'double'.
%                               default = 'single'
%
% Output:
%   [vX, vY, vZ] (vectors):     Coordinate (row) vectors (x, y and z).
%                               The order of these vectors corresponds directly to the input SIZE.
%                               Dimensions equal to 1 will output NaN.
%
% Note:
%   - Any axis with one pixel length will have a corresponding vector equal to NaN.
%     As such, SIZE=[1,N] will output vX = vZ = NaN and length(vY) = N.
%     On the other hand, SIZE=[N,1], will output vY = vZ = NaN and length(vX) = N.
%
% Example:
%   [x,y,z] = EMC_coordVectors([10,9,8], 'gpu', {});
%   [x,y]   = EMC_coordVectors([64,64], 'cpu', {'origin', -1; 'normalize',true});
%   [x]     = EMC_coordVectors([64, 1], 'cpu', {'origin', -1; 'normalize',true});
%   [~,y]   = EMC_coordVectors([1, 64], 'cpu', {'origin', -1; 'normalize',true});
%
% Other EMC-files required:
%   EMC_is3d, EMC_getOption, EMC_setPrecision
%
% See also EMC_is3d
%

% Created:  18Jan2020, R2019a
% Version:  v.1.0   New SIZE convention (see EMC_is3d).
%                   Now follow MATLAB convention for scalars and vectors (TF, 23Jan2020)
%           v.1.1   Rename (EMC_coordVectors to EMC_coordVectors) and unittest (TF, 24Jan2020).
%

%% checkIN
[is3d, SIZE, ndim] = EMC_is3d(SIZE);

% Extract optional parameters
OPTION = EMC_getOption(OPTION, {'origin', 'shift', 'normalize', 'isotrope', 'half', 'precision'}, false);

if isfield(OPTION, 'origin')
    if ~isscalar(OPTION.origin) || ~isnumeric(OPTION.origin)
        error('EMC:origin', 'OPTION.origin should be an integer, got %s of size: %s', ...
              class(OPTION.origin), mat2str(size(OPTION.origin)))
    elseif OPTION.origin ~= 1 && OPTION.origin ~= -1 && OPTION.origin ~= 0 && OPTION.origin ~= 2
        error('EMC:origin', 'OPTION.origin should be 0, 1, 2, or -1, got %d', OPTION.origin)
    end
else
    OPTION.origin = 1;  % default
end

if isfield(OPTION, 'half')
    if ~islogical(OPTION.half) || ~isscalar(OPTION.half)
        error('EMC:half', 'OPTION.half should be a boolean, got %s', class(OPTION.half));
    end
else
    OPTION.half = false;  % default
end

if isfield(OPTION, 'shift')
    if ~isnumeric(OPTION.shift) || ~isrow(OPTION.shift)
        error('EMC:shift', ...
              'OPTION.shift should be a vector of float|int, got %s', class(OPTION.shift))
    elseif any(isnan(OPTION.shift)) || any(isinf(OPTION.shift))
        error('EMC:shift', ...
              'OPTION.shift should not contain NaNs or Inf, got %s', mat2str(OPTION.shift, 2))
    elseif numel(OPTION.shift) ~= ndim
        error('EMC:shift', ...
              'For a %dd SIZE, OPTION.shift should be a vector of %d float|int, got %s', ...
              ndim, ndim, mat2str(OPTION.shift, 2))
    elseif (OPTION.half || OPTION.origin == -1) && any(OPTION.shift)
        error('EMC:shift', ...
              'OPTION.shifts are not allowed with half=true or origin=-1 , got %s', mat2str(OPTION.shift, 2))
    end
else
    OPTION.shift = zeros(1, ndim);  % default
end

if isfield(OPTION, 'normalize')
    if ~islogical(OPTION.normalize) || ~isscalar(OPTION.normalize)
        error('EMC:normalize', 'OPTION.normalize should be a boolean, got %s', class(OPTION.normalize))
    end
else
    OPTION.normalize = false;  % default
end

if isfield(OPTION, 'isotrope')
    if ~islogical(OPTION.isotrope) || ~isscalar(OPTION.isotrope)
        error('EMC:isotrope', 'OPTION.isotrope should be a boolean, got %s', class(OPTION.isotrope))
    end
else
    OPTION.isotrope = false;  % default
end

if isfield(OPTION, 'precision')
    if ~(ischar(OPTION.precision) || isstring(OPTION.precision)) || ...
       ~strcmpi(OPTION.precision, 'single') && ~strcmpi(OPTION.precision, 'double')
        error('EMC:precision', "OPTION.precision should be 'single' or 'double'")
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
if OPTION.origin ~= 0
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

% Any dimension with a size of 1 is not computed and has its corresponding vector equals to NaN.
isDim = SIZE > 1;

% real space
if OPTION.origin >= 0
    if isDim(1)
        if OPTION.half
            if OPTION.origin == 0 && ~mod(SIZE(1), 2)  % real center and even pixels
                vX = 0.5:EMC_setPrecision((SIZE(1)/2), OPTION.precision);
            else
                vX = 0:EMC_setPrecision(floor(SIZE(1)/2), OPTION.precision);
            end
        else
            vX = (limits(1, 1) - OPTION.shift(1)):(limits(2, 1) - OPTION.shift(1));
        end
    else
        vX = nan;
    end
    
    if isDim(2); vY = (limits(1, 2) - OPTION.shift(2)):(limits(2, 2) - OPTION.shift(2)); else; vY = nan; end
    if is3d;     vZ = (limits(1, 3) - OPTION.shift(3)):(limits(2, 3) - OPTION.shift(3)); else; vZ = nan; end
    
% reciprocal space
else
    if isDim(1)
        if (OPTION.half)
            vX = 0:EMC_setPrecision(floor(SIZE(1)/2), OPTION.precision);
        else
            vX = [0:limits(2,1), limits(1,1):-1];
        end
    else
        vX = nan;
    end
    
    if isDim(2); vY = [0:limits(2,2), limits(1,2):-1]; else; vY = nan; end
    if is3d;     vZ = [0:limits(2,3), limits(1,3):-1]; else; vZ = nan; end
end

% In my case [TF], it is faster to create on host and then push to device for vectors < 8000 elements.
if strcmpi(METHOD, 'gpu')
    if isDim(1); vX = gpuArray(vX); end
    if isDim(2); vY = gpuArray(vY); end
    if is3d;     vZ = gpuArray(vZ); end
elseif ~(ischar(METHOD) || isstring(METHOD)) || ~strcmpi(METHOD, 'cpu')
    error('EMC:METHOD', "METHOD must be 'gpu' or 'cpu'");
end

if (OPTION.isotrope)
    radius = max(abs(limits));
    radius_min = min(radius);
    if OPTION.normalize; radius = radius .* min(SIZE); end
    if isDim(1); vX = vX .* (radius_min / radius(1)); end
    if isDim(2); vY = vY .* (radius_min / radius(2)); end
    if is3d;     vZ = vZ .* (radius_min / radius(3)); end
elseif (OPTION.normalize)
    if isDim(1); vX = vX ./ SIZE(1); end
    if isDim(2); vY = vY ./ SIZE(2); end
    if is3d;     vZ = vZ ./ SIZE(3); end
end

end  % EMC_coordVectors
