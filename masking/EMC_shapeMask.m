function [MASK] = EMC_shapeMask(SHAPE, SIZE, RADIUS, METHOD, OPTION)
%
% EMC_maskMask(SHAPE, SIZE, RADIUS, METHOD, OPTION)
% Create a real space mask of a given shape.
%
% SHAPE (str | numeric):            'sphere', 'cylinder' or 'rectangle'.
%
% SIZE (vector):                    [x, y] or [x, y, z] dimension of mask.
%                                   Sould be a 2d/3d row vector of integers.
%
% RADIUS (vector):                  [x, y] or [x, y, z] radius in each dimension.
%                                   Should correspond to SIZE.
%
% METHOD (str):                     Device of the output MASK; 'gpu' or 'cpu'.
%
% OPTION (cell | struct):           Optional parameters.
%                                   If cell: {field,value ; ...}, note the ';' between parameters.
%                                   NOTE: Can be empty.
%                                   NOTE: Unknown fields will raise an error.
%
%	-> 'shift' (vector):            [x, y] or [x, y, z] translations to apply.
%                                   Should correspond to SIZE.
%                                   default = no shifts
%
%  	-> 'origin' (int):              Origin convention - Center of rotation.
%                                   0, 1 or 2; see EMC_multi_gridVectors for more details.
%                                   NOTE: origin=-1 is not allowed as this function is only for real
%                                         space masks.
%                                   NOTE: origin=0 is currently not allowed for SHAPE = 'rectangle' and
%                                         'cylinder'.
%                                   default = 1
%
%   -> 'taper' (bool|float|matrix): Gaussian taper to apply to the boolean mask.
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
%  	-> 'sym' (int):                 Restrict the mask to the first asymmetric unit.
%                                   Should correspond to the central symmetry (positive int).
%                                   If 1, return the full mask.
%                                   default = 1
%
%   -> 'precision' (str):           Precision of the output MASK; 'single' or 'double'.
%                                   default = 'single'
%
%---------
% RETURN:                           MASK: 2d/3d mask
%
%---------
% NOTE:                             If taper, the function forces the edges of the mask to be blurred,
%                                   which is equivalent to EMC_resize:force_taper option but with the
%                                   desired roll off.
%
%---------
% EXAMPLE:                          [MASK] = EMC_shapeMask('cylinder', [128,128], [30,30], 'cpu', {})
%                                   [MASK] = EMC_shapeMask('cylinder', [128,128], [30,30], 'cpu', ...
%                                                          {'shift',[0,10]})
%

%% 
[SIZE, RADIUS, OPTION, flg, ndim] = checkIN(SIZE, RADIUS, METHOD, OPTION);

cutoffLow = 0.001;  % everything below this value is set to 0.

if (flg.taper)
    if ~(flg.ownTaper)
        % Compute the size of the kernel in pixel.
        smallestDimension = min(SIZE);
        kernelSize = round((smallestDimension * OPTION.taper - 1) / 2) * 2 + 1;  % closest odd int
        if kernelSize < 9; kernelSize = 9; end  % at least 8 pixels.

        % Keep the gaussian at ~0.5 at the middle of the taper.
        middle = ceil(kernelSize / 2) / 2;
        sigma = sqrt(-1 * middle^2 / (2 * log(0.5)));
        kernel = EMC_multi_gaussianKernel(zeros(1,ndim) + kernelSize, sigma, ...
                                          {'precision', OPTION.precision});
    end
end

if strcmpi(SHAPE, 'sphere')
    % Adjust the radius for the convolution.
    RADIUS = RADIUS + ceil(kernelSize/2);
    
    [vX, vY, vZ] = EMC_multi_vectorCoordinates(SIZE, METHOD, {'origin', OPTION.origin; ...
                                                              'shift', OPTION.shift; ...
                                                              'precision', OPTION.precision});
    if (flg.is3d)
        MASK = (vX'./RADIUS(1)).^2 + (vY./RADIUS(2)).^2 + reshape(vZ./RADIUS(3),1,1,[]).^2;
    else
        MASK = (vX'./RADIUS(1)).^2 + (vY./RADIUS(2)).^2;
    end

    % Binary mask
    if strcmpi(OPTION.precision, 'single')
        MASK = single((MASK <= 1));
    else
        MASK = double((MASK <= 1));
    end

elseif strcmpi(SHAPE, 'cylinder')
    % Adjust the radius for the convolution.
    RADIUS(1:2) = RADIUS(1:2) + ceil(kernelSize/2);

    if (flg.is3d)
        % First compute the circle (the shifts are dealt afterwards)
        [vX, vY, ~] = EMC_multi_vectorCoordinates(SIZE(1:2), METHOD, {'origin', OPTION.origin; ...
                                                                      'precision', OPTION.precision});
        % Invariant in Z
        RADIUS(3) = RADIUS(3) + floor(kernelSize/2);  % adjust Z
        if (flg.gpu)
            vZ = zeros([1, 1, RADIUS(3) .*2 + 1], OPTION.precision, 'gpuArray');
        else
            vZ = zeros([1, 1, RADIUS(3) .*2 + 1], OPTION.precision);
        end
      	MASK = (vX' ./ RADIUS(1)).^2 + (vY ./ RADIUS(2)).^2 + vZ;  % broadcast

        % Binary mask
        if strcmpi(OPTION.precision, 'single')
            MASK = single((MASK <= 1));
        else
            MASK = double((MASK <= 1));
        end

        % Resize the cylinder to the desired SIZE in Z, taking into account the shifs.
        limits = EMC_multi_limits(size(MASK), SIZE, {'origin', OPTION.origin; 'shift', OPTION.shift});
        MASK = EMC_resize(MASK, limits, {'taper', false});

    else  % cylinder 2d
        [vX, vY] = EMC_multi_vectorCoordinates(SIZE, METHOD, {'origin', OPTION.origin; ...
                                                              'precision', OPTION.precision; ...
                                                              'shift', OPTION.shift});
        MASK = (vX' ./ RADIUS(1)).^2 + (vY./RADIUS(2)).^2;

        % Binary mask
        if strcmpi(OPTION.precision, 'single')
            MASK = single((MASK <= 1));
        else
            MASK = double((MASK <= 1));
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
    limits = EMC_multi_limits(size(MASK), SIZE, {'origin', OPTION.origin; 'shift', OPTION.shift});
    MASK = EMC_resize(MASK, limits, {'taper', false});
else
    error("SHAPE should be 'sphere', 'cylinder' or 'rectangle', got %s", SHAPE)
end

% Restrict the mask to the first symmetry pair.
if (flg.sym)
    [~, angles, ~] = EMC_multi_gridMasks('cylindrical', SIZE, METHOD, {'shift', OPTION.shift; ...
                                                                       'origin', OPTION.origin; ...
                                                                       'precision', OPTION.precision});
    sectorMax = 2 * pi / OPTION.sym * 1.025;  % small overlap
    angles = (angles > (2*pi-sectorMax/2) | angles < sectorMax/2);
    MASK = MASK .* angles;
end

% Finally, apply the kernel to the binary mask.
% I've tried doing the convolution in Fourier space,
% but in most situations it is considerably slower...
if (flg.taper)
    % Force taper at the edges.
    MASK = EMC_resize(MASK, zeros(1,ndim*2), {'force_taper', true; 'taper', zeros(1, ceil(kernelSize/2))});
    MASK = convn(MASK, kernel, 'same');  % this is the most expensive part
    MASK = MASK ./ max(MASK(:));  % rounding errors; max=1
    MASK(MASK <= cutoffLow) = 0;  % predictable end of the pass
end

end  % end EMC_shapeMask


function [SIZE, RADIUS, OPTION, flg, ndim] = checkIN(SIZE, RADIUS, METHOD, OPTION)

[SIZE, flg.is3d, ndim] = EMC_is3d(SIZE);

validateattributes(SIZE, {'numeric'}, {'vector', 'nonnegative', 'integer'}, 'checkIN', 'SIZE');
validateattributes(RADIUS, {'numeric'}, {'vector', 'nonnegative', 'numel', ndim}, 'checkIN', 'RADIUS');

if strcmpi(METHOD, 'gpu')
    flg.gpu = true;
elseif strcmpi(METHOD, 'cpu')
    flg.gpu = false;
else
    error("METHOD should be 'gpu' or 'cpu', got %s", METHOD)
end

OPTION = EMC_extract_option(OPTION, {'shift', 'origin', 'taper', 'sym', 'precision'}, false);

% shift
if isfield(OPTION, 'shift')
    validateattributes(OPTION.shift, {'numeric'}, ...
                       {'row', 'numel', ndim, 'finite', 'nonnan'}, 'checkIN', 'shift')
else
    OPTION.shift = zeros(1, ndim);  % default
end
    
% origin
if isfield(OPTION, 'origin')
    if ~(OPTION.origin == 0 || OPTION.origin == 1 || OPTION.origin == 2)
        % EMC_resize (used with 'rectangle' and 3d 'cylinders') will raise an error if origin=0
        error("origin should be 0, 1, or 2, got %.04f", OPTION.origin)
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
    if (isfloat(OPTION.sym) || isinteger(OPTION.sym)) && OPTION.sym > 0
        if OPTION.sym == 1
            flg.sym = false;
        else
            flg.sym = true;
        end
    else
        error('sym should be a positive, non zero integer')
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
