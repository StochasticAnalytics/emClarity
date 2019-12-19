function [MASK] = EMC_shapeMask(SHAPE, SIZE, RADIUS, METHOD, OPTIONAL)
%
% EMC_maskMask(SHAPE, SIZE, RADIUS, METHOD, OPTIONAL)
% Create a mask of a given shape.
%
% SHAPE (str | numeric):            'sphere', 'cylinder' or 'rectangle'.
%
% SIZE (vector):                    [x, y] or [x, y, z] dimension of mask.
%                                   Sould be a 2d/3d row vector of integers.
%
% RADIUS (vector):                  [x, y] or [x, y, z] radius in each dimension.
%                                   Should correspond to SIZE.
%
% METHOD (str):                     Device to use; 'gpu' or 'cpu'.
%
% OPTIONAL (cell | struct):         Optional parameters.
%                                   If cell: {field,value ; ...}, note the ';' between parameters.
%                                   NOTE: Can be empty.
%                                   NOTE: Unknown fields will raise an error.
%
%	-> 'shift' (vector):            [x, y] or [x, y, z] translations to apply.
%                                   Should correspond to SIZE.
%                                   default = no shifts
%
%  	-> 'origin' (int):              Origin convention - Center of rotation.
%                                   -1, 0, 1 or 2; see EMC_multi_gridVectors for more details.
%                                   NOTE: origin=0 is not allowed with SHAPE='rectangle'.
%                                   default = 1
%
%   -> 'taper' (bool|float):        Strength of the taper to apply at the mask edge.
%                                   If bool: apply or not the default taper.
%                                   If float: strength of the gaussian taper. Smaller values
%                                   gives shaper egdes. 'taper'=0 is equivalent to 'taper'=false.
%                                   NOTE: the size of the taper is scaled to the size of the
%                                         smallest axis, with a minimum of 7 pixels.
%                                   NOTE: the sigma of the gaussian is, sigma=taper*min(SIZE)/min(RADIUS).
%                                         As such, the taper is scaled to the size of the image and
%                                         invariant relative to the radius.
%                                   default = 0.02
%
%  	-> 'asym_restrict' (bool):      Turn on the asymmetric restriction.
%                                   Currently only for SHAPE = 'cylinder'.
%
%---------
% RETURN:                           MASK: 2d/3d mask
%
%---------
% EXAMPLE:                          [MASK] = EMC_shapeMask('cylinder', [128,128], [30,30], 'cpu', {})
%                                   [MASK] = EMC_shapeMask('cylinder', [128,128], [30,30], 'cpu', ...
%                                                          {'shift',[0,10]})
%

%% 
[SIZE, RADIUS, OPTIONAL, flg] = checkIN(SIZE, RADIUS, METHOD, OPTIONAL);
[vX, vY, vZ] = EMC_multi_vectorCoordinates(SIZE, METHOD, {'origin', OPTIONAL.origin; ...
                                                          'shift', OPTIONAL.shift});

cutoffLow = 0.025;
sigma = OPTIONAL.taper * min(SIZE) / min(RADIUS);

if strcmpi(SHAPE, 'sphere')
    if (flg.is3d)
        MASK = sqrt((vX'./RADIUS(1)).^2 + (vY./RADIUS(2)).^2 + reshape(vZ./RADIUS(3),1,1,[]).^2);
    else
        MASK = sqrt((vX'./RADIUS(1)).^2 + (vY./RADIUS(2)).^2);
    end
    solid = (MASK <= 1);
    if (flg.taper)
        MASK = exp(-1 .* ((MASK-1).^2) ./ (2*sigma^2));
        MASK(solid) = 1;
        MASK(MASK <= cutoffLow) = 0;
    else
        MASK(solid) = 1;
        MASK(~solid) = 0;
    end

elseif strcmpi(SHAPE, 'cylinder')
    if (flg.is3d)
        error('not supported yet.')
        vZ = reshape(abs(vZ), 1, 1, []);
        if (flg.gpu)
            gZ = zeros(SIZE, 'single', 'gpuArray') + abs(vZ);
            MASK = ones(SIZE, 'single', 'gpuArray');
        else
            gZ = zeros(SIZE, 'single') + abs(vZ);
            MASK = ones(SIZE, 'single');
        end
      	MASK = MASK .* sqrt((vX' ./ RADIUS(1)).^2 + (vY ./ RADIUS(2)).^2);
        solid = (MASK <= 1) & (gZ <= RADIUS(3));
    else
        MASK = sqrt((vX' ./ RADIUS(1)).^2 + (vY./RADIUS(2)).^2);
        solid = (MASK <= 1) ;
    end
    if (flg.taper)
        MASK = exp(-1 .* ((MASK-1).^2) ./ (2*sigma^2));
        MASK(solid) = 1;
        MASK(MASK <= cutoffLow) = 0;
    else
        MASK(solid) = 1;
        MASK(~solid) = 0;
    end

elseif strcmpi(SHAPE, 'rectangle')
    if (flg.taper)
        % taper should be large enough to go to the cutoffLow.
        taper = (1:ceil(min(RADIUS) * sqrt(-2 * log(cutoffLow) * sigma.^2))) ./ min(RADIUS);
        taper = exp(-1 .* (taper.^2) ./ (2*sigma^2));
        size_rectangle = (RADIUS + length(taper)) .* 2 + 1;
    else
        taper = false;
        size_rectangle = RADIUS .* 2 + 1;
    end
    
    if (flg.gpu)
        MASK = ones(size_rectangle, 'single', 'gpuArray');
    else
        MASK = ones(size_rectangle, 'single');
    end
    limits = EMC_multi_limits(size_rectangle, SIZE, {'origin', OPTIONAL.origin ; 'shift', OPTIONAL.shift});
    MASK = EMC_resize(MASK, limits, {'origin', OPTIONAL.origin; 'taper', taper});
    MASK(MASK <= cutoffLow) = 0;
else
    error("SHAPE should be 'sphere', 'cylinder' or 'rectangle', got %s", SHAPE)
end

end


function [SIZE, RADIUS, OPTIONAL, flg] = checkIN(SIZE, RADIUS, METHOD, OPTIONAL)

[flg.is3d, ndim] = EMC_is3d(SIZE);

validateattributes(SIZE, {'numeric'}, {'vector', 'nonnegative', 'integer'}, 'checkIN', 'SIZE');
validateattributes(RADIUS, {'numeric'}, {'vector', 'nonnegative', 'numel', ndim}, 'checkIN', 'RADIUS');

if strcmpi(METHOD, 'gpu')
    flg.gpu = true;
elseif strcmpi(METHOD, 'cpu')
    flg.gpu = false;
else
    error("METHOD should be 'gpu' or 'cpu', got %s", METHOD)
end

OPTIONAL = EMC_extract_optional(OPTIONAL, {'shift', 'origin', 'taper', 'asym_restrict'});

% shift
if isfield(OPTIONAL, 'shift')
    validateattributes(OPTIONAL.shift, {'numeric'}, ...
                       {'row', 'numel', ndim, 'finite', 'nonnan'}, 'checkIN', 'shift')
else
    OPTIONAL.shift = zeros(1, ndim);  % default
end
    
% origin
if isfield(OPTIONAL, 'origin')
    if ~(OPTIONAL.origin == -1 || OPTIONAL.origin == 0 || OPTIONAL.origin == 1 || OPTIONAL.origin == 2)
        error("origin should be 0, 1, 2, or -1, got %s", num2str(ORIGIN))
    end
else
    OPTIONAL.origin = 1;  % default
end

% taper
if isfield(OPTIONAL, 'taper')
    % bool
    if islogical(OPTIONAL.taper)
        if OPTIONAL.taper
            flg.taper = true;
            OPTIONAL.taper = 0.02;  % default
        else
            flg.taper = false;
            OPTIONAL.taper = 0;
        end
    elseif isfloat(OPTIONAL.taper) || isint(OPTIONAL.taper)
        if OPTIONAL.taper < 0
            error('taper should be positive, got %f', OPTIONAL.taper)
        end
        if OPTIONAL.taper == 0 || isnan(OPTIONAL.taper)
            flg.taper = false;
        else
            flg.taper = true;
        end
    else
        error('taper should be a boolean or a positive float, got %s', class(OPTIONAL.taper))
    end
else
     OPTIONAL.taper = 0.02;  % default
     flg.taper = true;
end

% asym_restrict
if isfield(OPTIONAL, 'asym_restrict')
    if ~islogical(OPTIONAL.asym_restrict)
        error('asym_restrict should be a boolean, got %s', class(OPTIONAL.asym_restrict))
    elseif OPTIONAL.asym_restrict && ~strcmpi(SHAPE, 'cylinder')
        error("asym_restrict only available for SHAPE = 'cylinder'")
    elseif OPTIONAL.asym_restrict
        error('asym_restrict not supported yet.')
    end
else
    OPTIONAL.asym_restrict = false;  % default
end

end  % checkIN