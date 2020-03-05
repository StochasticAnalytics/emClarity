function [gX, gY, gZ, vX, vY, vZ] = EMC_coordTransform(SIZE, METHOD, OPTION, varargin)
%
% [gX, gY, gZ, vX, vY, vZ] = EMC_coordTransform(SIZE, METHOD, OPTION, varargin)
% Compute flow-field grids (grid coordinates) and coordinate vectors used for interpolations.
%
% WARNING: The outputs are meant to be used for interpolation. To compute grid coordinates,
%          use EMC_coordGrids. To compute grid vectors, use EMC_coordVectors.
%
% Input:
%   SIZE (vector):           	Size (in pixel) of the grids to compute; [x, y, z] or [x, y].
%                               NOTE: [1, N] or [N, 1] is not allowed.
%
%   METHOD (str):            	Device to use; 'gpu' or 'cpu'.
%
%   OPTION (cell | struct):   	Optional parameters.
%                               If cell: {field,value ; ...}, note the ';' between parameters.
%                               NOTE: Can be empty.
%                               NOTE: Unknown fields will raise an error.
%
%     -> 'direction' (str):     Direction convention (see BH_defineMatrix.m for more details),
%                               'forward'|'fwd' or 'inverse'|'inv'. If the grids are the query
%                               of interpolation, the direction should be 'inverse' to produce
%                               CCW rotation on the final image.
%                               default = 'inverse'
%
%     -> 'rotm' (vector):       Rotation to apply; 3D:3x3 or 2D:2x2 rotation matrices.
%                               default = no rotation
%
%     -> 'shift' (vector):     	[x, y, z] or [x, y] translations to apply (should correspond to SIZE).
%                               Shifts are not allowed with half=true or origin=-1.
%                               NOTE: Shifts are applied BEFORE rotation and scaling.
%                               NOTE: Shifts do NOT change the center of rotation.
%                               NOTE: NaNs or Inf are not accepted.
%                               default = no shifts
%
%     -> 'mag' (vector|int):    [x, y, z] or [x, y] scaling to apply (should correspond to SIZE).
%                               If only one int|float: isotropic magnification applied.
%                               NOTE: scaling is applied BEFORE rotation.
%                               default = 1 (no scaling)
%
%     -> 'sym' (int):           Central symmetry. For each symmetry unit, a set of grids will be computed.
%                               default = 1 (no symmetry)
%
%     -> 'origin'(int):         Origin convention - Center of rotation.
%                               See EMC_coordVectors for more details.
%                               default = 1
%
%     -> 'offset' (vector):     [x, y, z] or [x, y] offset to apply (should correspond to SIZE).
%                               Offsets are used to adjust the center of rotation defined by 'origin'. 
%                               NOTE: this effectively apply a shift on both the vectors and the grids.
%                               NOTE: if there is no rotation or scaling to apply, this has no effect
%                                     on the final interpolated image.
%                               default = no offset
%
%     -> 'binary' (logical):    Binary mask (of size=SIZE) to apply on the grid coordinates before
%                               transformation.
%
%     -> 'normalize' (bool):    Normalize the vectors and the grids between -0.5 and 0.5.
%                               default = false
%
%     -> 'precision' (str):     Precision of the vectors and grids; 'single' or 'double'.
%                               default = 'single'
%
%   (optional)
%   varagin (3 row vectors):    Use 3 pre-created coordinate vectors corresponding to vX, vY and vZ.
%                               If SIZE correspond to a 3d array, vZ should be NaN.
%                               NOTE: Half vectors are not accepted.
%                               NOTE: For efficiency, this function do NOT check for NaNs or Infs values.
%
% Output:
%   gX, gY, gZ (num arrays):    Grid coordinates, either 2d or 3d.
%                               If 'sym' > 1, the grids are 1xn cells (n=sym) containing the grid
%                               coordinates for each symmetric unit.
%
%   vx, vY, vZ (num vectors):   Vector coordinates.
%
% Notes:
%   - To change the scale, the rotation matrix is scaled, nothing else.
%   - To apply a shift, the vectors are shifted BUT the grids are NOT. Therefore, shifts are applied
%     to the vectors only once the grids are created.
%   - To apply an offset (change the center of rotation), both the vectors and the grids are shifted.
%
% Example:
%   [gX, gY, ~, vX, vY, ~] = EMC_coordTransform([10,9], 'cpu', {'rotm', [0,-1; 1,0]; 'shift', [0, 2]});
%
% Other EMC-files required:
%   EMC_is3d, EMC_getOption, EMC_coordVectors
%
% See also EMC_coordVectors, EMC_coordGrids
%

% Created:  18Jan2020, R2019a
% Version:  v.1.0.  Rename (EMC_multi_gridCoordinates to EMC_coordTransform) and switch to new
%                   error identifier convention (TF, 30Jan2020).
%

%% MAIN
[SIZE, OPTION, flg, ndim] = checkIN(SIZE, METHOD, OPTION);

% Scale the rotation matrix.
% By default, the rotation is CCW, so the rotm must be transposed (it is done implicitely later).
% To effectively make the scaling BEFORE the rotation, flip the mag to [z, y, x] to match the transposed
% rotm.
if strcmpi(OPTION.direction, 'inverse')
    if any(OPTION.mag - 1)
        OPTION.rotm = (eye(ndim) ./ flip(OPTION.mag)) * OPTION.rotm;
    end
else  % forward
    if any(OPTION.mag - 1)
        OPTION.rotm = (eye(ndim) .* OPTION.mag) * OPTION.rotm;
    end
    % Negate the offsets/shifts if 'forward'.
    OPTION.shift = OPTION.shift .* -1;
    OPTION.offset = OPTION.offset .* -1;
end

% Option to use pre-created vectors, otherwise make them.
if ~isempty(varargin)
    if length(varargin) ~= 3
        error('EMC:varargin', ...
              'varargin should contain 3 row vectors, got %s elements', length(varargin))
    elseif ~flg.is3d
        if ~isscalar(varargin{3}) && ~isnan(varargin{3})
            error('EMC:varargin', 'For a 2d case, vZ should be NaN')
        else
            vX = EMC_setMethod(EMC_setPrecision(varargin{1}, OPTION.precision), METHOD);
            vY = EMC_setMethod(EMC_setPrecision(varargin{2}, OPTION.precision), METHOD);
        end
    else
        vX = EMC_setMethod(EMC_setPrecision(varargin{1}, OPTION.precision), METHOD);
       	vY = EMC_setMethod(EMC_setPrecision(varargin{2}, OPTION.precision), METHOD);
        vZ = EMC_setMethod(EMC_setPrecision(varargin{3}, OPTION.precision), METHOD);
    end

    if ~isnumeric(vX) || ~isrow(vX) || SIZE(1) ~= length(vX)
        error('EMC:varargin', 'varargin{1} (vX) should be a numeric row vector of %d elements', SIZE(1))
    elseif ~isnumeric(vY) || ~isrow(vY) || SIZE(2) ~= length(vY)
        error('EMC:varargin', 'varargin{2} (vY) should be a numeric row vector of %d elements', SIZE(2))
    elseif flg.is3d && ~isnumeric(vZ) || ~isrow(vZ) || SIZE(3) ~= length(vZ)
        error('EMC:varargin', 'varargin{3} (vZ) should be a numeric row vector of %d elements', SIZE(3))
    end

    % Apply offsets and|or normalize if whished. Note: shifts are not applied to vectors.
    if (flg.offset)
     	vX = vX - OPTION.offset(1);
      	vY = vY - OPTION.offset(2);
     	if flg.is3d; vZ = vZ - OPTION.offset(3); end
    end
    if (OPTION.normalize)
        vX = vX ./ SIZE(1);
      	vY = vY ./ SIZE(2);
        if flg.is3d; vZ = vZ ./ SIZE(3); end
    end

else  % varargin is empty
    OPTION = EMC_getOption(OPTION, {'offset', 'origin', 'normalize', 'precision'}, true);
    [vX, vY, vZ] = EMC_coordVectors(SIZE, METHOD, OPTION, false);
    % Note: shifts are not applied to vectors.
end

% The vectors are now ready.

% Compute the grids; Optionally evaluate only a smaller masked region.
if flg.is3d
    [X, Y, Z] = ndgrid(vX, vY, vZ);
    if flg.binary 
        X = X(binaryVol);
        Y = Y(binaryVol);
        Z = Z(binaryVol);
    end
else
    [X, Y] = ndgrid(vX, vY);
    if flg.binary
        X = X(binaryVol);
        Y = Y(binaryVol);
    end
end

% Apply the shifts on the vectors only once the grids are created.
if flg.shift
    vX = vX - OPTION.shift(1);
    vY = vY - OPTION.shift(2);
    if flg.is3d; vZ = vZ - OPTION.shift(3); end
end

% Compute a grid for each symmetry unit.
if flg.sym
    gX = {OPTION.sym, 1};
    gY = {OPTION.sym, 1};
    gZ = {OPTION.sym, 1};
    symInc = 360 ./ OPTION.sym;
end

% Rotate and scale the coordinates
for iSym = 1:OPTION.sym
    % Only in plane symmetries considered anywhere
    % so inverse|forward shouldn't matter.
    if iSym > 1
    	R = OPTION.rotm * BH_defineMatrix([iSym.*symInc,0,0],'Bah','inverse');
    else
      	R = OPTION.rotm;
    end
    
    if flg.transform || iSym > 1
        if flg.is3d
            % CCW rotation by default (rotm')
            XTrans = X.*R(1,1) + Y.*R(2,1) + Z.*R(3,1);
            YTrans = X.*R(1,2) + Y.*R(2,2) + Z.*R(3,2);
            ZTrans = X.*R(1,3) + Y.*R(2,3) + Z.*R(3,3);
        else
            % CCW rotation by default (rotm')
            XTrans = X.*R(1,1) + Y.*R(2,1);
            YTrans = X.*R(1,2) + Y.*R(2,2);
            ZTrans = nan;
        end
    else
        XTrans = X;
        YTrans = Y;
        if flg.is3d; ZTrans = Z; else; ZTrans = nan; end
    end

    % Only use as cell if symmetry is requested
    if flg.sym
        gX{iSym+1} = XTrans;  % I [TF] hope this doesn't generate a copy.
        gY{iSym+1} = YTrans;
        gZ{iSym+1} = ZTrans;
    else
        gX = XTrans;
        gY = YTrans;
        gZ = ZTrans;
    end
end % loop over symmetry units

clear X Y Z XTrans YTrans ZTrans
end  % EMC_coordTransform


function [SIZE, OPTION, flg, ndim] = checkIN(SIZE, METHOD, OPTION)

[flg.is3d, SIZE, ndim] = EMC_is3d(SIZE);

if ~(strcmpi(METHOD, 'gpu') || strcmpi(METHOD, 'cpu'))
    if isstring(METHOD) || ischar(METHOD) 
        error('EMC:METHOD', "SYSTEM should be 'gpu' or 'cpu', got %s", METHOD)
    else
        error('EMC:METHOD', "SYSTEM should be 'gpu' or 'cpu', got %s", clas(METHOD))
    end
end

% flags
flg.shift = false;
flg.transform = false;
flg.sym = false;
flg.binary = false;

% Extract optional parameters
OPTION = EMC_getOption(OPTION, {'rotm', 'shift', 'mag', 'sym', 'direction', ...
                                'origin', 'offset', 'binary', 'normalize', 'precision'}, false);

% rotm
if isfield(OPTION, 'rotm')
    if ~isnumeric(OPTION.rotm) || ~ismatrix(OPTION.rotm) 
        error('EMC:rotm', 'rotm should be a %dx%d numeric matrix, got %s', ...
              ndim, ndim, class(OPTION.rotm))
    elseif numel(OPTION.rotm) == ndim^2
        error('EMC:rotm', 'rotm should be a %dx%d numeric matrix, got size:%s', ...
              ndim, ndim, mat2str(size(OPTION.rotm)))
    end
    % Most of the time, it will not be an identity matrix, so don't check and do transformation anyway.
    flg.transform = true;
else
    OPTION.rotm = eye(ndim);  % default
end

% shift
if isfield(OPTION, 'shift')
    if ~isnumeric(OPTION.shift) || ~isvector(OPTION.shift)
        error('EMC:shift', ...
              'shift should be a vector of float|int, got %s', class(OPTION.shift))
    elseif any(isnan(OPTION.shift)) || any(isinf(OPTION.shift))
        error('EMC:shift', ...
              'shift should not contain NaNs or Inf, got %s', mat2str(OPTION.shift, 2))
    elseif numel(OPTION.shift) ~= ndim
        error('EMC:shift', ...
              'For a %dd SIZE, shift should be a vector of %d float|int, got %s', ...
              ndim, ndim, mat2str(OPTION.shift, 2))
    elseif any(OPTION.shift)
        flg.shift = true;
    end
else
    OPTION.shift = zeros(1, ndim);  % default
end

% mag
if isfield(OPTION, 'mag')
    if ~isnumeric(OPTION.mag)
        error('EMC:mag', ...
              'mag should be a numeric scalar or vector, got %s', class(OPTION.mag))
    elseif isvector(OPTION.mag)
        if length(OPTION.mag) ~= ndim
            error('EMC:mag', ...
                  'mag should be a vector of %d elements, got %d elements', ndim, length(OPTION.mag))
      	elseif any(isnan(OPTION.mag)) || any(isinf(OPTION.mag))
            error('EMC:mag', ...
                  'mag should not have any nan nor inf, got:%s', mat2str(OPTION.mag))
        end
        flg.transform = true;
    elseif isscalar(OPTION.mag) && ~any(isnan(OPTION.mag)) || ~any(isinf(OPTION.mag))
        OPTION.mag = zeros(1, ndim) + OPTION.mag;  % isotropic scaling
        flg.transform = true;
    else
       	error('EMC:mag', ...
              'mag should be a numeric scalar or a numeric vector of %d elements', ndim)
    end
else
    OPTION.mag = ones(1, ndim);  % default
end

% sym
if isfield(OPTION, 'sym')
    if ~isnumeric(OPTION.sym) || ~isscalar(OPTION.sym) || OPTION.sym < 1 || rem(OPTION.sym, 1)
        error('EMC:sym', ...
              'sym should be a positive integer')
    elseif OPTION.sym ~= 1
        flg.sym = true;
    end
else
    OPTION.sym = 1;  % default
end

% direction
if isfield(OPTION, 'direction')
    if any(strmcpi(['inverse', 'inv'], OPTION.direction))
        OPTION.direction = 'inverse';
    elseif any(strmcpi([ 'forward', 'fwd'], OPTION.direction))
        OPTION.direction = 'forward';
    else
        error('EMC:direction', "direction should be 'forward' or 'inverse'")
    end
else
    OPTION.direction = 'inverse';  % default
end

% origin
if isfield(OPTION, 'origin')
    if ~isnumeric(OPTION.origin) || ~isscalar(OPTION.origin) || ...
       ~(OPTION.origin == 1 || OPTION.origin == -1 || OPTION.origin == 0 || OPTION.origin == 2)
        error('EMC:origin', 'origin should be 0, 1, 2, or -1, got %d', OPTION.origin)
    end
else
    OPTION.origin = 1;  % default
end

% offset
if isfield(OPTION, 'offset')
    if ~isnumeric(OPTION.offset) || ~isvector(OPTION.offset)
        error('EMC:offset', ...
              'offset should be a vector of float|int, got %s', class(OPTION.offset))
    elseif any(isnan(OPTION.offset)) || any(isinf(OPTION.offset))
        error('EMC:offset', ...
              'offset should not contain NaNs or Inf, got %s', mat2str(OPTION.offset, 2))
    elseif numel(OPTION.offset) ~= ndim
        error('EMC:offset', ...
              'For a %dd SIZE, offset should be a vector of %d float|int, got %s', ...
              ndim, ndim, mat2str(OPTION.offset, 2))
    end
else
    OPTION.offset = zeros(1, ndim);  % default
end

% binary
if isfield(OPTION, 'binary')
    if ~islogical(OPTION.binary) || size(OPTION.binary) ~= SIZE
        error('EMC:binary', 'binary should be a logical of size:%s', mat2str(SIZE))
    else
        flg.binary = true;
    end
end

% normalize
if isfield(OPTION, 'normalize')
    if ~islogical(OPTION.normalize) || ~isscalar(OPTION.normalize)
        error('EMC:normalize', 'normalize should be a boolean,')
    end
else
    OPTION.normalize = false;  % default
end

% precision is checked by EMC_coordVectors or EMC_setPrecision
if ~isfield(OPTION, 'precision')
    OPTION.precision = 'single';  % default
end

end  % checkIN
