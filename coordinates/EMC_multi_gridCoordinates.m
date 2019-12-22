function [gX, gY, gZ, vX, vY, vZ] = EMC_multi_gridCoordinates(SIZE, METHOD, TRANS, varargin)
% [gX, gY, gZ, vX, vY, vZ] = EMC_multi_gridCoordinates(SIZE, METHOD, TRANS, varargin)
%
% Compute flow-field grids (gridCoordinates), for interpolations.
% WARNING: To compute grid masks, use EMC_multi_gridMasks.
% WARNING: To compute gridVectors, use EMC_multi_gridVectors.
%
% SIZE (vector):                Size of the grid to compute (x, y, z) or (x, y).
%                               Sould be a 2d/3d row vector of integers.
%
% METHOD (str):                 Device to use; 'gpu' or 'cpu'.
%
% TRANS (cell | str):           Transformation directives.
%   Syntax:
%       -> {}:                  no rotation, no shifts, no scaling
%       -> {field, value; ...}: Any optional fields. Fields that are not specified are set
%                               to their default value. Note the ';' between parameters.
%                               NOTE: unknown fields will raise an error.
%
%   Optional fields:
%       -> 'direction' (str):   direction convention (see BH_defineMatrix.m for more details),
%                               'forward'|'fwd' or 'inverse'|'inv'. If the grids are the query
%                               of interpolation, the direction should be 'inverse' to produce
%                               CCW rotation on the final image.
%                               default = 'inverse'
%
%       -> 'rotm' (vector):     Rotation to apply; 3D:3x3 or 2D:2x2 rotation matrices.
%                               If 2D: 3x3 rotation matrices are also accepted, but the extra
%                               row/column will be ignored.
%                               default = no rotation
%
%       -> 'shift' (vector):    (x, y, z) or (x, y) translations to apply; should correspond to SIZE.
%                               Shifts are applied BEFORE rotation and scaling.
%                               NOTE: Shifts do NOT change the center of rotation.
%                               default = no shifts
%
%       -> 'mag' (vector|int):  (x, y, z) or (x, y) scaling factor to apply; should correspond to SIZE.
%                               If only one int|float: isotropic magnification applied.
%                               NOTE: scaling applied BEFORE rotation.
%                               default = 1 (no scaling)
%
%       -> 'sym' (int):         Central symmetry. For each symmetry unit, a set of vectors
%                               and grids will be computed.
%                               default = 1 (no symmetry)
%
%       -> 'origin'(int):     	Origin convention - Center of rotation.
%                               -1: zero frequency first (fft output)
%                               0: real origin (if even nb of pixels, the center is in between 2 pixels)
%                               1: right origin (extra pixel to the left; ceil((N+1)/2))
%                               2: left origin (extra pixel to the right)
%                               emClarity uses the 'right' origin convention ('origin', 1).
%                               default = 1
%
%       -> 'offset' (vector):   Offset to apply; should correspond to SIZE. Offsets are used to adjust
%                               the center of rotation defined by 'origin'. 
%                               NOTE: this effectively apply a shift on both the vectors and the grids.
%                               NOTE: if there is no rotation or scaling to apply, this has no effect
%                                     on the final interpolated image.
%                               default = no offset
%
%       -> 'binary' (array):    Binary mask indicating the gridCoordinate pixels to ignore during
%                               transformation.
%
%       -> 'normalize' (bool):  Normalize the vectors and the grids between 0 and 1.
%                               default = false
%
% varagin (1x3 | 1x2 cell):     Use 3 pre-created gridVectors corresponding to vX, vY, vZ.
%                               NOTE: Half vectors are not accepted.
%
%---------
% RETURN:                       gX, gY, gZ are the gridCoordinates, vx, vY, vZ are the gridVectors.
%                               If TRANS.sym > 1, the grids are cells.
%
%---------
% EXAMPLE: [gX, gY, gZ, vX, vY, vZ] = EMC_multi_gridCoordinates([10,9], 'cpu', ...
%                                       {'rotm', [0,-1; 1,0]; 'shift', [0, 2]});
%
% See also EMC_multi_vectorCoordinates, EMC_multi_gridMasks

%% MAIN
% 1) To change the scale, the rotation matrix is scaled, nothing else.
% 2) To apply a shift, the vectors are shifted BUT the grids are NOT. As such, add shifts to the vectors
%    only once the grids are created.
% 3) To apply an offset (change the center of rotation), both the vectors and the grids are shifted.

if strcmpi(METHOD,'GPU')
  SIZE = gpuArray(single(SIZE));
else
  SIZE = single(SIZE);
end

[SIZE, METHOD, TRANS, flg, ndim] = checkIN(SIZE, METHOD, TRANS);

% Scale the rotation matrix.
% By default, the rotation is CCW, so the rotm must be transposed (it is not, but behaves like it).
% To effectively make the scaling BEFORE the rotation, flip the mag to [z, y, x] to match the transposed
% rotm.
if strcmpi(TRANS.direction, 'inverse')
    if any(TRANS.mag)
        TRANS.rotm = (eye(ndim) ./ flip(TRANS.mag)) * TRANS.rotm;
    end
else
    if any(TRANS.mag - 1)
        TRANS.rotm = (eye(ndim) .* TRANS.mag) * TRANS.rotm;
    end
    % Negate the offsets/shifts if 'forward'.
    TRANS.shift = TRANS.shift .* -1;
    TRANS.offset = TRANS.offset .* -1;
end

% Option to use pre-created vectors which is surprisingly expensive to create, otherwise make them.
if ~isempty(varargin)
    if isnumeric(varargin{1})
        vX = varargin{1};
        vY = varargin{2};
        if (flg.is3d)
            vZ = varargin{3};
        else
            vZ = nan;
        end
        
        % offsets
        if (flg.offset)
            vX = vX - TRANS.offset(1);
            vY = vY - TRANS.offset(2);
            if (flg.is3d); vZ = vZ - TRANS.offset(3); end
        end
        
        % normalize
        if (TRANS.normalize)
            vX = vX ./ SIZE(1);
            vY = vY ./ SIZE(2);
            if (flg.is3d); vZ = vZ ./ SIZE(3); end
        end
    else
        error('1st varargin should be numeric (gridVectors), got %s' + class(varargin{1}));
    end
else
    % Shift the vectors with the given offsets.
    [vX, vY, vZ] = EMC_multi_vectorCoordinates(...
        SIZE, METHOD, TRANS.offset, TRANS.origin, TRANS.normalize, false);
end

% The vectors are now generated, centered, and normalized if wished.

% gridCoordinates; Optionally evaluate only a smaller masked region.
if (flg.is3d)
    [X, Y, Z] = ndgrid(vX, vY, vZ);
    if (flg.binary)  
        X = X(binaryVol);
        Y = Y(binaryVol);
        Z = Z(binaryVol);
    end
else
    [X, Y] = ndgrid(vX, vY);
    if (flg.binary) 
        X = X(binaryVol);
        Y = Y(binaryVol);
    end
end

% Apply the shifts on the vectors only once the gridCoordinates are created.
if (flg.shift)
    vX = vX - TRANS.shift(1);
    vY = vY - TRANS.shift(2);
    if (flg.is3d); vZ = vZ - TRANS.shift(3); end
end

% Compute a grid for each symmetry unit.
if (flg.sym)
  	gX = {TRANS.sym, 1};
    gY = {TRANS.sym, 1};
    gZ = {TRANS.sym, 1};
    symInc = 360 ./ TRANS.sym;
end

% Rotate and scale the coordinates
for iSym = 1:TRANS.sym
   	if (flg.transform)
        % Only in plane symmetries considered anywhere
        % so inverse|forward shouldn't matter.
        if iSym > 1
            R = TRANS.rotm * BH_defineMatrix([iSym.*symInc,0,0],'Bah','inverse');
        else
            R = TRANS.rotm;
        end

        if (flg.is3d)
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
        YTrans = Y ;
        if (flg.is3d); ZTrans = Z; else; ZTrans = nan; end
    end

    % Only use as cell if symmetry is requested
    if (flg.sym)
        gX{iSym+1} = XTrans;
        gY{iSym+1} = YTrans;
        gZ{iSym+1} = ZTrans;
    else
        gX = XTrans;
        gY = YTrans;
        gZ = ZTrans;
    end
end % loop over symmetry units

clear X Y Z XTrans YTrans ZTrans
end  % EMC_multi_gridCoordinates


function [SIZE, METHOD, TRANS, flg, ndim] = checkIN(SIZE, METHOD, TRANS)
%% checkIN
% Sanity checks of BH_multi_gridCoordinates inputs.
[SIZE, flg.is3d, ndim] = EMC_is3d(SIZE);
validateattributes(SIZE, {'numeric'}, {'row', 'nonnegative', 'integer'}, 'checkIN', 'SIZE');

if ~(strcmpi(METHOD, 'gpu') || strcmpi(METHOD, 'cpu'))
     error("SYSTEM should be 'gpu' or 'cpu', got %s", METHOD)
end

% flags
flg.shift = false;
flg.transform = false;
flg.sym = false;
flg.binary = false;

% Extract optional parameters
TRANS = EMC_extract_option(TRANS, {'rotm', 'shift', 'mag', 'sym', 'direction', ...
                                   'origin', 'offset', 'binary', 'normalize'}, false);

if isfield(TRANS, 'rotm')
    validateattributes(TRANS.rotm, {'numeric'}, {'numel', ndim.^2, 'square'}, 'checkIN', 'rotm')
    % Most of the time, it will not be an identity matrix, so don't check and do transformation anyway.
    flg.transform = true;
else
    TRANS.rotm = eye(ndim);  % default
end

if isfield(TRANS, 'shift')
    validateattributes(TRANS.shift, {'numeric'}, ...
                       {'row', 'numel', ndim, 'finite', 'nonnan'}, 'checkIN', 'shift')
    if any(TRANS.shift)
        flg.shift = true;
    end
end

if isfield(TRANS, 'mag')
    if length(TRANS.mag) == 1
        TRANS.mag = zeros(1, ndim) + TRANS.mag;  % isotropic scaling
    else
        validateattributes(TRANS.mag, {'numeric'}, ...
                           {'vector', 'numel', ndim, 'finite', 'nonnan'}, 'checkIN', 'mag')
    end
    flg.transform = true;
else
    TRANS.mag = ones(1, ndim);  % default
end

if isfield(TRANS, 'sym')
    validateattributes(TRANS.sym, {'numeric'}, {'numel', ndim, 'integer', 'positive'}, 'checkIN', 'sym')
    if TRANS.sym ~= 1
        flg.sym = true;
    end
else
    TRANS.sym = 1;  % default
end

if isfield(TRANS, 'direction')
    if ~contains(['inverse', 'inv'], TRANS.direction)
        TRANS.direction = 'inverse';
    elseif ~contains([ 'forward', 'fwd'], TRANS.direction)
        TRANS.direction = 'forward';
    else
        error("direction should be 'forward' or 'inverse', got %s", TRANS.direction)
    end
else
    TRANS.direction = 'inverse';  % default
end

if isfield(TRANS, 'origin')
    if ~(TRANS.origin == -1 || TRANS.origin == 0 || TRANS.origin == 1 || TRANS.origin == 2)
        error("center should be 0, 1, 2, or -1, got %d", TRANS.origin)
    end
else
    TRANS.origin = 1;  % default
end

if isfield(TRANS, 'offset')
    validateattributes(TRANS.offset, {'numeric'}, ...
                       {'row', 'numel', ndim, 'finite', 'nonnan'}, 'checkIN', 'offset')
else
    TRANS.offset = zeros(1, ndim);  % default
end

if isfield(TRANS, 'binary')
    validateattributes(TRANS.offset, {'numeric'}, {'ndims', ndim, 'binary'}, 'checkIN', 'binary')
    flg.binary = true;
end

if isfield(TRANS, 'normalize')
    if ~islogical(TRANS.normalize)
        error('normalize should be a boolean, got %s', class(TRANS.normalize))
    end
else
    TRANS.normalize = false;  % default
end

end  % checkIN
