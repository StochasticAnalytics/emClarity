function [vX, vY, vZ] = EMC_multi_vectorCoordinates(SIZE, METHOD, SHIFT, ORIGIN, ...
                                                    flgNormalize, ...
                                                    flgHalf)
% [vX, vY, vZ] = EMC_multi_vectorCoordinates(SIZE, METHOD, SHIFT, ORIGIN, flgNormalize, flgHalf)
% Compute gridVectors.
%
% SIZE (int/float vector):      Size of the vectors to compute; [x, y, z] or [x, y]
%                              	Values correspond to a number of pixels, as such, the decimal
%                             	should be 0.
%
% METHOD (str):               	Device to use; 'gpu' or 'cpu'
%
% SHIFT (int/float vector):     [x, y, z] or [x, y] translations to apply (should correspond to size)
%
% ORIGIN (int):                	Origin convention
%                               -1: zero frequency first (fft output)
%                               0: real origin (if even nb of pixel, the center is in between 2 pixels)
%                               1: right origin (extra pixel to the left; ceil((N+1)/2))
%                               2: left origin (extra pixel to the right)
%                               emClarity uses the right origin convention (ORIGIN=1).
%
% flgNormalize (bool):        	Normalize the vectors between 0 and 1
% flgHalf (bool)              	Compute half of the vectors (rfft)
%
% EXAMPLE: [x,y,z] = EMC_multi_vectorCoordinates([10,9,8], 'cpu', [0,0,0], 1, true, false)

%% TODO
% 1) For Ben, to check:
%       Are you sure about the vector coordinates when flgOrigin == 0?
%       Here is what you do:
%       flgOrigin = 0: 6pixels(0, 1, 2, 3, -2, -1) but shouldn't it be 6pixels(0, 1, 2, -3, -2, -1)?
%       For now, I did as I think is correct. I'll wait for your answer.
%
% 2) For Ben, to check:
%       flgOrigin = -1 and -2: I don't understand the vectors in that case. Could you explain?
%
% 3) Mask - use less memory:
%       flgHalf: Keep the option for rfft.
%       Add an option (like flgHalf = -1) to compute only 1/4 (or 1/8 if 3d) of the grids/vectors. 
%       This could be useful for masks and filters that usually have a C4 symmetry.
%       To regenerate the full grids, we can then have a function that takes this 1/4|8 grid and
%       the ORIGIN/flgShiftOrigin, to compute the entire grid with the desired origin. It should
%       be faster, specially for sphere/cylinder masks and filters where computing the taper can
%       be expensive...

%% Check inputs
if numel(SIZE) == 3
  flg3d = 1;
  ndim = 3;
elseif numel(SIZE) == 2
  flg3d = 0;
  ndim = 2;
else
  error('Only 2D or 3D grid vectors are supported, got ' + str(numel(SIZE)) + 'D');
end
validateattributes(SIZE, {'numeric'}, {'vector', 'numel', ndim, 'nonnegative', 'integer'}, '', 'SIZE');
validateattributes(SHIFT, {'numeric'}, {'vector', 'numel', ndim, 'finite', 'nonnan'}, '', 'SHIFT');

if ~ismember(ORIGIN, [-1, 0, 1, 2])
  error('ORIGIN should be -1, 0, 1, 2, got %s', num2str(ORIGIN))
elseif ~islogical(flgNormalize)
  error('flgNormalize should be a boolean, got %s', class(flgNormalize))
elseif ~islogical(flgHalf)
  error('flgHalf should be a boolean, got %s', class(flgHalf))
end

%% Create vectors with defined origin and shifts.
% For efficiency, shifts are applied to the boundaries and not to the vectors directly.
% On the other hand, normalization is done on the vectors, as rounding errors on the 
% boundaries might lead to significative errors on the vectors.

offset = zeros(1, ndim);
if ORIGIN > 0
    if ORIGIN == 2
        direction = -1;
    else
        direction = 1;
    end
    for dim = 1:ndim
        if ~mod(SIZE(dim), 2)  % even dimensions
            offset(dim) = offset(dim) + direction * 0.5;  % shift half a pixel
        end
    end
end

% centered
if ORIGIN >= 0
    if (flgHalf)
        % there is probably a more elegant way, but I can't bother.
        if ORIGIN == 0
            vX = (rem((SIZE(1)-1)/2, 1)):((SIZE-1)/2);
        else
            vX = SHIFT(1):(ceil((SIZE(1)-1)/2) - SHIFT(1));
        end
    else
        vX = (-SIZE(1)/2 + 0.5 - offset(1) - SHIFT(1)):(SIZE(1)/2 - 0.5 - offset(1) - SHIFT(1));
    end
    vY = (-SIZE(2)/2 + 0.5 - offset(2) - SHIFT(2)):(SIZE(2)/2 - 0.5 - offset(2) - SHIFT(2));
    if (flg3d)
        vZ = (-SIZE(3)/2 + 0.5 - offset(3) - SHIFT(3)):(SIZE(3)/2 - 0.5 - offset(3) - SHIFT(3));
    else
        vZ = nan;
    end

% fft output
else
    if any(SHIFT)
        error('Shifts are not allowed with ORIGIN = -1 (fft output), got %s', mat2str(SHIFT, 2))
    end
    if (flgHalf)
        vX = 0:floor(SIZE(1)/2);
    else
        vX = [0:(floor((SIZE(1)-1)/2)), (-1*floor(SIZE(1)/2)):-1];
    end
    vY = [0:floor((SIZE(2)-1)/2), (-1*floor(SIZE(2)/2)):-1];
    if (flg3d)
        vZ = [0:floor((SIZE(3)-1)/2), (-1*floor(SIZE(3)/2)):-1];
    else
        vZ = nan;
    end
end

if strcmpi(METHOD, 'gpu')
    vX = gpuArray(vX);
    vY = gpuArray(vY);
    if (flg3d); vZ = gpuArray(vZ); end
elseif ~strcmpi(METHOD, 'cpu')
    error("METHOD must be 'gpu' or 'cpu', got %s", METHOD);
end

if (flgNormalize)
    vX = vX ./ SIZE(1);
    vY = vY ./ SIZE(2);
    if (flg3d); vZ = vZ ./ SIZE(3); end
end
end  % EMC_multi_vectorCoordinates
