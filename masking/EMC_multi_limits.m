function [LIMITS] = EMC_multi_limits(CURRENT, DESIRED, OPTIONAL)
% [LIMITS] = EMC_limits(CURRENT, DESIRED, OPTIONAL)
%
% Compute the LIMITS: pixels to add/remove to the CURRENT size to
% resize it to the DESIRED size while preserving the origin.
%
% CURRENT (vector):         Current size; [x, y(, z)]
%                           NOTE: the number of dimensions is actually not limited to 3.
%
% DESIRED (vector):         Desired size; [x, y(, z)]
%                           Should correspond to CURRENT.
%
% OPTIONAL (cell|struct):   Optional parameters.
%                           If cell: {field, value; ...}, note the ';' between parameters.
%                           NOTE: Can be empty.
%                           NOTE: Unknown fields will raise an error.
%
%   -> 'origin' (int):      Origin convention - Center of rotation.
%                           -1, 1 or 2; see EMC_multi_gridVectors for more details.
%                           defaut = 1
%
%   -> 'shift' (vector):    [x, y] or [x, y, z] translations to apply.
%                           Should correspond to CURRENT.
%                           NOTE: Shifts the CURRENT origin (defined by the 'origin' parameter)
%                                 by this amount of pixels before calculating the LIMITS.
%                           default = no shifts
%
%-------
% NOTE:                     This function operates in pixel space (integer), therefore all inputs
%                           must be in pixel space and 'origin' = 0 is not allowed.
%
%-------
% RETURN:                   LIMITS
%                           Format is as follows:
%                           [z_left, z_right, y_left, y_right, x_left, x_right] if 3D; shape:(6,)
%                           [y_left, y_right, x_left, x_right] if 2D; shape:(4,)
% 
% See also EMC_resize.m

%% checkIN
ndim = numel(CURRENT);
validateattributes(CURRENT, {'numeric'}, {'vector', 'integer', 'positive'}, '', 'CURRENT');
validateattributes(DESIRED, {'numeric'}, {'vector', 'integer', 'positive', 'numel', ndim}, '', 'DESIRED');

OPTIONAL = EMC_extract_optional(OPTIONAL, {'shift', 'origin'});

if isfield(OPTIONAL, 'origin')
    if ~(OPTIONAL.origin == -1 || OPTIONAL.origin == 1 || OPTIONAL.origin == 2)
        error("origin should be 1, 2, or -1, got %d", OPTIONAL.origin)
    end
else
    OPTIONAL.origin = 1;  % default
end

if isfield(OPTIONAL, 'shift')
    validateattributes(OPTIONAL.shift, {'numeric'}, ...
                       {'row', 'numel', ndim, 'finite', 'nonnan'}, '', 'shift');
else
    OPTIONAL.shift = zeros(1, ndim);  % default
end

%% MAIN
size_difference = DESIRED - CURRENT;
remain = mod(size_difference, 2);
odd_input = mod(CURRENT, 2);

LIMITS = zeros(1, ndim .* 2);
left = floor((size_difference ./ 2)) + odd_input .* remain + OPTIONAL.shift;
right = floor((size_difference ./ 2)) + ~(odd_input) .* remain - OPTIONAL.shift;

if OPTIONAL.origin == 1 || OPTIONAL.origin == -1
    LIMITS(1:2:end) = left;
    LIMITS(2:2:end) = right;
elseif OPTIONAL.origin == 2
    LIMITS(1:2:end) = right;
    LIMITS(2:2:end) = left;
else
    error('origin should be -1, 1 or 2, got %d', OPTIONAL.origin)
end

end  % EMC_limits
