function [LIMITS] = EMC_multi_limits(CURRENT, DESIRED, OPTION)
% [LIMITS] = EMC_limits(CURRENT, DESIRED, OPTION)
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
% OPTION (cell|struct):     Optional parameters.
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

OPTION = EMC_extract_option(OPTION, {'shift', 'origin'}, false);

if isfield(OPTION, 'origin')
    if ~(OPTION.origin == -1 || OPTION.origin == 1 || OPTION.origin == 2)
        error("origin should be 1, 2, or -1, got %d", OPTION.origin)
    end
else
    OPTION.origin = 1;  % default
end

if isfield(OPTION, 'shift')
    validateattributes(OPTION.shift, {'numeric'}, ...
                       {'row', 'numel', ndim, 'finite', 'nonnan'}, '', 'shift');
else
    OPTION.shift = zeros(1, ndim);  % default
end

%% MAIN
size_difference = DESIRED - CURRENT;
remain = mod(size_difference, 2);
odd_input = mod(CURRENT, 2);

LIMITS = zeros(1, ndim .* 2);
left = floor((size_difference ./ 2)) + odd_input .* remain + OPTION.shift;
right = floor((size_difference ./ 2)) + ~(odd_input) .* remain - OPTION.shift;

if OPTION.origin == 1 || OPTION.origin == -1
    LIMITS(1:2:end) = left;
    LIMITS(2:2:end) = right;
elseif OPTION.origin == 2
    LIMITS(1:2:end) = right;
    LIMITS(2:2:end) = left;
else
    error('origin should be -1, 1 or 2, got %d', OPTION.origin)
end

end  % EMC_limits
