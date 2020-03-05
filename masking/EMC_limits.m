function [LIMITS] = EMC_limits(CURRENT, DESIRED, OPTION)
%
% [LIMITS] = EMC_limits(CURRENT, DESIRED, OPTION)
% Compute the LIMITS: pixels to add/remove to the CURRENT size to
% resize it to the DESIRED size while preserving the origin.
%
% Input:
%   CURRENT (vector):           Current size in pixel; [x, y, z] or [x, y]
%                               NOTE: the number of dimensions is actually not limited to 3.
%
%   DESIRED (vector):           Desired size in pixel; [x, y, z] or [x, y]
%                               Should correspond to CURRENT.
%
%   OPTION (cell|struct):       Optional parameters.
%                               If cell: {field, value; ...}, note the ';' between parameters.
%                               NOTE: Can be empty.
%                               NOTE: Unknown fields will raise an error.
%
%     -> 'origin' (int):        Origin convention - Center of rotation.
%                               NOTE: origin=0 is not allowed. See EMC_resize for more details.
%                               defaut = 1
%
%     -> 'shift' (int vector):  [x, y] or [x, y, z] shifts to apply (in pixel).
%                               Should correspond to CURRENT.
%                               WARNING: As the LIMITS are meant to be used with EMC_resize and
%                                        to be consistent with EMC_coordVectors (which defines most
%                                        coordinates), positive shifts translate the origin to the
%                                        LEFT which has the effect to shift the output image of
%                                        EMC_resize to the right. Negative shifts translate the origin
%                                        to the RIGHT.
%                               NOTE: shifts are not allowed with origin=-1.
%                               default = no shifts
%
% Outputs:
%   LIMITS (row vector):        Computed limits. Correspond to EMC_resize:LIMITS (see example below).
%                               numel: numel(CURRENT) * 2.
%                               Format - if 3d: [z_left, z_right, y_left, y_right, x_left, x_right]
%                                        if 2d: [y_left, y_right, x_left, x_right]
% Note:
%   - This function operates in pixel space (integer), therefore all inputs must be in pixel space
%     and 'origin' = 0 is not allowed.
%
% Other EMC-files required:
%   EMC_getOption.m
%
% Example:
%   - Shift an image by 5 pixel in x to the right:
%     img = randn(128,128);
%     limits = EMC_limits(size(img), size(img), {'shift', [5,0]});
%     imgShifted = EMC_resize(img, limits, {});
%
%   - Pad img to the desired size:
%     img = randn(128,128);
%     limits = EMC_limits(size(img), [150,155], {});
%     imgPadded = EMC_resize(img, limits, {}); 
%
% See also EMC_resize
%

% Created:  18Jan2020
% Version:  v.1.0.1 more explicit checkIN
%           v.1.1   force shifts to be integers, correct shifts with origin=2
%                   and unittest (TF, 31Jan2020).
%

%% checkIN
if ~isnumeric(CURRENT) || ~isrow(CURRENT) || any(isinf(CURRENT)) ||  ~all(CURRENT > 0) || any(rem(CURRENT, 1)) 
    error('EMC:LIMITS', 'CURRENT should be a row vector of positive integers')
elseif ~isnumeric(DESIRED) || ~isrow(DESIRED) || any(isinf(DESIRED)) ||  ~all(DESIRED > 0) || any(rem(DESIRED, 1)) 
    error('EMC:LIMITS', 'DESIRED should be a row vector of positive integers')
elseif ~isequal(size(CURRENT), size(DESIRED))
    error('EMC:LIMITS', 'CURRENT and DESIRED should have the same size, got %s ~= %s', ...
          mat2str(size(CURRENT)), mat2str(size(DESIRED)))
end

OPTION = EMC_getOption(OPTION, {'shift', 'origin'}, false);

if isfield(OPTION, 'origin')
    if ~isscalar(OPTION.origin) || ~isnumeric(OPTION.origin) || ...
       ~(OPTION.origin == -1 || OPTION.origin == 1 || OPTION.origin == 2)
        error('EMC:origin', "OPTION.origin should be 1, 2, or -1")
    end
else
    OPTION.origin = 1;  % default
end

if isfield(OPTION, 'shift')
    if ~isnumeric(OPTION.shift) || ~isrow(OPTION.shift)
        error('EMC:shift', 'OPTION.shift should be a vector of float|int, got %s', class(OPTION.shift))
    elseif any(isnan(OPTION.shift)) || any(isinf(OPTION.shift)) || any(rem(OPTION.shift, 1))
        error('EMC:shift', 'OPTION.shift should only contain integers, got %s', mat2str(OPTION.shift, 2))
    elseif ~isequal(size(OPTION.shift), size(CURRENT))
        error('EMC:shift', 'CURRENT and OPTION.shift should have the same size, got %s - %s', ...
              mat2str(size(CURRENT)), mat2str(size(OPTION.shift)))
    elseif OPTION.origin == -1 && any(OPTION.shift)
      	error('EMC:shift', 'OPTION.shifts are not allowed with half=true or origin=-1, got %s', ...
              mat2str(OPTION.shift, 2))
    end
else
    OPTION.shift = zeros(1, numel(CURRENT));  % default
end

%% MAIN
sizeDiff = DESIRED - CURRENT;
remain = mod(sizeDiff, 2);
isOddInput = mod(CURRENT, 2);

LIMITS = zeros(1, numel(CURRENT) .* 2);
left = floor((sizeDiff ./ 2)) + isOddInput .* remain;
right = floor((sizeDiff ./ 2)) + ~(isOddInput) .* remain;

if OPTION.origin == 1 || OPTION.origin == -1
    LIMITS(1:2:end) = left + OPTION.shift;
    LIMITS(2:2:end) = right - OPTION.shift;
else  % OPTION.origin == 2
    LIMITS(1:2:end) = right + OPTION.shift;
    LIMITS(2:2:end) = left - OPTION.shift;
end

end  % EMC_limits
