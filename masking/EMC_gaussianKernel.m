function KERNEL = EMC_gaussianKernel(SIZE, SIGMA, METHOD, OPTION)
%
% [KERNEL] = EMC_gaussianKernel(SIZE, SIGMA, METHOD, OPTION)
% Compute a 2d/3d real gaussian kernel (zeroth order).
%
% Input:
%   SIZE (int vector):      Size (in pixel) of the kernel to compute; [x, y, z] or [x, y].
%                           NOTE: [1,1] is not accepted.
%
%   SIGMA (float|vector):   Standard deviation of the gaussian.
%                           Should be positive.
%                           If float: isotropic kernel
%                           If vector: anisotropic kernel (one sigma per axis); should correspond to SIZE.
%
%   METHOD (str):           Device to compute the kernel; 'gpu' or 'cpu'
%
%   OPTION (cell|struct):   Optional parameters.
%                           If cell: {field, value; ...}, note the ';' between parameters.
%                           NOTE: Can be empty.
%                           NOTE: Unknown fields will raise an error.
%
%     -> 'precision' (str): Precision of the kernel; 'single' or 'double'
%                           default = 'single'
%
% Output:
%   KERNEL (num):           Gaussian kernel of SIZE.
%
% Other EMC-files required:
%   EMC_is3d.m, EMC_getOption.m, EMC_coordVectors.m
%
% See also EMC_convnSeparable
%

% Created:  18Jan2020, R2019a
% Version:  v.1.0.  switch from sum(kernel(:)) to sum(kernel, 'all').
%           v.1.1.0 unittest (TF, 4Feb2020).
%           v.1.1.1 odd SIZEs are now supported (TF, 4Feb2020).
%           v.1.1.2 returning a scalar is not possible now (TF, 4Feb2020).
%           v.1.2   METHOD is now a required input (TF, 9Mar2020).
%

%% checkIN
[is3d, SIZE, ndim] = EMC_is3d(SIZE);
if all(SIZE == 1)
    error('EMC:SIZE', 'The output kernel cannot be a scalar (got SIZE: %s)', mat2str(SIZE))
end

if ~isnumeric(SIGMA) || any(isinf(SIGMA)) || ~all(SIGMA > 0)
    error('EMC:SIGMA', 'SIGMA should be a positive numeric, without NaN/Inf elements, got %s', class(SIGMA))
elseif isscalar(SIGMA)
    SIGMA = zeros(1,ndim) + SIGMA;  % isotropic
elseif ~isvector(SIGMA) || length(SIGMA) ~= ndim  % anisotropic
    error('EMC:SIGMA', 'SIGMA should be a positive float or a vector of size %d', ndim)
end

if ~(ischar(METHOD) || isstring(METHOD)) || ~(strcmpi(METHOD, 'gpu') || strcmpi(METHOD, 'cpu'))
    error('EMC:METHOD', "METHOD should be 'gpu' or 'cpu'")
end

OPTION = EMC_getOption(OPTION, {'precision'}, false);

% checks in EMC_coordVectors
if ~isfield(OPTION, 'precision')
    OPTION.precision = 'single';
end

%% zeroth order derivative
% Use the real center (useful if even dimension).
[vX, vY, vZ] = EMC_coordVectors(SIZE, METHOD, {'precision', OPTION.precision; 'origin', 0});

if is3d
    KERNEL = exp(-0.5 .* ((vX'./SIGMA(1)).^2 + (vY./SIGMA(2)).^2 + reshape((vZ./SIGMA(3)).^2,1,1,[])));
else
    if SIZE(1) == 1  % row vector
        KERNEL = exp(-0.5 .* (vY./SIGMA(2)).^2);
    elseif SIZE(2) == 1
        KERNEL = exp(-0.5 .* (vX'./SIGMA(1)).^2);
    else
        KERNEL = exp(-0.5 .* ((vX'./SIGMA(1)).^2 + (vY./SIGMA(2)).^2));
    end
end
KERNEL = KERNEL ./ sum(KERNEL, 'all');  % normalized to have sum of weights to 1.

end
