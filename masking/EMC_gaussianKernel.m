function KERNEL = EMC_gaussianKernel(SIZE, SIGMA, OPTION)
%
% [KERNEL] = EMC_gaussianKernel(SIZE, SIGMA, OPTION)
% Compute a 2d/3d real/fft gaussian kernel (zeroth order).
%
% Input:
%   SIZE (int vector):    	Size (in pixel) of the kernel to compute; [x, y, z] or [x, y].
%
%   SIGMA (float|vector):   Standard deviation of the gaussian.
%                           If float: isotropic kernel
%                           If vector: anisotropic kernel (one sigma per axis); should correspond to SIZE.
%
%   OPTION (cell|struct):  	Optional parameters.
%                           If cell: {field, value; ...}, note the ';' between parameters.
%                           NOTE: Can be empty.
%                           NOTE: Unknown fields will raise an error.
%
%     -> 'origin' (bool): 	Origin convention; -1, 0, 1 or 2.
%                           If -1: switch to a reciprocal kernel. See 'Notes' for more details.
%                           default = 1
%
%	  -> 'method' (str):  	Device to compute the kernel; 'gpu' or 'cpu'
%                           default = 'cpu'
%
% 	  -> 'precision' (str):	Precision of the kernel; 'single' or 'double'
%                           default = 'single'
%
% Output:
%   KERNEL (num):           Gaussian kernel of SIZE.
%
% Note:
%   - If origin=-1, the output kernel is meant to be multiplied directly to the fft of the image to blur
%     and not convolved in real space. Even though it gives similar results compared with real space
%     convolution, it doesn't give exactly the same results... To be continued...
%
% Other EMC-files required:
%   EMC_is3d.m, EMC_getOption.m, EMC_coordVectors.m
%
% See also EMC_coordVectors
%

% Created:  18Jan2020, R2019a
% Version:  v.1.1   fft kernel (origin=-1) can use more than one sigma; 1 per axis (TF, 27Jan2020).
%

%% checkIN
[is3d, SIZE, ndim] = EMC_is3d(SIZE);

if ~isnumeric(SIGMA)
    error('EMC_gaussianKernel:SIGMA', 'SIGMA should be numeric, got %s', class(SIGMA))
elseif isscalar(SIGMA) && SIGMA > 0
    SIGMA = zeros(1,ndim) + SIGMA;  % isotropic
elseif ~(isvector(SIGMA) && all(SIGMA > 0) && length(SIGMA) == ndim)  % anisotropic
    error('EMC_gaussianKernel:SIGMA', 'SIGMA should be a nonnegative float or a vector of size %d', ndim)
end

OPTION = EMC_getOption(OPTION, {'method', 'precision', 'origin'}, false);

% checks in EMC_coordVectors
if ~isfield(OPTION, 'method')
    OPTION.method = 'cpu';
end
if ~isfield(OPTION, 'precision')
    OPTION.precision = 'single';
end

% origin
if isfield(OPTION, 'origin')
    if ~isnumeric(OPTION.origin) || ~isscalar(OPTION.origin)
        error('EMC_gaussianKernel:origin', 'origin should be an integer')
    elseif OPTION.origin == 0 || OPTION.origin == 1 || OPTION.origin == 2
        isfft = false;
        normalizeCoords = false;
    elseif OPTION.origin == -1
        isfft = true;
        normalizeCoords = true;
    else
        error('EMC_gaussianKernel:origin', "origin should be -1, 0, 1, or 2, got %.02f", OPTION.origin)
    end
else
    OPTION.origin = 1;  % default
    isfft = false;
    normalizeCoords = false;
end

%% zeroth order derivative
[vX, vY, vZ] = EMC_coordVectors(SIZE, OPTION.method, {'precision', OPTION.precision;
                                                      'origin', OPTION.origin;
                                                      'normalize', normalizeCoords});

if isfft
    if is3d
        KERNEL = exp(-2 * (SIGMA(1)*pi.*vX').^2 + (SIGMA(2)*pi.*vY).^2 + reshape((SIGMA(3)*pi.*vZ).^2,1,1,[]));
    else
        KERNEL = exp(-2 * (SIGMA(1)*pi.*vX').^2 + (SIGMA(2)*pi.*vY).^2);
    end
else
    if is3d
        KERNEL = exp(-0.5 .* ((vX./SIGMA(1))'.^2 + (vY./SIGMA(2)).^2 + reshape((vZ./SIGMA(3)).^2,1,1,[])));
    else
        KERNEL = exp(-0.5 .* ((vX./SIGMA(1))'.^2 + (vY./SIGMA(2)).^2));
    end
    KERNEL = KERNEL ./ sum(KERNEL(:));  % normalized to have sum of weights to 1.
end

end
