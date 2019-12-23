function [KERNEL] = EMC_multi_gaussianKernel(SIZE, SIGMA, OPTION)
%
% Compute a 2d or 3d real space gaussian kernel.
%
% SIZE (vector|int):        [x, y] or [x, y, z] dimensions.
%                           Should be odd integers.
%
% SIGMA (float):            Standard deviation of the gaussian.
%
% OPTION (cell|struct):     Optional parameters.
%                           If cell: {field,value ; ...}, note the ';' between parameters.
%                           NOTE: Can be empty.
%                           NOTE: Unknown fields will raise an error.
%
%   -> 'method' (str):      Device to compute the kernel; 'gpu' or 'cpu'
%                           default = 'cpu'
%
%   -> 'precision' (str):   Precision of the kernel; 'single' or 'double'
%                           default = 'single'
%
%   -> 'power' (int):       Power of the gaussian; should be higher than 2.
%                           default = 2
%

%% MAIN
[SIZE, flg3d] = EMC_is3d(SIZE);
validateattributes(SIZE, {'numeric'}, {'integer', 'odd'}, '', 'SIZE')  % more checks in vectorCoordinates
validateattributes(SIGMA, {'numeric'}, {'numel', 1, 'nonnegative'}, '', 'SIGMA')

OPTION = EMC_extract_option(OPTION, {'method', 'precision', 'power'}, false);

if ~isfield(OPTION, 'method')
    OPTION.method = 'cpu';
end
if ~isfield(OPTION, 'precision')
    OPTION.precision = 'single';
end

if isfield(OPTION, 'power')
    pow = OPTION.power;
    if ~isinteger(pow) || pow < 2
        error('power should be higher than 2, got %d', pow)
    else
        error('power should be an integer, got %s', class(pow))
    end
else
    pow = 2;  % default
end

% zeroth order derivative
[vX, vY, vZ] = EMC_multi_vectorCoordinates(SIZE, OPTION.method, {'precision', OPTION.precision});
if (flg3d)
    KERNEL = exp(-0.5 .* (vX'.^pow + vY.^pow + reshape(vZ.^pow, 1, 1, [])) ./ (SIGMA.^pow));
else
    KERNEL = exp(-0.5 .* (vX'.^pow + vY.^pow) ./ (SIGMA.^pow));
end
KERNEL = KERNEL ./ sum(KERNEL(:));  % normalized to have sum of weights to 1.

end
