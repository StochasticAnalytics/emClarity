function IMAGE = EMC_applyBandpass(IMAGE, BANDPASS, OPTION)
%
% IMAGE = EMC_applyBandpass(IMAGE, BANDPASS, OPTION)
% Apply a BANDPASS filter to a real IMAGE.
%
% Input:
%   IMAGE (numeric):           	2d/3d image to Fourier filter.
%                              	Should be real.
%
%   BANDPASS (numeric):         Bandpass filter to apply to the IMAGE.
%                             	NOTE: The bandpass should be zero frequency first (origin=-1)
%                                     as it is directely multiplied by the DFT of the IMAGE.
%                               NOTE: Should have the same precision (single/doule) and same
%                                     method (cpu/gpu) than IMAGE.
%
%   OPTION (cell|struct):       Optional parameters.
%                              	If cell: {field, value; ...}, note the ';' between parameters.
%                              	NOTE: Can be empty.
%                             	NOTE: Unknown fields will raise an error.
%
%     -> 'ifft' (bool):         By default (ifft=true), the filtered IMAGE is switched back to
%                              	real space (ifft) before returning. If false, the IMAGE is left in
%                              	Fourier space (zero frequency first).
%                              	default = true
%
%     -> 'uniform' (bool):      Set the mean to ~0 and variance to ~1 of the output real space IMAGE.
%                             	NOTE: the mean and variance are changed in frequency space for efficiency.
%                             	default = true
%
% Output:
%   IMAGE (single|double):    	Filtered IMAGE; precision and method is unchanged.
%
% Notes:
%   - Once EMC_rffn and EMC_irfftn are ready, add the possibility to use half bandpass. If the size
%     of the bandpass is the size of the expected half grid, trigger the option.
%
% Examples:
%   - img = randn(128,128,128,'single','gpuArray');
%     bandpass = EMC_getBandpass(size(img), 1, nan, 15, 'gpu', {});  % low pass 15A
%     imgFiltered = EMC_applyBandpass(myImg, [10,10,-5,0], {});
%
% Other EMC-files required:
%   EMC_is3d, EMC_sharePrecision, EMC_shareMethod, EMC_getOption
%
% See also EMC_getBandpass.
%

% Created:  18Jan2020, R2019a
% Version:  v.1.0.  symplify the inputs - use EMC_resize directly if you want to apply
%                   a taper and/or pad|crop before applying the bandpass (TF, 2Feb2020).
%           v.1.1.  explicit checks for IMAGE and BANDPASS; unittest (TF, 3Feb2020).
%

%% checkIN
[~, imgSize] = EMC_is3d(size(IMAGE));
if isvector(IMAGE)
    error('EMC:IMAGE', 'IMAGE should be a 2d or 3d matrix, got vector')
elseif ~isequal(size(BANDPASS), imgSize)
    error('EMC:IMAGE', 'IMAGE (size:%s) and BANDPASS (size:%s) should have the same size', ...
          mat2str(imgSize), mat2str(size(BANDPASS)))
elseif ~isreal(IMAGE)
    error('EMC:IMAGE', 'IMAGE should be real, got complex')
elseif ~EMC_sharePrecision(IMAGE, BANDPASS) || ~EMC_shareMethod(IMAGE, BANDPASS)
    error('EMC:IMAGE', 'IMAGE and BANDPASS should have the same precision and method')
end

OPTION = EMC_getOption(OPTION, {'ifft', 'uniform'}, false);

if isfield(OPTION, 'ifft')
    if ~islogical(OPTION.ifft) && ~isscalar(OPTION.ifft)
        error('EMC:ifft', 'OPTION.ifft should be a boolean')
    end
else
    OPTION.ifft = true;  % default
end

if isfield(OPTION, 'uniform')
    if ~islogical(OPTION.uniform) && ~isscalar(OPTION.uniform)
        error('EMC:uniform', 'OPTION.uniform should be a boolean')
    end
else
    OPTION.uniform = true;  % default
end

%% Apply bandpass and uniform if wished.
% IMAGE = EMC_rfftn(IMAGE) .* BANDPASS;
IMAGE = fftn(IMAGE) .* BANDPASS;

if OPTION.uniform
    IMAGE(1) = 0;  % Mean ~0
    IMAGE = IMAGE ./ (sqrt(sum(abs(IMAGE).^2, 'all')) ./ numel(IMAGE));  % Variance ~1
end

if OPTION.ifft
    % IMAGE = EMC_irfftn(IMAGE);
    IMAGE = ifftn(IMAGE, 'symmetric');
end

end  % EMC_applyBandpass
