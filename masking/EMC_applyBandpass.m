function IMAGE = EMC_applyBandpass(IMAGE, BANDPASS, OPTION)
%
% IMAGE = EMC_applyBandpass(IMAGE, BANDPASS, OPTION)
% Apply a BANDPASS filter to a real IMAGE.
%
% Input:
%   IMAGE (numeric):           	2d/3d real image to Fourier filter.
%
%   BANDPASS (numeric):         Bandpass filter to apply to the IMAGE.
%                               Half bandpass are also accepted.
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
%     -> 'standardize' (bool):  Set the mean to ~0 and variance to ~1 of the output real space IMAGE.
%                             	NOTE: the mean and variance are changed in frequency space.
%                             	default = true
%
% Output:
%   IMAGE (single|double):    	Filtered IMAGE; precision and method is unchanged.
%
% Examples:
%   - img = randn(128,128,128,'single','gpuArray');
%     bandpass = EMC_getBandpass(size(img), 1, nan, 15, 'gpu', {});  % low pass 15A
%     imgFiltered = EMC_applyBandpass(myImg, [10,10,-5,0], {});
%
% Other EMC-files required:
%   EMC_is3d, EMC_sharePrecision, EMC_shareMethod, EMC_getOption, EMC_rfftn, EMC_irfftn
%
% See also EMC_getBandpass.
%

% Created:  18Jan2020, R2019a
% Version:  v.1.0.  symplify the inputs - use EMC_resize directly if you want to apply
%                   a taper and/or pad|crop before applying the bandpass (TF, 2Feb2020).
%           v.1.1.  explicit checks for IMAGE and BANDPASS; unittest (TF, 3Feb2020).
%           v.1.2.  half bandpass are now supported (TF, 8Feb2020).
%           v.1.2.1 bug fix: standardization with half grid is now properly done. (TF, 22Feb2020).
%

%% checkIN
[~, imgSize] = EMC_is3d(size(IMAGE));
if isvector(IMAGE)
    error('EMC:IMAGE', 'IMAGE should be a 2d or 3d matrix, got vector')
elseif ~isreal(IMAGE)
    error('EMC:IMAGE', 'IMAGE should be real, got complex')
elseif ~EMC_sharePrecision(IMAGE, BANDPASS) || ~EMC_shareMethod(IMAGE, BANDPASS)
    error('EMC:IMAGE', 'IMAGE and BANDPASS should have the same precision and method')
end

OPTION = EMC_getOption(OPTION, {'ifft', 'standardize'}, false);

if isfield(OPTION, 'ifft')
    if ~islogical(OPTION.ifft) && ~isscalar(OPTION.ifft)
        error('EMC:ifft', 'OPTION.ifft should be a boolean')
    end
else
    OPTION.ifft = true;  % default
end

if isfield(OPTION, 'standardize')
    if ~islogical(OPTION.standardize) && ~isscalar(OPTION.standardize)
        error('EMC:standardize', 'OPTION.standardize should be a boolean')
    end
else
    OPTION.standardize = true;  % default
end

%% Apply bandpass and uniform if wished.
bandSize = size(BANDPASS);

% full bandpass
if isequal(bandSize, imgSize)
    IMAGE = fftn(IMAGE) .* BANDPASS;

    if OPTION.standardize
        IMAGE(1) = 0;
        IMAGE = IMAGE ./ (sqrt(sum(abs(IMAGE).^2, 'all')) / numel(IMAGE));
    end

    if OPTION.ifft
        IMAGE = ifftn(IMAGE, 'symmetric');
    end

% half bandpass
elseif isequal([floor(imgSize(1)/2)+1, imgSize(2:end)], bandSize)
    IMAGE = EMC_rfftn(IMAGE) .* BANDPASS;

    if OPTION.standardize
        IMAGE(1) = 0;

        % Capture every independant chunk (same for 2d or 3d)
        cD = ceil(imgSize(1)/2);  % center donor
        factor = sum(abs(IMAGE(1,:,:)).^2, 'all');  % unique row/plane
        factor = factor + 2*sum(abs(IMAGE(2:cD,:,:)).^2, 'all');  % common chunk
        if ~mod(imgSize(1),2); factor = factor + sum(abs(IMAGE(cD+1,:,:)).^2, 'all'); end  % unique row/plane

        IMAGE = IMAGE ./ (sqrt(factor) / prod(imgSize, 'native'));
    end

    if OPTION.ifft
        IMAGE = EMC_irfftn(IMAGE, imgSize);
    end
else
    error('EMC:IMAGE', 'IMAGE (size:%s) and BANDPASS (size:%s) should have the same size', ...
          mat2str(imgSize), mat2str(size(BANDPASS)))
end

end  % EMC_applyBandpass
