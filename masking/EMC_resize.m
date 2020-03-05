function [OUT] = EMC_resize(IMAGE, LIMITS, OPTION)
%
% [OUT] = EMC_resize(IMAGE, LIMITS, OPTION)
% Pad and/or crop an IMAGE.
% 
% Inputs:
%   IMAGE (single|double):          2d/3d IMAGE to pad and/or crop.
%
%   LIMITS (row vector):            Number of pixels to pad or crop, for each axis.
%                                   If 2d IMAGE: [xleft, xright, yleft, yright].
%                                   If 3d IMAGE: [xleft, xright, yleft, yright, zleft, zright].
%                                   See EMC_limits for more details.
%                                   NOTE: Positive values indicate the number of pixels to pad, while
%                                         negative values indicate the number of pixels to crop.
%
%   OPTION (cell|struct):           Optional parameters.
%                                   If cell: {param, value ; ...}, note the ';' between parameters.
%                                   NOTE: Can be empty.
%                                   NOTE: Unknown fields will raise an error.
%
%     -> 'origin' (int):            Origin convention - Center of rotation.
%                                   -1, 1 or 2; see EMC_multi_gridVectors for more details.
%                                   NOTE: 1|2 produce identical results because LIMITS is not relative
%                                         to the origin but to the edges of the IMAGE. On the other hand,
%                                         'origin'=-1 is different as it specifies a not-centered IMAGE
%                                         (zero frequency first).
%                                   NOTE: This function opperates in pixel space: origin=0 is not allowed.
%                                   defaut = 1
%
%     -> 'value' (float|int|str):   If float|int: value to pad with. Casted to desired 'precision'.
%                                   If 'uniform': pad with white gaussian noise.
%                                   If 'mean': pad with the mean of the IMAGE.
%                                   default = 0
%
%     -> 'taper' (bool|cell|vector|float|int):
%                                   Apply a taper to the IMAGE before padding.
%                                   If bool: Apply or not the default taper; equivalent to {'cosine', 7}.
%                                   If cell: {type, size} with type = 'linear' or 'cosine'
%                                            and with size = size of the taper (in pixel).
%                                   If vector: column vector used as taper (left to right <-> center to edge).
%                                   If float|int: equivalent to [float]; one pixel taper.
%                                   NOTE: Every IMAGE dimension should be larger than the taper, with the
%                                         exception of origin=-1, which requires the dimensions to be at
%                                         least 2 times larger than the taper.
%                                   default = {'cosine', 7}
%
%     -> 'force_taper' (bool):      By default, (force_taper=false) only the edges that are padded are
%                                   tapered. If true, apply the taper to every edges even if they are not
%                                   padded. If cropping is required, apply the taper AFTER cropping.
%                                   NOTE: this has no effect if 'taper' = false.
%                                   default = false
%
%     -> 'precision' (str):         Precision of the padded/cropped IMAGE; 'single' or 'double'.
%                                   default = same as IMAGE
%
% Output:
%   OUT (single|double):            Output 2d/3d image.
%
% Examples:
%   - [OUT] = EMC_resize(randn(64,64), [10,10,-5,0], {})
%   - [OUT] = EMC_resize(randn(64,64), [10,10,-5,0], {'value', 2; 'origin', -1})
%
% Other EMC-files required:
%   EMC_setPrecision.m, EMC_getClass.m, EMC_is3d.m, EMC_getOption.m, EMC_isOnGpu, EMC_taper
%
% See also EMC_limits.
%

% Created:  18Jan2020, R2019a
% Version:  v.1.0.1 better error message display (TF, 20Jan2020).
%           v.1.1   unittest (TF, 21Jan2020).
%           v.1.1.1 add isscalar when expecting a scalar (TF, 21Jan2020).
%           v.1.1.2 switch from mean(a(:)) to mean('all'), and std(a, 0, 'all') (TF, 1Feb2020).
%

%% MAIN
[OPTION, flg] = checkIN(IMAGE, LIMITS, OPTION);

% Short cut: nothing to do to the IMAGE, except maybe changing its precision.
if ~flg.pad && ~flg.crop && ~(flg.taper && OPTION.force_taper)
    OUT = EMC_setPrecision(IMAGE, OPTION.precision);
    return
end

inSize = size(IMAGE);
outSize = inSize + sum(reshape(LIMITS, 2, []));
if any(outSize < 1)
    outSize = size(IMAGE) + sum(reshape(LIMITS, 2, []));
    error('EMC:LIMITS', 'Output dimensions should be at least 1 pixel, got output size: %s', mat2str(outSize))
end

% Option to pad with white gaussian noise.
if flg.uniform
    val = EMC_setPrecision(mean(IMAGE, 'all'), OPTION.precision);
    std_ = EMC_setPrecision(std(IMAGE, 0, 'all'), OPTION.precision);
else
    val = OPTION.value;
    std_ = nan;
end

% Allocate memory for the output image.
% At this point, if nothing to pad or crop, it means the goal
% is only to apply a taper and maybe also change the precision.
% In this case, no need to allocate new memory.
if flg.pad || flg.crop
    if flg.gpu
        if flg.uniform
            OUT = randn(outSize, OPTION.precision, 'gpuArray') .* std_ + val;
        else
            OUT = zeros(outSize, OPTION.precision, 'gpuArray');
            if OPTION.value ~= 0
                OUT = OUT + OPTION.value;
            end
        end
    else
        if flg.uniform
            OUT = randn(outSize, OPTION.precision) .* std_ + val;
        else
            OUT = zeros(outSize, OPTION.precision);
            if OPTION.value ~= 0
                OUT = OUT + OPTION.value;
            end
        end
    end
end

x=1; y=2; z=3;
if flg.is3d
    if flg.fft  % not-centered IMAGE: 'origin' = -1
        crop = reshape(LIMITS .* (LIMITS < 0), 2, []);
        l = floor((inSize+1) / 2) + crop(1, :);  % left side
        r = floor((inSize-2) / 2) + crop(2, :);  % right side
        
        % Taper
        if flg.taper; IMAGE = applyTaper_fft3d(IMAGE, LIMITS, OPTION, l, r, val); end
        
        % Pad/Crop
        if flg.pad || flg.crop
            OUT(1:l(x),       1:l(y),       1:l(z))       = IMAGE(1:l(x),       1:l(y),       1:l(z));
            OUT(end-r(x):end, 1:l(y),       1:l(z))       = IMAGE(end-r(x):end, 1:l(y),       1:l(z));
            OUT(1:l(x),       end-r(y):end, 1:l(z))       = IMAGE(1:l(x),       end-r(y):end, 1:l(z));
            OUT(end-r(x):end, end-r(y):end, 1:l(z))       = IMAGE(end-r(x):end, end-r(y):end, 1:l(z));
            OUT(1:l(x),       1:l(y),       end-r(z):end) = IMAGE(1:l(x),       1:l(y),       end-r(z):end);
            OUT(end-r(x):end, 1:l(y),       end-r(z):end) = IMAGE(end-r(x):end, 1:l(y),       end-r(z):end);
            OUT(1:l(x),       end-r(y):end, end-r(z):end) = IMAGE(1:l(x),       end-r(y):end, end-r(z):end);
            OUT(end-r(x):end, end-r(y):end, end-r(z):end) = IMAGE(end-r(x):end, end-r(y):end, end-r(z):end);
        else
            OUT = EMC_setPrecision(IMAGE, OPTION.precision);  % force_taper without pad or crop
        end
    else  % centered IMAGE: 'origin' = 0|1|2
        crop = abs(LIMITS .* (LIMITS < 0));
        pad  = LIMITS .* (LIMITS > 0);

        % Taper
        if flg.taper && (any(pad) || OPTION.force_taper)
            IMAGE = applyTaper_real3d(IMAGE, LIMITS, OPTION, pad, crop, val);
        end

        % Pad/Crop
        if flg.pad && flg.crop
            OUT(1+pad(1):end-pad(2), ...
                1+pad(3):end-pad(4), ...
                1+pad(5):end-pad(6)) = IMAGE(1+crop(1):end-crop(2), ...
                                                   1+crop(3):end-crop(4), ...
                                                   1+crop(5):end-crop(6));
        elseif flg.pad
            OUT(1+pad(1):end-pad(2), 1+pad(3):end-pad(4), 1+pad(5):end-pad(6)) = IMAGE;
        elseif flg.crop
            OUT(:,:,:) = IMAGE(1+crop(1):end-crop(2), 1+crop(3):end-crop(4), 1+crop(5):end-crop(6));
        else
            OUT = EMC_setPrecision(IMAGE, OPTION.precision);  % force_taper without pad or crop
        end
    end
else  % 2d
    if flg.fft  % not-centered IMAGE: 'origin' = -1
        crop = reshape(LIMITS .* (LIMITS < 0), 2, []);
        l = floor((inSize+1) / 2) + crop(1, :);  % left side
        r = floor((inSize-2) / 2) + crop(2, :);  % right side

        % Taper
        if flg.taper; IMAGE = applyTaper_fft2d(IMAGE, LIMITS, OPTION, l, r, val); end

        % Pad/Crop
        if flg.pad || flg.crop
            OUT(1:l(x),       1:l(y))       = IMAGE(1:l(x),       1:l(y));
            OUT(end-r(x):end, 1:l(y))       = IMAGE(end-r(x):end, 1:l(y));
            OUT(1:l(x),       end-r(y):end) = IMAGE(1:l(x),       end-r(y):end);
            OUT(end-r(x):end, end-r(y):end) = IMAGE(end-r(x):end, end-r(y):end);
        else
            OUT = EMC_setPrecision(IMAGE, OPTION.precision);  % force_taper without pad or crop
        end
    else  % centered IMAGE: 'origin' = 0|1|2
        crop = abs(LIMITS .* (LIMITS < 0));
        pad  = LIMITS .* (LIMITS > 0);

        % Taper
        if flg.taper && (any(pad) || OPTION.force_taper)
            IMAGE = applyTaper_real2d(IMAGE, LIMITS, OPTION, pad, crop, val);
        end

        % Pad/Crop
        if flg.pad && flg.crop
            OUT(1+pad(1):end-pad(2), 1+pad(3):end-pad(4)) = IMAGE(1+crop(1):end-crop(2), 1+crop(3):end-crop(4));
        elseif flg.pad
            OUT(1+pad(1):end-pad(2), 1+pad(3):end-pad(4)) = IMAGE;
        elseif flg.crop
            OUT(:,:) = IMAGE(1+crop(1):end-crop(2), 1+crop(3):end-crop(4));
        else
            OUT = EMC_setPrecision(IMAGE, OPTION.precision);  % force_taper without pad or crop
        end
    end
end

end  % EMC_resize


function IMAGE = applyTaper_fft3d(IMAGE, LIMITS, OPTION, l, r, val)

s = numel(OPTION.taper);
pad = LIMITS > 0;
try
    if any(pad(1:2) > 0) || OPTION.force_taper
        t = OPTION.taper';
        tf = flip(t);
        IMAGE(l(1)-s+1:l(1), :, :)         = IMAGE(l(1)-s+1:l(1), :, :)         .*t  +val.*(1-t);
        IMAGE(end-r(1):end-r(1)+s-1, :, :) = IMAGE(end-r(1):end-r(1)+s-1, :, :) .*tf +val.*(1-tf);
    end
    if any(pad(3:4) > 0) || OPTION.force_taper
        t = OPTION.taper;
        tf = flip(t);
        IMAGE(:, l(2)-s+1:l(2), :)         = IMAGE(:, l(2)-s+1:l(2), :)         .*t  +val.*(1-t);
        IMAGE(:, end-r(2):end-r(2)+s-1, :) = IMAGE(:, end-r(2):end-r(2)+s-1, :) .*tf +val.*(1-tf);
    end
    if any(pad(5:6) > 0) || OPTION.force_taper
        t = reshape(OPTION.taper, 1, 1, []);
        tf = flip(t);
        IMAGE(:, :, l(3)-s+1:l(3))         = IMAGE(:, :, l(3)-s+1:l(3))         .*t  +val.*(1-t);
        IMAGE(:, :, end-r(3):end-r(3)+s-1) = IMAGE(:, :, end-r(3):end-r(3)+s-1) .*tf +val.*(1-tf);
    end
catch
    raiseError(IMAGE, LIMITS, OPTION);
end

end  % applyTaper_fft3d


function IMAGE = applyTaper_fft2d(IMAGE, LIMITS, OPTION, l, r, val)

s = numel(OPTION.taper);
pad = LIMITS > 0;
try
    if any(pad(1:2) > 0) || OPTION.force_taper
        t = OPTION.taper';
        tf = flip(t);
        IMAGE(l(1)-s+1:l(1), :)         = IMAGE(l(1)-s+1:l(1), :)         .*t  +val.*(1-t);
        IMAGE(end-r(1):end-r(1)+s-1, :) = IMAGE(end-r(1):end-r(1)+s-1, :) .*tf +val.*(1-tf);
    end
    if any(pad(3:4) > 0) || OPTION.force_taper
        t = OPTION.taper;
        tf = flip(t);
        IMAGE(:, l(2)-s+1:l(2))         = IMAGE(:, l(2)-s+1:l(2))         .*t  +val.*(1-t);
        IMAGE(:, end-r(2):end-r(2)+s-1) = IMAGE(:, end-r(2):end-r(2)+s-1) .*tf +val.*(1-tf);
    end
catch
    raiseError(IMAGE, LIMITS, OPTION);
end

end  % applyTaper_fft2d


function IMAGE = applyTaper_real3d(IMAGE, LIMITS, OPTION, pad, crop, val)

l = numel(OPTION.taper) + crop;  % broadcast
ty = OPTION.taper;
tx = ty';
tz = reshape(ty, 1, 1, []);

try
    % top, down, left, right, bottom, front
    if pad(1) || OPTION.force_taper
        IMAGE(1+crop(1):l(1),:,:)         = IMAGE(1+crop(1):l(1),:,:) .* flip(tx) + val.*(1-flip(tx));
    end
    if pad(2) || OPTION.force_taper
        IMAGE(end-l(2)+1:end-crop(2),:,:) = IMAGE(end-l(2)+1:end-crop(2),:,:) .* tx + val.*(1-tx);
    end
    if pad(3) || OPTION.force_taper
        IMAGE(:,1+crop(3):l(3),:)         = IMAGE(:,1+crop(3):l(3),:) .* flip(ty) + val.*(1-flip(ty));
    end
    if pad(4) || OPTION.force_taper
        IMAGE(:,end-l(4)+1:end-crop(4),:) = IMAGE(:,end-l(4)+1:end-crop(4),:) .* ty + val.*(1-ty);
    end
    if pad(5) || OPTION.force_taper
        IMAGE(:,:,1+crop(5):l(5))         = IMAGE(:,:,1+crop(5):l(5)) .* flip(tz) + val.*(1-flip(tz));
    end
    if pad(6) || OPTION.force_taper
        IMAGE(:,:,end-l(6)+1:end-crop(6)) = IMAGE(:,:,end-l(6)+1:end-crop(6)) .* tz + val.*(1-tz);
    end
catch
    raiseError(IMAGE, LIMITS, OPTION);
end

end  % applyTaper_real3d


function IMAGE = applyTaper_real2d(IMAGE, LIMITS, OPTION, pad, crop, val)

l = numel(OPTION.taper) + crop;  % broadcast
ty = OPTION.taper;
tx = ty';

try
    % top, down, left, right
    if pad(1) || OPTION.force_taper
        IMAGE(1+crop(1):l(1), :)        = IMAGE(1+crop(1):l(1), :) .* flip(tx) + val.*(1-flip(tx));
    end
    if pad(2) || OPTION.force_taper
        IMAGE(end-l(2)+1:end-crop(2),:) = IMAGE(end-l(2)+1:end-crop(2),:) .* tx + val.*(1-tx);
    end
    if pad(3) || OPTION.force_taper
        IMAGE(:, 1+crop(3):l(3))        = IMAGE(:, 1+crop(3):l(3)) .* flip(ty) + val.*(1-flip(ty));
    end
    if pad(4) || OPTION.force_taper
        IMAGE(:,end-l(4)+1:end-crop(4)) = IMAGE(:,end-l(4)+1:end-crop(4)) .* ty + val.*(1-ty);
    end
catch
    raiseError(IMAGE, LIMITS, OPTION);
end

end  % applyTaper_real2d


function [OPTION, flg, ndim] = checkIN(IMAGE, LIMITS, OPTION)
% IMAGE
if ~isnumeric(IMAGE)
    error('EMC:IMAGE', 'IMAGE should be numeric, got %s', class(IMAGE))
elseif isvector(IMAGE) || isscalar(IMAGE)
    error('EMC:IMAGE', 'IMAGE should be 2d or 3d, got size: %s', mat2str(size(IMAGE)))
end

[flg.is3d, ~, ndim] = EMC_is3d(size(IMAGE));

% LIMITS
if ~isnumeric(LIMITS)
    error('EMC:LIMITS', 'For a %dd IMAGE, LIMITS should be a row vector of %d integers, got %s', ...
          ndim, ndim*2, class(LIMITS));
elseif ~isrow(LIMITS) || any(isnan(LIMITS)) || any(isinf(LIMITS)) || ...
       any(rem(LIMITS,1)) || numel(LIMITS) ~= ndim * 2
    error('EMC:LIMITS', 'For a %dd IMAGE, LIMITS should be a row vector of %d integers, got %s', ...
          ndim, ndim*2, mat2str(LIMITS));
end

if any(LIMITS > 0)
    flg.pad = true;
else
    flg.pad = false;
end

if any(LIMITS < 0)
    flg.crop = true;
else
    flg.crop = false;
end

% Optinal parameters.
OPTION = EMC_getOption(OPTION, {'origin', 'value', 'taper', 'force_taper', 'precision'}, false);

% origin
if isfield(OPTION, 'origin')
    if isscalar(OPTION.origin) || isnumeric(OPTION.origin)
        if OPTION.origin == -1
            flg.fft = true;
        elseif OPTION.origin == 1 || OPTION.origin == 2
            flg.fft = false;
        else
            error('EMC:origin', 'OPTION.origin should be 1, 2, or -1')
        end
    else
        error('EMC:origin', 'OPTION.origin should be 1, 2, or -1')
    end
else
    OPTION.origin = 1;  % default
    flg.fft = false;
end

% precision
[precision, flg.gpu, ~] = EMC_getClass(IMAGE);
if isfield(OPTION, 'precision')
    if ~(ischar(OPTION.precision) || isstring(OPTION.precision)) || ...
       ~(strcmpi('single', OPTION.precision) || strcmpi('double', OPTION.precision))
      	error('EMC:precision', "OPTION.precision should be 'single' or 'double'")
    end
else
    OPTION.precision = precision;
end

% value
if isfield(OPTION, 'value')
    if strcmpi(OPTION.value, 'uniform')
        flg.uniform = true;
    elseif strcmpi(OPTION.value, 'mean')
        OPTION.value = EMC_setPrecision(mean(IMAGE, 'all'), OPTION.precision);
        flg.uniform = false;
    elseif isnumeric(OPTION.value) && isscalar(OPTION.value)
        flg.uniform = false;
        OPTION.value = EMC_setPrecision(OPTION.value, OPTION.precision);
        if ~flg.gpu && EMC_isOnGpu(OPTION.value)
            OPTION.value = gather(OPTION.value);
        end
    else
        error('EMC:value', "OPTION.value should be a float|int, 'uniform' or 'mean'")
    end
else
    OPTION.value = 0;  % default
    flg.uniform = false;
end

% taper
if isfield(OPTION, 'taper')
    % bool
    if islogical(OPTION.taper) && isscalar(OPTION.taper)
        if OPTION.taper
            OPTION.taper = EMC_taper('cosine', 7, {});  % default
            flg.taper = true;
        else
            flg.taper = false;
        end
    % {type, size}
    elseif iscell(OPTION.taper)
        if numel(OPTION.taper) ~= 2
            error('EMC:taper', 'OPTION.taper should be {type(str), size(int)}.')
        else
            OPTION.taper = EMC_taper(OPTION.taper{1}, OPTION.taper{2}, {});
            flg.taper = true;
        end
    % vector|int|float: own taper
    elseif isnumeric(OPTION.taper) && isrow(OPTION.taper)
        if ~flg.gpu && EMC_isOnGpu(OPTION.taper)
            OPTION.taper = gather(OPTION.taper);
        end
        flg.taper = true;
    else
        error('EMC:taper', 'OPTION.taper should be bool, cell or a row vector, got %s', class(OPTION.taper))
    end
else
     OPTION.taper = EMC_taper('cosine', 7, {});  % default
     flg.taper = true;
end

% force_taper
if isfield(OPTION, 'force_taper')
    if ~islogical(OPTION.force_taper) || ~isscalar(OPTION.force_taper)
        error('EMC:force_taper', ...
              'OPTION.force_taper should be a boolean, got %s', class(OPTION.force_taper))
    end
else
    OPTION.force_taper = false;  % default
end

end  % checkIN


function raiseError(IMAGE, LIMITS, OPTION)
% Error handling when applying taper: the input IMAGE should have a minimum size.
% Generate an useful error message for debugging and tests.
%

inSize = size(IMAGE);
ndim = length(inSize);

% CASE 1: the taper is too large given the input IMAGE.
if OPTION.origin == -1 && any(numel(OPTION.taper) > floor(inSize/2))
    error('EMC:taper', 'For a size of %s, the maximum taper for a fft IMAGE is %s, got %d', ...
          mat2str(inSize), mat2str(floor(inSize/2)), numel(OPTION.taper));
elseif any(numel(OPTION.taper) > inSize)
    error('EMC:taper', 'For a size of %s, the maximum taper for an IMAGE is %s, got %d', ...
          mat2str(inSize), mat2str(inSize), numel(OPTION.taper));
end

extendedSize = reshape(ones(2, ndim) .* inSize, 1, []);

if OPTION.force_taper
    taperToApply = ones(2, ndim*2) .* numel(OPTION.taper);
else
    if OPTION.origin == -1
        taperToApply = reshape(repmat(any(reshape(LIMITS > 0, 2, [])), 2, 1), 1, []);
        taperToApply = taperToApply .* numel(OPTION.taper);
    else
        taperToApply = (LIMITS > 0) .* numel(OPTION.taper);
    end
end

% CASE 2: For an fft IMAGE, the size(taper + cropping) should be smaller or equal
%         than half of the size of the IMAGE.
if OPTION.origin == -1
    halfSize = reshape(floor(extendedSize / 2), 2, []);
    halfSize(1, :) = halfSize(1, :) + mod(inSize, 2);  % count extra pixel if odd;
    halfSize = reshape(halfSize, 1, []);
    maxCrop = -1 .* (halfSize - taperToApply);
    
% CASE 3: For a real space IMAGE, if the size of the IMAGE is smaller than size (cropping + tapter),
%         applying the taper will raise an index error.
else
    maxCrop = -1 .* (extendedSize - taperToApply);
end

error('EMC:taper', ['One axis is too small given the cropping and taper required.\n\nGiven the inputs ', ...
      '(IMAGE size: %s and taper size: %s),\nthe maximum cropping for each edges is %s, but got %s.'], ...
      mat2str(inSize), mat2str(taperToApply), mat2str(maxCrop), mat2str((LIMITS < 0) .* LIMITS));

end % raiseError
