function [OUT] = EMC_resize(IMAGE, LIMITS, OPTIONAL)
% [OUT] = EMC_resize(IMAGE, LIMITS, OPTIONAL)
%
% Pad and/or crop an IMAGE.
%
% IMAGE (single | double):          2d/3d image to pad an/or crop.
%
% LIMITS (vector):                  Number of pixels to pad or crop, for each axis.
%                                   If 2d IMAGE: [xleft, xright, yleft, yright].
%                                   If 3d IMAGE: [xleft, xright, yleft, yright, zleft, zright].
%                                   See EMC_multi_limits for more details.
%                                   NOTE: Positive values indicate the number of pixels to pad, while
%                                         negative values indicate the number of pixels to crop.
%
% OPTIONAL (cell | struct):         Optional parameters.
%                                   If cell: {field,value ; ...}, note the ';' between parameters.
%                                   NOTE: Can be empty.
%                                   NOTE: Unknown fields will raise an error.
%
%   -> 'origin' (int):              Origin convention - Center of rotation.
%                                   -1, 1 or 2; see EMC_multi_gridVectors for more details.
%                                   NOTE: 1|2 produce identical results because LIMITS is not relative
%                                         to the origin but to the edges of the IMAGE. On the other hand,
%                                         'origin'=-1 is different as it specifies a not-centered IMAGE
%                                         (zero frequency first). As such, the 'edges' of IMAGE are at
%                                         the center of the array.
%                                   NOTE: This function opperates in pixel space. As such, origin=0 is not
%                                         allowed.
%                                   defaut = 1
%
%   -> 'value' (float | str):       If float: value to pad with.
%                                   If 'uniform': pad with white gaussian noise.
%                                   NOTE: this parameter is ignored if no padding is applied.
%                                   default = 0
%
%   -> 'taper' (bool|cell|vector):  Apply a taper to the IMAGE padded edges before padding.
%                                   If bool: Apply or not the default taper; equivalent to {'cosine', 7}
%                                   If cell: {type, size} with type = 'linear' or 'cosine' and with 
%                                   size = size of the taper (in pixel).
%                                   If vector: ROW vector used as taper (left to right <-> center to edge).
%                                   NOTE: Only the edges that are padded are tapered; this parameter is
%                                         ignored if no padding is applied.
%                                   default = {'cosine', 7}
%
%   -> 'precision' (str):           Precision of the padded/cropped image; 'single' or 'double'.
%                                   NOTE: this is ignored if ~any(LIMITS) <=> no cropping nor padding.
%                                   default = same as IMAGE
%
%---------
% EXAMPLE:                          [OUT] = EMC_resize(randn(64,64), [10,10,-5,0], {})
%                                   [OUT] = EMC_resize(randn(64,64), [10,10,-5,0], ...
%                                                      {'value', 2; 'origin', -1})
%
% See also EMC_multi_limits.m

%% MAIN
[OPTIONAL, flg] = checkIN(IMAGE, LIMITS, OPTIONAL);

% Short cut: nothing to pad|crop
if ~flg.pad && ~flg.crop
    OUT = IMAGE;
    return
end

% Cropping and padding are done simultaneously.
IN_size = size(IMAGE);
OUT_size = IN_size + sum(reshape(LIMITS, 2, []));

% Option to pad with white gaussian noise.
if (flg.uniform)
    val = mean(IMAGE(:));
    std_ = std(IMAGE(:));
else
    val = OPTIONAL.value;
    std_ = nan;
end

% Allocate memory for the output image.
if (flg.gpu)
    if (flg.uniform)
       	OUT = randn(OUT_size, OPTIONAL.precision, 'gpuArray') .* std_ + val;
    else
      	OUT = zeros(OUT_size, OPTIONAL.precision, 'gpuArray');
        if OPTIONAL.value ~= 0
            OUT = OUT + OPTIONAL.value;
        end
   	end
else
  	if (flg.uniform)
       	OUT = randn(OUT_size, OPTIONAL.precision) .* std_ + val;
    else
      	OUT = zeros(OUT_size, OPTIONAL.precision);
        if OPTIONAL.value ~= 0
            OUT = OUT + OPTIONAL.value;
        end
  	end
end

x=1; y=2; z=3;
if (flg.is3d)
    if (flg.fft)  % not-centered IMAGE: 'origin' = -1
        crop = reshape(LIMITS .* (LIMITS < 0), 2, []);
        l = floor((IN_size+1) / 2) + crop(1);  % left side
        r = floor((IN_size-2) / 2) + crop(2);  % right side

        % Taper
        if (flg.taper)
            s = numel(OPTIONAL.taper);
            pad = LIMITS > 0;

            if any(pad(1:2) > 0)
                t = OPTIONAL.taper';
                tf = flip(t);
                IMAGE(l(x)-s+1:l(x), :, :)         = IMAGE(l(x)-s+1:l(x), :, :)         .*t  +val.*(1-t);
                IMAGE(end-r(x):end-r(x)+s-1, :, :) = IMAGE(end-r(x):end-r(x)+s-1, :, :) .*tf +val.*(1-tf);
            end
            if any(pad(3:4) > 0)
                t = OPTIONAL.taper;
                tf = flip(t);
                IMAGE(:, l(y)-s+1:l(y), :)         = IMAGE(:, l(y)-s+1:l(y), :)         .*t  +val.*(1-t);
                IMAGE(:, end-r(y):end-r(y)+s-1, :) = IMAGE(:, end-r(y):end-r(y)+s-1, :) .*tf +val.*(1-tf);
            end
            if any(pad(5:6) > 0)
                t = reshape(OPTIONAL.taper, 1, 1, []);
                tf = flip(t);
                IMAGE(:, :, l(z)-s+1:l(z))         = IMAGE(:, :, l(z)-s+1:l(z))         .*t  +val.*(1-t);
                IMAGE(:, :, end-r(z):end-r(z)+s-1) = IMAGE(:, :, end-r(z):end-r(z)+s-1) .*tf +val.*(1-tf);
            end
        end

        % Pad/Crop
        OUT(1:l(x),       1:l(y),       1:l(z))       = IMAGE(1:l(x),       1:l(y),       1:l(z));
        OUT(end-r(x):end, 1:l(y),       1:l(z))       = IMAGE(end-r(x):end, 1:l(y),       1:l(z));
        OUT(1:l(x),       end-r(y):end, 1:l(z))       = IMAGE(1:l(x),       end-r(y):end, 1:l(z));
        OUT(end-r(x):end, end-r(y):end, 1:l(z))       = IMAGE(end-r(x):end, end-r(y):end, 1:l(z));
        OUT(1:l(x),       1:l(y),       end-r(z):end) = IMAGE(1:l(x),       1:l(y),       end-r(z):end);
        OUT(end-r(x):end, 1:l(y),       end-r(z):end) = IMAGE(end-r(x):end, 1:l(y),       end-r(z):end);
        OUT(1:l(x),       end-r(y):end, end-r(z):end) = IMAGE(1:l(x),       end-r(y):end, end-r(z):end);
        OUT(end-r(x):end, end-r(y):end, end-r(z):end) = IMAGE(end-r(x):end, end-r(y):end, end-r(z):end);

    else  % centered IMAGE: 'origin' = 0|1|2
        crop = abs(LIMITS .* (LIMITS < 0));
        pad  = LIMITS .* (LIMITS > 0);

        if (flg.taper) && any(pad)
            l = numel(OPTIONAL.taper);
            ty = OPTIONAL.taper;
            tx = ty';
            tz = reshape(ty, 1, 1, []);

            % top, down, left, right, bottom, front
            if pad(1); IMAGE(1:l,:,:)         = IMAGE(1:l,:,:)         .* flip(tx) + val.*(1-flip(tx)); end
            if pad(2); IMAGE(end-l+1:end,:,:) = IMAGE(end-l+1:end,:,:) .* tx       + val.*(1-tx); end
            if pad(3); IMAGE(:,1:l,:)         = IMAGE(:,1:l,:)         .* flip(ty) + val.*(1-flip(ty)); end
            if pad(4); IMAGE(:,end-l+1:end,:) = IMAGE(:,end-l+1:end,:) .* ty       + val.*(1-ty); end
            if pad(5); IMAGE(:,:,1:l)         = IMAGE(:,:,1:l)         .* flip(tz) + val.*(1-flip(tz)); end
            if pad(6); IMAGE(:,:,end-l+1:end) = IMAGE(:,:,end-l+1:end) .* tz       + val.*(1-tz); end
        end

        % Pad/Crop
        if flg.pad && flg.crop
            OUT(1+pad(1):end-pad(2), ...
                1+pad(3):end-pad(4), ...
                1+pad(5):end-pad(6)) = IMAGE(1+crop(1):end-crop(2), ...
                                                   1+crop(3):end-crop(4), ...
                                                   1+crop(5):end-crop(6));
        elseif flg.pad  % pad only
            OUT(1+pad(1):end-pad(2), 1+pad(3):end-pad(4), 1+pad(5):end-pad(6)) = IMAGE;
        elseif flg.crop  % crop only
            OUT(:,:,:) = IMAGE(1+crop(1):end-crop(2), 1+crop(3):end-crop(4), 1+crop(5):end-crop(6));
        end
    end
else  % 2d
    if (flg.fft)  % not-centered IMAGE: 'origin' = -1
        crop = reshape(LIMITS .* (LIMITS < 0), 2, []);
        l = floor((IN_size+1) / 2) + crop(1);  % left side
        r = floor((IN_size-2) / 2) + crop(2);  % right side

        % Taper.
        if (flg.taper)
            s = numel(OPTIONAL.taper);
            pad = LIMITS > 0;
            if any(pad(1:2) > 0)
                t = OPTIONAL.taper';
                tf = flip(t);
                IMAGE(l(x)-s+1:l(x), :)         = IMAGE(l(x)-s+1:l(x), :)         .*t  +val.*(1-t);
                IMAGE(end-r(x):end-r(x)+s-1, :) = IMAGE(end-r(x):end-r(x)+s-1, :) .*tf +val.*(1-tf);
            end
            if any(pad(3:4) > 0)
                t = OPTIONAL.taper;
                tf = flip(t);
                IMAGE(:, l(y)-s+1:l(y))         = IMAGE(:, l(y)-s+1:l(y))         .*t  +val.*(1-t);
                IMAGE(:, end-r(y):end-r(y)+s-1) = IMAGE(:, end-r(y):end-r(y)+s-1) .*tf +val.*(1-tf);
            end
        end
        
        % Pad/Crop
        OUT(1:l(x),       1:l(y))       = IMAGE(1:l(x),       1:l(y));
        OUT(end-r(x):end, 1:l(y))       = IMAGE(end-r(x):end, 1:l(y));
        OUT(1:l(x),       end-r(y):end) = IMAGE(1:l(x),       end-r(y):end);
        OUT(end-r(x):end, end-r(y):end) = IMAGE(end-r(x):end, end-r(y):end);

    else  % centered IMAGE: 'origin' = 0|1|2
        crop = abs(LIMITS .* (LIMITS < 0));
        pad  = LIMITS .* (LIMITS > 0);

        if (flg.taper)
            l = numel(OPTIONAL.taper);
            ty = OPTIONAL.taper;
            tx = ty';
            % top, down, left, right
            if pad(1); IMAGE(1:l, :)         = IMAGE(1:l, :)         .* flip(tx) + val.*(1-flip(tx)); end
            if pad(2); IMAGE(end-l+1:end, :) = IMAGE(end-l+1:end, :) .* tx       + val.*(1-tx); end
            if pad(3); IMAGE(:, 1:l)         = IMAGE(:, 1:l)         .* flip(ty) + val.*(1-flip(ty)); end
            if pad(4); IMAGE(:, end-l+1:end) = IMAGE(:, end-l+1:end) .* ty       + val.*(1-ty); end
        end

        if flg.pad && flg.crop  % pad and crop
            OUT(1+pad(1):end-pad(2), 1+pad(3):end-pad(4)) = IMAGE(1+crop(1):end-crop(2), ...
                                                                  1+crop(3):end-crop(4));
        elseif flg.pad  % pad only
            OUT(1+pad(1):end-pad(2), 1+pad(3):end-pad(4)) = IMAGE;
        elseif flg.crop  % crop only
            OUT(:,:) = IMAGE(1+crop(1):end-crop(2), 1+crop(3):end-crop(4));
        end
    end
end

end  % EMC_resize


function [OPTIONAL, flg, ndim] = checkIN(IMAGE, LIMITS, OPTIONAL)
% Standard sanity check.

% LIMITS
edges = numel(LIMITS);
if edges == 6
    flg.is3d = true;
    ndim = 3;
elseif edges ==4
    flg.is3d = false;
    ndim = 2;
else
    error('LIMITS should be of size 4 (2d) or 6 (3d), got %d', edges)
end
validateattributes(LIMITS, {'numeric'}, {'integer', 'vector'}, 'checkIN', 'LIMITS')
if ndim ~= ndims(IMAGE)
    error('IMAGE (%fD) and LIMITS (%fD) do not correspond', ndims(IMAGE), ndim)
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

% Extract optional parameters
OPTIONAL = EMC_extract_optional(OPTIONAL, {'origin', 'value', 'taper', 'precision'});

if isfield(OPTIONAL, 'origin')
    if OPTIONAL.origin == -1
        flg.fft = true;
    elseif OPTIONAL.origin == -1 || OPTIONAL.origin == 1 || OPTIONAL.origin == 2
        flg.fft = false;
    else
        error("origin should be 1, 2, or -1, got %d", OPTIONAL.origin)
    end
else
    OPTIONAL.origin = 1;  % default
    flg.fft = false;
end

if isfield(OPTIONAL, 'value')
    if strcmpi(OPTIONAL.value, 'uniform')
        flg.uniform = true;
    elseif isinteger(OPTIONAL.value) || isfloat(OPTIONAL.value)
        flg.uniform = false;
    else
        error("value should be a float|int or 'uniform'")
    end
else
    OPTIONAL.value = 0;  % default
    flg.uniform = false;
end

if isfield(OPTIONAL, 'taper')
    % bool
    if islogical(OPTIONAL.taper)
        if OPTIONAL.taper
            OPTIONAL.taper = EMC_taper('cosine', 1, 0, 7);  % default
            flg.taper = true;
        else
            flg.taper = false;
        end
    % [type, size]
    elseif iscell(OPTIONAL.taper)
        if numel(OPTIONAL.taper) ~= 2
            error('taper not recognized.')
        else
            OPTIONAL.taper = EMC_taper(OPTIONAL.taper{1}, 1, 0, OPTIONAL.taper{2});
            flg.taper = true;
        end
    % vector: own taper
    else
        validateattributes(OPTIONAL.taper, {'numeric'}, {'row'})
        flg.taper = true;
    end
else
     OPTIONAL.taper = EMC_taper('cosine', 1, 0, 7);  % default
     flg.taper = true;
end

% extract precision and method.
if isfield(OPTIONAL, 'precision')
    if ~contain(['single', 'double'], OPTIONAL.precision)
        error("presision should be 'single' or 'double'")
    end
    if isa(IMAGE, 'gpuArray')
        flg.gpu = true;
    else
        flg.gpu = false;
    end
else
    OPTIONAL.precision = class(IMAGE);
    if strcmp(OPTIONAL.precision, 'gpuArray')
        flg.gpu = true;
        OPTIONAL.precision = classUnderlying(IMAGE);
    else
        flg.gpu = false;
    end
end

end  % checkIN