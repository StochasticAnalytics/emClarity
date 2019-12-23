function [IMAGE] = EMC_FilterCenterNormalize(IMAGE, BANDPASS, OPTION)
%
% Apply a BANDPASS filter to a real (real space and real|not-complex) IMAGE, set mean=0 and variance=1.
%
% Before applying the filter to the IMAGE, EMC_resize can be called to pad/crop
% or apply a taper to the IMAGE (see OPTION).
%
% IMAGE (numeric):                  2d/3d image to filter.
%
% BANDPASS (numeric|nan):           Bandpass filter to apply to the IMAGE.
%                                   Should correspond to IMAGE.
%                                   If nan, do not apply any bandpass filter.
%                                   NOTE: The bandpass should be not-centered (origin=-1) as it is
%                                         directely applied to the FT of the IMAGE.
%                                   NOTE: Must be on the same device than IMAGE.
%
% OPTION (cell | struct):           Optional parameters.
%                                   If cell: {field,value ; ...}, note the ';' between parameters.
%                                   NOTE: Can be empty.
%                                   NOTE: Unknown fields will raise an error.
%
%   -> 'do_ifft' (bool):            By default (do_ifft=true), the filtered IMAGE is switched back to
%                                   real space (ifft) before returning. If false, the IMAGE is left in
%                                   Fourier space (zero frequency first).
%                                   default = true
%
%   -> 'center_normalize' (bool):   Set the mean to 0 and variance to 1 of the output real space IMAGE.
%                                   default = true
%
%   -> 'limits' (vector):           IMAGE limits.
%                                   See EMC_resize:LIMITS for more details.
%                                   Padding/cropping happens BEFORE Fourier transformation (on the real
%                                   space input IMAGE).
%                                   default = no padding/cropping
%
%   -> 'value' (float|str):         Value to pad with.
%                                   See EMC_resize:OPTION:value for more details.
%                                   NOTE: This is ignored if nothing to pad.
%                                   default = 'mean'
%
%   -> 'taper' (bool|cell|vector):  Apply a taper to the edges of the IMAGE that are padded.
%                                   See EMC_resize:OPTION:taper for more details.
%                                   NOTE: This is ignored if nothing to pad and 'force_taper'=false.
%                                   default = true
%
%   -> 'force_taper' (bool):        Apply the taper to all edges of IMAGE, whether or not they are
%                                   padded/cropped. See EMC_resize:OPTION:force_taper.
%                                   default = false
%
%   -> 'precision' (str):           Precision of the output IMAGE.
%                                   NOTE: This is changed before Fourier transformation.
%                                   default = same as input IMAGE.
%

%% checkIN

% Extract optional parameters
OPTION = EMC_extract_option(OPTION, {'do_ifft', 'center_normalize', 'limits', ...
                                     'value', 'taper', 'force_taper', 'precision'}, false);

if isfield(OPTION, 'do_ifft')
    if ~islogical(OPTION.do_ifft)
        error('do_ifft should be a boolean, got %s', class(OPTION.do_ifft))
    end
else
    OPTION.do_ifft = true;  % default
end

if isfield(OPTION, 'center_normalize')
    if ~islogical(OPTION.center_normalize)
        error('center_normalize should be a boolean, got %s', class(OPTION.center_normalize))
    end
else
    OPTION.center_normalize = true;  % default
end

% checkIN will be done by EMC_resize.
if ~isfield(OPTION, 'limits')
    OPTION.limits = zeros(1, ndims(IMAGE)*2);  % default
end

% 'value', 'taper', 'force_taper' and 'precision' are options for EMC_resize.
option_resize = EMC_extract_option(OPTION, {'value', 'taper', 'force_taper', 'precision'}, true);

%% Resize, apply bandpass and center/normalize if wished.
if ~isnan(BANDPASS)
    IMAGE = fftn(EMC_resize(IMAGE, OPTION.limits, option_resize)) .* BANDPASS;
else
    IMAGE = fftn(EMC_resize(IMAGE, OPTION.limits, option_resize));
end

if (OPTION.center_normalize)
    IMAGE(1) = 0;  % Mean = 0
    IMAGE = IMAGE ./ (sqrt(sum(sum(sum(abs(IMAGE).^2)))) ./ numel(IMAGE));  % Variance = 1
end

if (OPTION.do_ifft)
    IMAGE = real(ifftn(IMAGE));
end

end  % end EMC_applyFilter
