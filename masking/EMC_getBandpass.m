function BANDPASS = EMC_getBandpass(SIZE, PIXEL, HIGHPASS, LOWPASS, METHOD, OPTION)
%
% BANDPASS = EMC_getBandpass(SIZE, PIXEL, HIGHPASS, LOWPASS, METHOD, OPTION)
% Create a 2D/3D filter; either a lowpass, highpass or bandpass filter.
%
% WARNING: EMC functions usually set the default origin to 1. This function has the default set to -1
%          as bandpass filters are usually directly applied to fft outputs (not centered).
%
% Input:
%   SIZE (vector):                      Size (in pixel) of the filter to compute; [x, y, z] or [x, y].
%                                       NOTE: [1, N] or [N, 1] is not allowed.
%
%   PIXEL (float):                      Pixel size in A/pix.
%
%   HIGHPASS (float|nan):               Cutoff, in Angstrom, of the high pass filter.
%                                       This cutoff corresponds to resolution where the pass is
%                                       fully recovered. Must be higher (lower frequency) than LOWPASS.
%                                       If nan or 0, the filter will be a low pass filter.
%
%   LOWPASS (float|str|nan):            Cutoff, in Angstrom, of the low pass filter.
%                                       This cutoff corresponds to resolution where the pass starts
%                                      	to roll off. Must be lower (higher frequency) than HIGHPASS.
%                                       If nan or 0, the filter will be a high pass filter.
%                                       If 'nyquist', starts the roll off at the Nyquist frequency.
%                                       NOTE: should be higher (lower resolution) or equal to the
%                                             Nyquist frequency.
%
%   METHOD (str):                       Device to use; 'gpu' or 'cpu'.
%
%   OPTION (cell|struct):               Optional parameters.
%                                       If cell: {field, value; ...}, note the ';' between parameters.
%                                       NOTE: Can be empty.
%                                       NOTE: Unknown fields will raise an error.
%
%     -> 'origin' (int):                Origin convention. See EMC_coordVectors for more details.
%                                       default = -1
%
%     -> 'half' (bool):                 Compute only the half mask. See EMC_coordVectors for more details.
%                                       default = false
%
%     -> 'highpassRoll' (float|str):	Standard deviation of the HIGHPASS gaussian roll.
%                                       If 'extended': extend the gaussian roll from the zero frequency
%                                                      to HIGHPASS.
%                                       NOTE: This is ignored if there is no HIGHPASS filter.
%                                       NOTE: The radial grids used to compute the gaussian are normalized
%                                             between [0,0.5]. The sigma should correspond to this range.
%                                       default = 0.02
%
%     -> 'highpassThresh'(float):   	Sets the strength of the zero frequency, between 0 and 1.
%                                       NOTE: This is ignored if there is no HIGHPASS filter and if
%                                             'highpassRoll' is not equal to 'extended'.
%                                       NOTE: If 0, the Gaussian will go to Inf, therefore the roll off
%                                             becomes smaller than one pixel, which is equivalent to
%                                             highpassRoll = 0.
%                                       default = 1e-3
%
%     -> 'lowpassRoll' (float|str):     Standard deviation of the LOWPASS gaussian roll.
%                                       If 'extended': extend the gaussian roll from the LOWPASS cutoff
%                                                      to Nyquist.
%                                       If LOWPASS is at Nyquist, the default gaussian roll is used anyway.
%                                       NOTE: This is ignored if there is no LOWPASS filter.
%                                       NOTE: The radial grids used to compute the gaussian are normalized
%                                             between [0,0.5]. The sigma should correspond to this range.
%                                       default = 0.02
%
%     -> 'lowpassThresh'(float):        Sets the strength of the Nyquist frequency, between 0 and 1.
%                                       NOTE: This is ignored if there is no LOWPASS filter and if
%                                             'lowpassRoll' is not equal to 'extended'.
%                                       NOTE: If 0, the Gaussian will go to Inf, therefore the roll off
%                                             becomes smaller than one pixel, which is equivalent to
%                                             lowpassRoll = 0.
%                                       default = 1e-3
%
%     -> 'precision' (str):             Precision of the output MASK; 'single' or 'double'.
%                                       default = 'single'
%
% Output:
%   Bandpass (numeric):                 Bandpass filter of desired SIZE, METHOD and precision.
%
% Note:
%   - The default gaussian std (lowroll and highroll) of 0.02 and window (lowthresh and highthresh)
%     of 0.001 gives a roll off large enough to 'remove' the Gibbs Phenomenon (ringing effect).
%     These filters are far from perfect filter. One improvement to reduce unwanted frequency,
%     could be to shift the cutoffs at ~0.8 or ~0.7 and not at 1 as they are currently.
%
%   - The size of the gaussian roll is relative to the size of the filter. It might cause an issue
%     if the filter is very small (<100pixel wide) where the roll off is only over a few pixels.
%     The current solution is to increase the low|high roll. One other solution could be to make
%     sure the taper is at least ~10pixels.
%
% Example:
%   - BANDPASS = EMC_bandpass([3740,3838], 3, nan, 8, 'cpu', {});  % lowpass filter with cutoff at 8A
%   - BANDPASS = EMC_bandpass([3740,3838], 3, 10, 8, 'cpu', {'origin', 1});
%
% Other EMC-files required:
%   EMC_getOption, EMC_coordVectors
%
% See also EMC_applyBandpass.
%

% Created:  18Jan2020, R2019a
% Version:  v.1.0.  unittest (TF, 2Feb2020).
%           v.1.0.1 LOWPASS can be 0, which is now equivalent to NaN (no lowpass filter).
%                   If no highpass and no lowpass, compute a full pass (TF, 17Feb2020).
%

%% checkIN
defaultSigma = 0.02;
[SIZE, OPTION, flg] = checkIN(SIZE, PIXEL, HIGHPASS, LOWPASS, METHOD, OPTION, defaultSigma);

% If no lowpass and no highpass, return ones.
if ~flg.lowpass && ~flg.highpass
    if strcmpi(METHOD, 'gpu')
        BANDPASS = ones(SIZE, OPTION.precision, 'gpuArray');
    else
        BANDPASS = ones(SIZE, OPTION.precision);
    end
    return
end

[vX, vY, vZ] = EMC_coordVectors(SIZE, METHOD, {'origin', OPTION.origin; ...
                                               'half', OPTION.half; ...
                                               'normalize', true; ...
                                               'precision', OPTION.precision});

% The radial grid is the same for both highpass and lowpass filters.
if flg.is3d
    radius = sqrt(vX'.^2 + vY.^2 + reshape(vZ,1,1,[]).^2);
else
    radius = sqrt(vX'.^2 + vY.^2);
end

% Normal gaussian
gaussian = @(x,m,s) exp(-1.*(x-m).^2 ./ (2.*s.^2));

%% High pass filter
if flg.highpass
    highpassCut = PIXEL / HIGHPASS;  % [1/pix]

    % Gaussian roll
    if strcmpi(OPTION.highpassRoll, 'extended')  % roll, up to the zero frequency
        sigma = sqrt(-1 * highpassCut^2 / (2 * log(OPTION.highpassThresh)));
        % Make sure the extended roll off is larger than default.
        if sigma < defaultSigma
            sigma = defaultSigma;
        end
    else  % classic gaussian roll
        sigma = OPTION.highpassRoll;
    end

    BANDPASS = (radius > highpassCut) + (radius <= highpassCut) .* gaussian(radius, highpassCut, sigma);
end

%% Lowpass filter
if flg.lowpass
    if strcmpi(LOWPASS, 'nyquist')
        lowpassCut = 0.5;
    else
        lowpassCut = PIXEL / LOWPASS;  % [1/pix]
    end
    
    % Gaussian roll
    if strcmpi(OPTION.lowpassRoll, 'extended')  % roll, up to nyquist
        sigma = sqrt(-1 * ((0.5 - lowpassCut).^2) / (2 * log(OPTION.lowpassThresh)));
        % Make sure the extended roll off is larger than default.
        if sigma < defaultSigma
            sigma = defaultSigma;
        end
    else  % classic gaussian roll
        sigma = OPTION.lowpassRoll;
    end

    if flg.highpass
        BANDPASS = (radius < lowpassCut) .* BANDPASS + ...
                   (radius >= lowpassCut) .* gaussian(radius, lowpassCut, sigma);
    else
        BANDPASS = (radius < lowpassCut) + ...
                   (radius >= lowpassCut) .* gaussian(radius, lowpassCut, sigma);
    end
end

end  % EMC_bandpass


function [SIZE, OPTION, flg] = checkIN(SIZE, PIXEL, HIGHPASS, LOWPASS, METHOD, OPTION, defaultSigma)

[flg.is3d, SIZE] = EMC_is3d(SIZE);
if SIZE(1) == 1 || SIZE(2) == 1
    error('EMC:SIZE', 'SIZE should be the size of a 2d or 3d array, got size %s', mat2str(SIZE))
end

if strcmpi(METHOD, 'gpu')
   	flg.gpu = true;
elseif strcmpi(METHOD, 'cpu')
    flg.gpu = false;
else
    error('EMC:METHOD', "METHOD should be 'cpu' or 'gpu'")
end

% PIXEL
if ~isscalar(PIXEL) || ~isnumeric(PIXEL) || isinf(PIXEL) || ~(PIXEL > 0)
    error('EMC:PIXEL', 'PIXEL should be a positive float|int')
end

OPTION = EMC_getOption(OPTION, {'origin', 'half', 'highpassRoll', 'highpassThresh', ...
                                'lowpassRoll', 'lowpassThresh', 'precision'}, false);

% let EMC_coordVectors do the checkIN, just set the default.
if ~isfield(OPTION, 'origin')
    OPTION.origin = -1;  % default
end
if ~isfield(OPTION, 'half')
    OPTION.half = false;
end
if ~isfield(OPTION, 'precision')
    OPTION.precision = 'single';
end

% HIGHPASS cutoff
if ~isscalar(HIGHPASS)
    error('EMC:HIGHPASS', 'HIGHPASS should be a scalar, got a %s of size %s', ...
          class(HIGHPASS), mat2str(size(HIGHPASS)))
elseif isnan(HIGHPASS) || HIGHPASS == 0
    flg.highpass = false;
elseif isnumeric(HIGHPASS) && ~isinf(HIGHPASS) && HIGHPASS > 0
    if HIGHPASS < PIXEL * 2
        error('EMC:HIGHPASS', ...
              'HIGHPASS should be greater (lower resolution) or equal to Nyquist (%.3f)', PIXEL * 2)
    end
    flg.highpass = true;
else
    error('EMC:HIGHPASS', 'HIGHPASS should be a positive float|int or nan')
end

% LOWPASS cutoff
if isscalar(LOWPASS) && isnumeric(LOWPASS) && ~isinf(LOWPASS)
    if isnan(LOWPASS) || LOWPASS == 0
        flg.highpass = false;
    elseif LOWPASS < PIXEL * 2
        error('EMC:LOWPASS', 'LOWPASS should be greater (lower resolution) or equal to Nyquist (%.3f)', ...
              PIXEL * 2)
    elseif flg.highpass && HIGHPASS < LOWPASS
        error('EMC:LOWPASS', 'LOWPASS should be smaller (high frequency) than HIGHPASS')
    end
    flg.lowpass = true;
elseif strcmpi(LOWPASS, 'nyquist')
    flg.lowpass = true;
else
    error('EMC:LOWPASS', "LOWPASS should be a non negative float|int, nan or 'nyquist'")
end

% highpassRoll
if isfield(OPTION, 'highpassRoll')
    if ~(isscalar(OPTION.highpassRoll) && isnumeric(OPTION.highpassRoll) ...
         && OPTION.highpassRoll >= 0 && ~isinf(OPTION.highpassRoll)) ...
       && ...
       ~((isstring(OPTION.highpassRoll) || ischar(OPTION.highpassRoll)) && ...
         strcmpi(OPTION.highpassRoll, 'extended'))
        error('EMC:highpassRoll', "OPTION.highpassRoll should be a positive float|int or 'extended'")
    end
else
    OPTION.highpassRoll = defaultSigma;  % default
end

% highpassThresh
if isfield(OPTION, 'highpassThresh')
    if ~isscalar(OPTION.highpassThresh) && ~isnumeric(OPTION.highpassThresh)
        error('EMC:highpassThresh', 'OPTION.highpassThresh should be a float|int, got %s of size %s', ...
              class(OPTION.highpassThresh), mat2str(size(OPTION.highpassThresh)))
    elseif OPTION.highpassThresh < 0 || OPTION.highpassThresh >= 1 || isnan(OPTION.highpassThresh)
        error('EMC:highpassThresh', 'OPTION.highpassThresh should be between 0 and 1, got %.3f', ...
              OPTION.highpassThresh)
    end
else
    OPTION.highpassThresh = 1e-3;  % default
end

% lowpassRoll
if isfield(OPTION, 'lowpassRoll')
    if ~(isscalar(OPTION.lowpassRoll) && isnumeric(OPTION.lowpassRoll) ...
         && OPTION.lowpassRoll >= 0 && ~isinf(OPTION.lowpassRoll)) ...
       && ...
       ~((isstring(OPTION.lowpassRoll) || ischar(OPTION.lowpassRoll)) && ...
         strcmpi(OPTION.lowpassRoll, 'extended'))
        error('EMC:lowpassRoll', "OPTION.lowpassRoll should be a positive float|int or 'extended'")
    end
else
    OPTION.lowpassRoll = defaultSigma;  % default
end

% lowpassThresh
if isfield(OPTION, 'lowpassThresh')
    if ~isscalar(OPTION.lowpassThresh) && ~isnumeric(OPTION.lowpassThresh)
        error('EMC:lowpassThresh', 'OPTION.lowpassThresh should be a float|int, got %s of size %s', ...
              class(OPTION.lowpassThresh), mat2str(size(OPTION.lowpassThresh)))
    elseif OPTION.lowpassThresh < 0 || OPTION.lowpassThresh >= 1 || isnan(OPTION.lowpassThresh)
        error('EMC:lowpassThresh', 'OPTION.lowpassThresh should be between 0 and 1, got %.3f', ...
              OPTION.lowpassThresh)
    end
else
    OPTION.lowpassThresh = 1e-3;  % default
end

end  % checkIN
