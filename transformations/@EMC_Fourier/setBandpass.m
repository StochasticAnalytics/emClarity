function obj = setBandpass(obj, PIXEL, HIGHPASS, LOWPASS, OPTION)
%
% obj = obj.setBandpass(PIXEL, HIGHPASS, LOWPASS, OPTION)
% Create a 2D/3D filter; either a lowpass, highpass or bandpass filter.
%
% Input:
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
%   OPTION (cell|struct):               Optional parameters.
%                                       If cell: {field, value; ...}, note the ';' between parameters.
%                                       NOTE: Can be empty.
%                                       NOTE: Unknown fields will raise an error.
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
% Output:
%   Bandpass (numeric):                 Bandpass filter.
%                                       If obj.half, the bandpass is non-redundant.
%                                       If obj.centered, the banpass is centered.
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
%   - >> ft = EMC_Fourier([3740,3838], 'gpu', {});
%     >> ft = ft.setBandpass(3, nan, 8, {});  % lowpass filter with cutoff at 8A
%     >> img = ft.applyBandpass(img, {});
%
% Other EMC-files required:
%   EMC_getOption, EMC_coordVectors
%
% See also EMC_Fourier.applyBandpass.
%

% Created:  18Jan2020, R2019a
% Version:  v.1.0.  unittest (TF, 2Feb2020).
%           v.1.0.1 LOWPASS can be 0, which is now equivalent to NaN (no lowpass filter).
%                   If no highpass and no lowpass, compute a full pass (TF, 17Feb2020).
%           v.1.1.  EMC_getBandpass is integrated into EMC_Fourier.setBandpass (TF, 8Mar2020).
%

%% checkIN
defaultSigma = 0.02;
[OPTION, flg] = checkIN(PIXEL, HIGHPASS, LOWPASS, OPTION, defaultSigma);

% If no lowpass and no highpass, return ones.
if ~flg.lowpass && ~flg.highpass
    if obj.isOnGpu
        obj.bandpass = ones(obj.size_freq, obj.precision, 'gpuArray');
    else
        obj.bandpass = ones(obj.size_freq, obj.precision);
    end
    return
end

if obj.centered; origin = -1; else; origin = 1; end
[vX, vY, vZ] = EMC_coordVectors(obj.size_freq, obj.method, {'origin', origin; ...
                                                            'half', obj.half; ...
                                                            'normalize', true; ...
                                                            'precision', obj.precision});

% The radial grid is the same for both highpass and lowpass filters.
if obj.is3d
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

    obj.bandpass = (radius > highpassCut) + (radius <= highpassCut) .* gaussian(radius, highpassCut, sigma);
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
        obj.bandpass = (radius < lowpassCut) .* obj.bandpass + ...
                   (radius >= lowpassCut) .* gaussian(radius, lowpassCut, sigma);
    else
        obj.bandpass = (radius < lowpassCut) + ...
                   (radius >= lowpassCut) .* gaussian(radius, lowpassCut, sigma);
    end
end

end  % setBandpass


function [OPTION, flg] = checkIN(PIXEL, HIGHPASS, LOWPASS, OPTION, defaultSigma)

% PIXEL
if ~isscalar(PIXEL) || ~isnumeric(PIXEL) || isinf(PIXEL) || ~(PIXEL > 0)
    error('EMC:PIXEL', 'PIXEL should be a positive float|int')
end

OPTION = EMC_getOption(OPTION, {'highpassRoll', 'highpassThresh', ...
                                'lowpassRoll', 'lowpassThresh'}, false);

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
