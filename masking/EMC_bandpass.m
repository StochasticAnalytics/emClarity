function [BANDPASS] = EMC_bandpass(SIZE, PIXEL, LOWCUT, HIGHCUT, METHOD, OPTION)
% [BANDPASS] = EMC_bandpass(SIZE, PIXEL_SIZE, LOWCUT, HIGHCUT, METHOD, OPTION)
%
% Create a 2D/3D filter; either lowpass, highpass or bandpass filter.
%
% WARNING: EMC functions usually set the default origin to 1. This function has the default set to -1
%          as bandpass filters are usually directly applied to fft outputs (not centered).
%
% SIZE (vector):                Size of the filter; (x, y, z) or (x, y).
%                               z=1 is still considered has 2d.
%
% PIXEL (float):                Sampling frequency (A/pix).
%
% LOWCUT (float|nan):           Resolution (in Angstrom) where the pass starts to recover.
%                               Must be higher (lower frequency) than HIGHCUT.
%                               If nan, act as a low pass filter.
%
% HIGHCUT (float|str|nan):      Resolution (in Angstrom) where the pass starts to roll off.
%                               Must be lower (higher frequency) than LOWCUT.
%                               If nan, act as a high pass filter.
%                               If 'nyquist', starts the roll off at the Nyquist frequency.
%                               NOTE: should be higher (lower resolution) or equal to the Nyquist frequency.
%
% METHOD (str):                 Device to use; 'gpu' or 'cpu'.
%
% OPTION (cell|struct):         Optional parameters.
%                               If cell: {field,value ; ...}, note the ';' between parameters.
%                               NOTE: Can be empty.
%                               NOTE: Unknown fields will raise an error.
%
% -> 'origin' (int):            Origin convention - Center of rotation.
%                               -1, 0, 1 or 2; see EMC_multi_gridVectors for more details.
%                               default = -1
%
% -> 'half' (bool):             Compute only half mask.
%                               default = false
%
% -> 'lowroll' (float|str):     Standard deviation of the gaussian roll at LOWCUT.
%                               See EMC_shapeMask for more details.
%                               If 'extended', extend the gaussian roll from the zero frequency to LOWCUT.
%                               NOTE: This is ignored if LOWCUT is nan or 0.
%                               NOTE: 0 is equivalent to no roll off.
%                               NOTE: The radial grids are normalized between [-0.5,0.5].
%                               default = 0.02
%
% -> 'lowthresh'(float):        If 'lowroll' = 'extended', sets the strength of the zero frequency,
%                               between 0 and 1.
%                               NOTE: This is ignored if LOWCUT is nan or 0 and if 'lowroll' ~= 'extended'.
%                               NOTE: If 0, the Gaussian will go to Inf, therefore the roll off becomes
%                                     smaller than one pixel (no roll off).
%                                     This is equivalent to lowroll = 0.
%                               default = 1e-3
%
% -> 'highroll' (float|str):    Standard deviation of the gaussian roll at HIGHCUT.
%                               See EMC_shapeMask for more details.
%                               If 'extended', extend the gaussian roll from HIGHCUT to Nyquist. 
%                               If HIGHCUT is at Nyquist, the default gaussian roll is used anyway.
%                               NOTE: This is ignored if HIGHCUT is nan.
%                               NOTE: 0 is equivalent to no roll off.
%                               NOTE: The radial grids are normalized between [-0.5,0.5].
%                               default = 0.02
%
% -> 'highthresh'(float):       If 'highroll' = 'extended', sets the strength of the Nyquist frequency,
%                               between 0 and 1.
%                               NOTE: This is ignored if HIGHCUT is nan or 0 and if 'highroll' ~= 'extended'.
%                               NOTE: If 0, the Gaussian will go to Inf, therefore the roll off becomes
%                                     smaller than one pixel (no roll off).
%                                     This is equivalent to highroll = 0.
%                               default = 1e-3
%
% -> 'precision' (str):         Precision of the output MASK; 'single' or 'double'.
%                               default = 'single'
%
%-------
% TODO 1:                       The gaussian taper is invariant to the size of the filter. It might cause
%                               an issue if the filter is very small (<100pixel wide). The current solution
%                               is to increase the low|high roll. If it is necessary, make sure the taper
%                               is at least ~7pixels.
%
%-------
% EXAMPLE:                      [BANDPASS] = EMC_bandpass([3740,3838], 3, 10, 8, 'cpu', {});
%                               [BANDPASS] = EMC_bandpass([3740,3838], 3, 10, 8, 'cpu', ...
%                                                         {'origin', 1; 'lowroll','extended'});
%
% See also EMC_shapeMask.

%% checkIN
default_lowRoll = 0.02;
default_highRoll = 0.02;

[OPTION, flg] = checkIN(SIZE, PIXEL, LOWCUT, HIGHCUT, METHOD, OPTION, ...
                        default_lowRoll, default_highRoll);

[vX, vY, vZ] = EMC_multi_vectorCoordinates(SIZE, METHOD, {'origin', OPTION.origin; ...
                                                          'half', OPTION.half; ...
                                                          'normalize', true; ...
                                                          'precision', OPTION.precision});
% The radial grid is the same for highpass and lowpass.
if (flg.is3d)
    radius = sqrt(vX'.^2 + vY.^2 + reshape(vZ,1,1,[]).^2);
else
    radius = sqrt(vX'.^2 + vY.^2);
end

% Normal gaussian
gaussian = @(x,m,s) exp(-1.*(x-m).^2 ./ (2.*s.^2));

if (flg.highpass)
    lowCut = PIXEL / LOWCUT;  % [1/pix]
    % Prevents the lowCut to be too close from the 0 which allows a >= 7pixel roll off.
    lowCut = max(7/min(SIZE), lowCut);

    % Low roll
    if strcmpi(OPTION.lowroll, 'extended')  % roll up to the zero frequency
        % set zero frequency to the lowthresh
        sigma = sqrt(-1 * lowCut^2 / (2 * log(OPTION.lowthresh)));
    else  % classic gaussian roll
        sigma = OPTION.lowroll;
    end

    BANDPASS = (radius > lowCut) + (radius <= lowCut) .* gaussian(radius, lowCut, sigma);
else
    if (flg.gpu)
        BANDPASS = ones(SIZE, OPTION.precision, 'gpuArray');
    else
        BANDPASS = ones(SIZE, OPTION.precision);
    end
end

if (flg.lowpass)
    if strcmpi(HIGHCUT, 'nyquist')
        highCut = 0.5;
    else
        highCut = PIXEL / HIGHCUT;  % [1/pix]
    end
    
    % High roll
    if strcmpi(OPTION.highroll, 'extended')  % roll up to nyquist
        % set nyquist to the highthresh
        sigma = sqrt(-1 * ((0.5 - highCut).^2) / (2 * log(OPTION.highthresh)));
        sigma_classic = default_highRoll;
        % Make sure the extended roll off is larger than default.
        if sigma < sigma_classic
            sigma = sigma_classic;
        end
    else  % classic gaussian roll
        sigma = OPTION.highroll;
    end

    BANDPASS = (radius < highCut) .* BANDPASS + (radius >= highCut) .* gaussian(radius, highCut, sigma);
end

clear radius
end  % EMC_bandpass


function [OPTION, flg] = checkIN(SIZE, PIXEL, LOWCUT, HIGHCUT, METHOD, OPTION, def_lowRoll, def_highRoll)

[SIZE, flg.is3d, ndim] = EMC_is3d(SIZE);
if strcmpi(METHOD, 'gpu')
   	flg.gpu = true;
elseif strcmpi(METHOD, 'cpu')
    flg.gpu = false;
else
    error("METHOD should be 'cpu' or 'gpu'")
end

validateattributes(SIZE, {'numeric'}, {'vector', 'numel', ndim, 'integer', 'positive'}, '', 'SIZE');
validateattributes(PIXEL, {'numeric'}, {'numel', 1, 'positive'}, '', 'PIXEL');

OPTION = EMC_extract_option(...
            OPTION, ...
            {'origin', 'half', 'lowroll', 'lowthresh', 'highroll', 'highthresh', 'precision'}, ...
            false);

% let EMC_multi_vectorCoordinates do the checkIN, just set the default.
if ~isfield(OPTION, 'origin')
    OPTION.origin = -1;  % default
end
if ~isfield(OPTION, 'half')
    OPTION.half = false;
end
if ~isfield(OPTION, 'precision')
    OPTION.precision = 'single';
end

% LOWCUT
if isnan(LOWCUT)
    flg.highpass = false;
elseif (isfloat(LOWCUT) || isinteger(LOWCUT)) && LOWCUT >= 0
    if LOWCUT == 0
        flg.highpass = false;
    end
    flg.highpass = true;
else
    error('LOWCUT should be a positive float|int or nan')
end

% HIGHCUT
if isnan(HIGHCUT)
    flg.lowpass = false;
elseif (isfloat(HIGHCUT) || isinteger(HIGHCUT)) && HIGHCUT >= 0
    if HIGHCUT == 0
        flg.lowpass = false;
    end
    if HIGHCUT < PIXEL * 2
        error('HIGHCUT should be higher (lower resolution) or equal to Nyquist (%.3f)', PIXEL * 2)
    elseif HIGHCUT > LOWCUT
        error('HIGHCUT should be lower (high frequency) than LOWCUT')
    end
    flg.lowpass = true;
elseif strcmpi(HIGHCUT, 'nyquist')
    flg.lowpass = true;
else
    error("HIGHCUT should be a positive float|int, nan or 'nyquist'")
end

% Low frequencies
if isfield(OPTION, 'lowroll')
    if isfloat(OPTION.lowroll) || isinteger(OPTION.lowroll)
        if OPTION.lowroll < 0
            error('lowroll should be positve, got %f', OPTION.lowroll)
        end
    elseif ~strcmpi(OPTION.lowroll, 'extended')
        error("lowroll should be a postive float|int or 'extended', got %s", OPTION.lowroll)
    end
else
    OPTION.lowroll = def_lowRoll;  % default
end

if isfield(OPTION, 'lowthresh')
    if ~(isfloat(OPTION.lowthresh) || isinteger(OPTION.lowthresh))
        error('lowthresh should be a float|int between 0 and 1, got %s', class(OPTION.lowthresh))
    elseif ~(0 <= OPTION.lowthresh <= 1)
        error('lowthresh should be between 0 and 1, got %f', OPTION.lowthresh)
    end
else
    OPTION.lowthresh = 1e-3;  % default
end

% High frequencies
if isfield(OPTION, 'highroll')
    if isfloat(OPTION.highroll) || isinteger(OPTION.highroll)
        if OPTION.highroll < 0
            error('highroll should be positve, got %f', OPTION.highroll)
        end
    elseif ~strcmpi(OPTION.highroll, 'extended')
        error("highroll should be a postive float|int or 'extended', got %s", OPTION.highroll)
    end
else
    OPTION.highroll = def_highRoll;  % default
end

if isfield(OPTION, 'highthresh')
    if ~(isfloat(OPTION.highthresh) || isinteger(OPTION.highthresh))
        error('highthresh should be a float|int between 0 and 1, got %s', class(OPTION.highthresh))
    elseif ~(0 <= OPTION.highthresh <= 1)
        error('highthresh should be between 0 and 1, got %f', OPTION.highthresh)
    end
else
    OPTION.highthresh = 1e-3;  % default
end

end  % checkIN
