function [BANDPASS] = EMC_bandpass(SIZE, PIXEL_SIZE, LOWCUT, HIGHCUT, METHOD, OPTIONAL)
%
% Create a 2D/3D filter; either lowpass, highpass or bandpass filter.
%
%
% SIZE (vector):                Size of the filter; (x, y, z) or (x, y).
%                               z=1 is still considered has 2d.
%
% PIXEL_SIZE (float):           Sampling frequency (A/pix).
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
% OPTIONAL (cell|struct):       Optional parameters.
%                               If cell: {field,value ; ...}, note the ';' between parameters.
%                               NOTE: Can be empty.
%                               NOTE: Unknown fields will raise an error.
%
% -> 'origin' (int):            Origin convention - Center of rotation.
%                               -1, 0, 1 or 2; see EMC_multi_gridVectors for more details.
%                               default = 1
%
% -> 'half' (bool):             Compute only half mask.
%                               default = false
%
% -> 'lowroll' (float|str):     Strength of the gaussian roll at LOWCUT.
%                               See EMC_shapeMask for more details.
%                               If 'extended', extend the gaussian roll from the zero frequency to LOWCUT.
%                               NOTE: This is ignored if LOWCUT is nan or 0.
%                               NOTE: 0 is equivalent to no roll off.
%                               default = 0.02
%
% -> 'lowthresh'(float):        If 'lowroll' = 'extended', sets the strength of the low frequencies
%                               (lower than LOWCUT), between 0 and 1.
%                               NOTE: This is ignored if LOWCUT is nan or 0 and if 'lowroll' ~= 'extended'.
%                               default = 0
%
% -> 'highroll' (float|str):    Strength of the gaussian roll at HIGHCUT.
%                               See EMC_shapeMask for more details.
%                               If 'extended', extend the gaussian roll from HIGHCUT to Nyquist. 
%                               If HIGHCUT is at Nyquist, the default gaussian roll is used anyway.
%                               NOTE: This is ignored if HIGHCUT is nan.
%                               default = 0.02
%
%-------
% TODO 1:                       The gaussian taper is invariant to the size of the filter. It might cause
%                               an issue if the filter is very small (<100pixel wide). The current solution
%                               is to increase the low|high roll. If it is necessary, make sure the taper
%                               is at least ~7pixels.
% TODO 2:                       At the moment, for a bandpass, 2 radial grids are computed. To make it
%                               more effiecient, compute the grid one time only and adjust the low|high
%                               cut and sigma.
%
%-------
% EXAMPLE:                      [BANDPASS] = EMC_bandpass([3740,3838], 3, 10, 8, 'cpu', {});
%                               [BANDPASS] = EMC_bandpass([3740,3838], 3, 10, 8, 'cpu', ...
%                                                         {'origin',-1 ; 'lowroll','extended'});
%
% See also EMC_shapeMask.

%% checkIN
default_lowRoll = 0.02;
default_highRoll = 0.02;

[flg3d, ndim] = EMC_is3d(SIZE);

validateattributes(SIZE, {'numeric'}, {'vector', 'numel', ndim, 'integer', 'positive'}, '', 'SIZE');
validateattributes(PIXEL_SIZE, {'numeric'}, {'numel', 1, 'positive'}, '', 'PIXEL_SIZE');

OPTIONAL = EMC_extract_optional(OPTIONAL, {'origin', 'half', 'lowroll', 'lowthresh', 'highroll'});

% let EMC_multi_vectorCoordinates do the checkIN, just set the default.
if ~isfield(OPTIONAL, 'origin')
    OPTIONAL.origin = 1;  % default
end
if ~isfield(OPTIONAL.origin, 'half')
    OPTIONAL.half = false;
end

% Low frequencies
if isfield(OPTIONAL, 'lowroll')
    if isfloat(OPTIONAL.lowroll) || isinteger(OPTIONAL.lowroll)
        if OPTIONAL.lowroll < 0
            error('lowroll should be positve, got %f', OPTIONAL.lowroll)
        end
    elseif ~strcmpi(OPTIONAL.lowroll, 'extended')
        error("lowroll should be a postive float|int or 'extended', got %s", OPTIONAL.lowroll)
    end
else
    OPTIONAL.lowroll = default_lowRoll;  % default
end

if isfield(OPTIONAL, 'lowthresh')
    if ~(isfloat(OPTIONAL.lowthresh) || isinteger(OPTIONAL.lowthresh))
        error('lowthresh should be a float|int between 0 and 1, got %s', class(OPTIONAL.lowthresh))
    elseif ~(0 <= OPTIONAL.lowthresh <= 1)
        error('lowthresh should be between 0 and 1, got %f', OPTIONAL.lowthresh)
    end
else
    OPTIONAL.lowthresh = 0;  % default
end

% High frequencies
if isfield(OPTIONAL, 'highroll')
    if isfloat(OPTIONAL.highroll) || isinteger(OPTIONAL.highroll)
        if OPTIONAL.highroll < 0
            error('highroll should be positve, got %f', OPTIONAL.highroll)
        end
    elseif ~strcmpi(OPTIONAL.highroll, 'extended')
        error("highroll should be a postive float|int or 'extended', got %s", OPTIONAL.highroll)
    end
else
    OPTIONAL.highroll = default_highRoll;  % default
end
        
%% Compute the bandpass
% to make sure the resolution of the filter is isotrope, set 'isotrope'=true
[vX, vY, vZ] = EMC_multi_vectorCoordinates(SIZE, METHOD, {'origin', OPTIONAL.origin; ...
                                                          'half', OPTIONAL.half; ...
                                                          'isotrope', true});
% Normal gaussian
gaussian = @(x,s) exp(-1.*(x-1).^2 ./ (2.*s.^2));
cutoff = 0.025;  % everything below 0.025 is set to 0.
smallest_dimension = min(SIZE);

% Low cut
if isfloat(LOWCUT) || isinteger(LOWCUT)
    lowCut = smallest_dimension * PIXEL_SIZE / LOWCUT;  % [1/pix]
    
    if (flg3d)
        radius = sqrt((vX'./lowCut).^2 + (vY./lowCut).^2 + reshape(vZ./lowCut,1,1,[]).^2);
    else
        radius = sqrt((vX'./lowCut).^2 + (vY./lowCut).^2);
    end

    % Low roll
    if strcmpi(OPTIONAL.lowroll, 'extended')  % roll up to the zero frequency
        if OPTIONAL.lowthresh > 0
            sigma = sqrt(-1 / (2 * log(OPTIONAL.lowthresh)));  % set zero frequency to the input threshold
        else
            sigma = sqrt(-1 / (2 * log(cutoff)));  % set zero frequency to the cutoff
        end
    else  % classic gaussian roll
        sigma = OPTIONAL.lowroll * smallest_dimension / lowCut;
    end

    BANDPASS = (radius > 1) + (radius <= 1) .* gaussian(radius, sigma);

% No low cut
elseif isnan(LOWCUT)
    if strcmpi(METHOD, 'GPU')
        BANDPASS = ones(SIZE, 'single', 'gpuArray');
    else
        BANDPASS = ones(SIZE, 'single');
    end
end

% High cut
if isnan(HIGHCUT)
    return  % no high cut
elseif isfloat(HIGHCUT) || isinteger(HIGHCUT)
    if HIGHCUT < PIXEL_SIZE * 2
        error('HIGHCUT should be higher (lower resolution) or equal to Nyquist (%.3f)', PIXEL_SIZE * 2)
    end
    highCut = smallest_dimension * PIXEL_SIZE / HIGHCUT;  % [1/pix]
elseif strcmpi(HIGHCUT, 'nyquist')
    highCut = smallest_dimension * 0.5;
else
    error("HIGHCUT should be a float|int, nan or 'nyquist', got %s", HIGHCUT)
end

if (flg3d)
    radius = sqrt((vX'./highCut).^2 + (vY./highCut).^2 + reshape(vZ./highCut,1,1,[]).^2);
else
    radius = sqrt((vX'./highCut).^2 + (vY./highCut).^2);
end

% High roll
% Either an extended roll off (from HIGHCUT to Nyquist) or the classic roll off.
% If extended, Nyquist will have the cutoff value. Also, make sure the extended
% roll off is at least as large as the classic one.
if strcmpi(OPTIONAL.highroll, 'extended')  % rolls up to nyquist
    sigma = sqrt(-1 * ((HIGHCUT / (2*PIXEL_SIZE) - 1).^2) / (2 * log(cutoff)));
    sigma_classic = default_highRoll * smallest_dimension / highCut;
    if sigma < sigma_classic
        sigma = sigma_classic;
    end
else  % classic gaussian roll
    sigma = OPTIONAL.highroll * smallest_dimension / highCut;
end

BANDPASS = (radius < 1) .* BANDPASS + (radius >= 1) .* gaussian(radius, sigma);
BANDPASS(BANDPASS < cutoff) = 0;
clear radius
end