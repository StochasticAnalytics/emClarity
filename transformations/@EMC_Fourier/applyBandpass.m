function ARRAY = applyBandpass(obj, ARRAY, OPTION)
%
% ARRAY = obj.applyBandpass(ARRAY, OPTION)
% Apply the current bandpass filter (obj.bandpass) to an ARRAY.
%
% Input:
%   ARRAY (numeric):           	2d/3d array to filter.
%
%   OPTION (cell|struct):       Optional parameters.
%                              	If cell: {field, value; ...}, note the ';' between parameters.
%                              	NOTE: Can be empty.
%                             	NOTE: Unknown fields will raise an error.
%
%     -> 'fft' (bool):          Wheter or not the ARRAY is in real space.
%                               If true:  compute the fft of ARRAY before applying the bandpass.
%                               If false: apply the bandpass directly to the ARRAY. The size of the
%                                         ARRAY should match the size of obj.bandpass.
%                               default = true
%
%     -> 'ifft' (bool):         Wheter or not the ARRAY should be switched back to real space after
%                               applying the bandpass. If false, the ARRAY is left in frequency space.
%                              	default = true
%
%     -> 'standardize' (bool):  Set the real space mean to ~0 and real space variance to ~1 of the ARRAY.
%                             	NOTE: the mean and variance are changed in frequency space.
%                             	default = true
%
% Output:
%   ARRAY (numeric):            Filtered ARRAY.
%
% Property used:
%   obj.bandpass
%
% Method used:
%   obj.fft
%   obj.ifft
%   obj.standardize
%
% Examples:
%   - >> ft = EMC_Fourier([3740,3838], 'gpu', {});
%     >> ft = ft.setBandpass(3, nan, 8, {});  % lowpass filter with cutoff at 8A
%     >> img = ft.applyBandpass(img, {});
%
% Other EMC-files required:
%   EMC_getOption
%
% See also obj.setBandpass.
%

% Created:  18Jan2020, R2019a
% Version:  v.1.0.  symplify the inputs - use EMC_resize directly if you want to apply
%                   a taper and/or pad|crop before applying the bandpass (TF, 2Feb2020).
%           v.1.1.  explicit checks for IMAGE and BANDPASS; unittest (TF, 3Feb2020).
%           v.1.2.  half bandpass are now supported (TF, 8Feb2020).
%           v.1.2.1 bug fix: standardization with half grid is now properly done. (TF, 22Feb2020).
%           v.1.3.  EMC_applyBandpass is integrated into EMC_Fourier.applyBandpass; can now used
%                   ARRAYs already in Fourier space (OPTION.fft) (TF, 8Mar2020).
%

%% checkIN
OPTION = EMC_getOption(OPTION, {'fft', 'ifft', 'standardize'}, false);

% isfreq
if isfield(OPTION, 'fft')
    if ~islogical(OPTION.fft) && ~isscalar(OPTION.fft)
        error('EMC:Fourier', 'OPTION.fft should be a boolean')
    end
else
    OPTION.fft = true;  % default
end

% ifft
if isfield(OPTION, 'ifft')
    if ~islogical(OPTION.ifft) && ~isscalar(OPTION.ifft)
        error('EMC:Fourier', 'OPTION.ifft should be a boolean')
    end
else
    OPTION.ifft = true;  % default
end

% standardize
if isfield(OPTION, 'standardize')
    if ~islogical(OPTION.standardize) && ~isscalar(OPTION.standardize)
        error('EMC:Fourier', 'OPTION.standardize should be a boolean')
    end
else
    OPTION.standardize = true;  % default
end

%% Apply filter.
if OPTION.fft
    ARRAY = obj.fft(ARRAY) .* obj.bandpass;
else
    ARRAY = ARRAY .* obj.bandpass;
end

if OPTION.standardize
    ARRAY = obj.standardize(ARRAY);
end

if OPTION.ifft
    ARRAY = obj.ifft(ARRAY);
end

end  % applyBandpass
