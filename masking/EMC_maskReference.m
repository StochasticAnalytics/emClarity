function [MASK, varargout] = EMC_maskReference(IMAGE, PIXEL, OPTION)
%
% [MASK, COM]                      = EMC_maskReference(IMAGE, PIXEL, OPTION) if 'fsc'=false.
% [MASK, MASK_CORE, FRACTION, COM] = EMC_maskReference(IMAGE, PIXEL, OPTION) if 'fsc'=true.
% Compute a real-space reference mask of a 3d/2d particle.
%
% Input:
%   IMAGE (numeric):                    2d/3d (real space) image to make a reference mask of.
%
%   PIXEL (float):                  	Pixel size (A/pix) of the IMAGE.
%
%   OPTION (cell | struct):             Optional parameters.
%                                       If cell: {field, value; ...}, note the ';' between parameters.
%                                       NOTE: Can be empty.
%                                       NOTE: Unknown fields will raise an error.
%
%     -> 'lowpass' (float):             Cutoff, in Angstrom, of the low pass filter to apply
%                                       to the IMAGE before computing the reference mask.
%                                       NOTE: If lower than Nyquist, adjust it to Nyquist.
%                                       default = 14
%
%     -> 'threshold' (float):           The first part of the algorithm is based on a connectivity expansion
%                                       scheme (dilation weigthed by the density strengh). The pixels/voxels
%                                       higher than a given threshold will be used as seeds for expansion.
%                                       This threshold is empirical and applied on the IMAGE after
%                                       preprocessing (median filter, lowpass, taper and standardization).
%                                       default = 2.4
%
%     -> 'com' (bool):                  Compute the center of mass of the reference mask (MASK).
%                                       default = false
%
%     -> 'origin' (int):                Origin of the IMAGE.
%                                       0, 1 or 2; see EMC_coordVectors for more details.
%                                       NOTE: This option is only used to compute the center of mass and
%                                             is therefore ignored if 'com'=false.
%                                       NOTE: 'origin'=-1 is not allowed as the input IMAGE should be in
%                                             real space.
%                                       defaut = 1
%
%     -> 'fsc' (bool):                  Whether or not the mask is meant for FSC calculation?
%                                       If true, use a stronger dilation for FSC masking compared to the
%                                       more stringent dilation used for alignment. Additionnaly, computes
%                                       the particle fraction (FRACTION).
%                                       default = false
%
%     -> 'hydration_scaling'(float):    Estimate of the molecular volume at the hydratation radius
%                                       of the underlying atoms. The default value (calculated empirically)
%                                       depends on the mask resolution, therefore, the 'lowpass' parameter.
%                                       NOTE: This option is only used with 'fsc'=true.
%                                       default = see checkIN.
%
% Output:
%   MASK (numeric):                     Reference mask. Its method and precision correspond to IMAGE.
%
%   varargout:                          If 'fsc' is true, varargout = {MASKCORE, FRACTION, COM}.
%                                       If 'fsc' is false, varargout = {COM}.
%
%                                       -> MASKCORE (numeric):  Core of the reference mask used
%                                                               for a tight fsc calculation, as
%                                                               in RELION, or as visual inspection.
%
%                                       -> FRACTION (float):    Estimated particle/background ratio
%                                                               in the MASK.
%
%                                       -> COM (vector):        Center of mass of the 2d/3d reference
%                                                               mask. If 'com'=false, return nan.
%
% Note:
%   - I [TF] don't understand this: "When the map resolution is lower than the masking
%     resolution, this will underestimate the scaling, artificialy *de*pressing the FSC".
%
%   - The mask is very dependent of OPTION.threshold, specially for alignment masks.
%     I [TF] am supprised this hasn't cause problems.
%
% Other EMC-files required:
%   EMC_getOption, EMC_is3d, EMC_resize, EMC_applyBandpass, EMC_getBandpass, EMC_mediaFilter,
%   EMC_gaussianKernel, EMC_convn, EMC_centerOfMass.
%

% Created:  18Jan2020, R2019a
% Version:  v.1.0. unittest; fix dilation; fix lowpass and COM with fsc=true;
%                  add MASK_CORE (TF, 9Mar2020).
%

% Global verbosity.
global EMC_gp_verbose
if isempty(EMC_gp_verbose)
    EMC_gp_verbose = true;
end

%% MAIN
[SIZE, OPTION, flg, tofloat] = checkIN(IMAGE, PIXEL, OPTION);

% Use this mask to suppress densities at the edges of IMAGE.
if flg.isOnGpu
    minusEdges = EMC_resize(ones(SIZE, flg.precision, 'gpuArray'), nan, {'force_taper', true});
else
    minusEdges = EMC_resize(ones(SIZE, flg.precision), nan, {'force_taper', true});
end

% Get the IMAGE ready: median filter (remove hot pixels), force edge to go to
% zeros to suppress edge artifacts and lowpass filter. For alignment masks, apply 'soft' lowpass.
if OPTION.fsc
    optBandpass = {'precision', flg.precision};  % default lowpassRoll
else
    optBandpass = {'precision', flg.precision; 'lowpassRoll', 'extended'};
end

% To have identical outputs with BH_mask3d, replace the bandpass with BH_bandpass3d.
if flg.is3d
    if flg.isOnGpu
        IMAGE = EMC_applyBandpass(...
                    gpuArray(medfilt3(gather(IMAGE), [3,3,3])) .* minusEdges, ...
                    EMC_getBandpass(SIZE, PIXEL, nan, OPTION.lowpass, flg.method, optBandpass), ...
                    {});
    else
        IMAGE = EMC_applyBandpass(...
                    medfilt3(IMAGE, [3,3,3]) .* minusEdges, ...
                    EMC_getBandpass(SIZE, PIXEL, nan, OPTION.lowpass, flg.method, optBandpass), ...
                    {});
    end
else
    IMAGE = EMC_applyBandpass(...
                medfilt2(IMAGE, [3,3]) .* minusEdges, ...
                EMC_getBandpass(SIZE, PIXEL, nan, OPTION.lowpass, flg.method, optBandpass), ...
                {});
end

% Make sure no wrap-around artifacts.
IMAGE = IMAGE .* minusEdges;

% Decrease progressively the threshold to select regions with lower densities.
if OPTION.fsc
    dilationThresholds = [0.9, 0.85, 0.75, 0.7, 0.65, 0.5, 0.35, 0.2, 0.1];  % loose
else
    dilationThresholds = [1.0, 0.9];  % tight
end

% Seeds: select regions of the IMAGE with highest values while making sure there is at least one pixel.
maxThreshold = min(OPTION.threshold .* std(IMAGE(IMAGE > 0), [], 'all'), max(IMAGE, [], 'all'));
binaryMask = IMAGE > maxThreshold;  % seeds

% The dilatation kernel is one of the key part of the algorithm as it restricts the selection of
% region with lowest density to regions that are in close proximity to the already selected
% regions: connectivity-based expansion.
dilationKernel = EMC_gaussianKernel([1,3], 3, flg.method, {'precision', flg.precision});

% Connectivity-based expansion of the current mask.
for threshold = dilationThresholds .* maxThreshold
    for i = 1:ceil(threshold.^2 ./ 3)
        binaryMask = (~binaryMask .* EMC_convn(tofloat(binaryMask), dilationKernel) > 0.00) + binaryMask;
    end

    % This part is crucial as it restricts the expansion to the close pixel higher than the current
    % threshold. This is where the 'connectivity' really happens.
    binaryMask = (IMAGE .* binaryMask) > threshold;
end

% At this stage, the currentMask is considered as the particle volume (non-hydrated).
% Save this volume to compute the particle fraction.
if OPTION.fsc
    particleVolEstimate = sum(binaryMask, 'all');

    taperKernel = EMC_gaussianKernel([1,4], 1.75, flg.method, {'precision', flg.precision});
    MASK_CORE = EMC_convn(tofloat(binaryMask), taperKernel);
    MASK_CORE = EMC_convn(sqrt(MASK_CORE),taperKernel);
end

% Expand the currentMask (particle volume) by at least 10 Angstroms.
rad = max(3, floor(10 ./ PIXEL));
for i = 1:size(binaryMask, 3)
    binaryMask(:,:,i) = binaryMask(:,:,i) + bwdist(binaryMask(:,:,i)) < (rad / 2);
end

% At this point, currentMask is probably a slight over estimate of the protein envelope, however,
% we want to leave room for regions which aren't yet aligned well and thus appear weak. Otherwise,
% the centering on just the strong density may prevent them from ever improving.

% For FSC masks, add an additional dilation. By expanding the mask to include more solvent,
% one can arbitrarily decrease the measured SNR. Even though this expansion is almost certainly
% containing surrounding solvent, estimate the signal reduction in the taper as the mask could
% cut through densities.
if OPTION.fsc
    smoothKernel = EMC_gaussianKernel([1, rad], rad/2, flg.method, {'precision', flg.precision});
    fscMask = EMC_convn(tofloat(binaryMask), convn(smoothKernel, smoothKernel));
    fscMask = fscMask ./ max(fscMask(:));
    maskVolume = sum(fscMask>0, 'all');

    powerReduction = sum(IMAGE.^2 .* (fscMask>0), 'all') ./ sum(IMAGE.^2 .* fscMask, 'all');

    % Scale the particle volume; 'remove' its hydration volume.
    % TODO: so we assume particleVolEstimate contains some solvent then?
    particleVolEstimate = particleVolEstimate ./ OPTION.hydration_scaling;

    % Finally, estimate the fraction of the mask taken by the particle, by
    % comparing the particle volume (not-hydrated) to the mask volume (scalled down
    % to take into account the power reduction due to the taper).
    FRACTION = maskVolume ./ (particleVolEstimate .* powerReduction);

    if EMC_gp_verbose
        fprintf(['FSC mask: Estimated particule volume : %d voxels\n', ...
                 '               Estimated mask volume : %d voxels\n', ...
                 '                     Power reduction : %2.3f\n', ...
                 '                   Particle fraction : %2.3f\n'], ...
                 particleVolEstimate, maskVolume, powerReduction, 1/FRACTION);
    end

    % Make sure no wrap-around artifacts.
    MASK = fscMask .* minusEdges;

    % To compute the COM, restrict to the region most likely to contain only the protein.
    if OPTION.com; COM = EMC_centerOfMass(MASK .* binaryMask, OPTION.origin); else; COM = nan; end

    % Setting varargout.
    varargout = {MASK_CORE, FRACTION, COM};
else
    % Make sure no wrap-around artifacts.
    MASK = binaryMask .* minusEdges;

    if OPTION.com; COM = EMC_centerOfMass(MASK, OPTION.origin); else; COM = nan; end

    % Setting varargout.
    varargout = {COM};
end

end  % EMC_maskReference


function [SIZE, OPTION, flg, tofloat] = checkIN(IMAGE, PIXEL, OPTION)

if ~isnumeric(IMAGE)
    error('EMC:IMAGE', 'IMAGE should be numeric, got %s', class(IMAGE))
elseif ~isreal(IMAGE)
    error('EMC:IMAGE', 'IMAGE should be real, got complex')
end
[flg.is3d, SIZE, ~] = EMC_is3d(size(IMAGE));


% PIXEL
if ~isscalar(PIXEL) || ~isnumeric(PIXEL) || isinf(PIXEL) || ~(PIXEL > 0)
    error('EMC:PIXEL', 'PIXEL should be a positive float|int')
end

OPTION = EMC_getOption(OPTION, {'origin', 'fsc', 'com', 'lowpass', 'threshold', ...
                                'hydration_scaling', 'precision', 'method'}, false);

% origin
if isfield(OPTION, 'origin')
    if ~isscalar(OPTION.origin) || ~isnumeric(OPTION.origin)
        error('EMC:origin', 'OPTION.origin should be an integer, got %s of size: %s', ...
              class(OPTION.origin), mat2str(size(OPTION.origin)))
    elseif OPTION.origin ~= 1 && OPTION.origin ~= 0 && OPTION.origin ~= 2
        error('EMC:origin', 'OPTION.origin should be 0, 1 or 2, got %d', OPTION.origin)
    end
else
    OPTION.origin = 1;  % default
end

% fsc
if isfield(OPTION, 'fsc')
    if ~isscalar(OPTION.fsc) || ~islogical(OPTION.fsc)
        error('EMC:fsc', 'OPTION.fsc should be a boolean, got %s', class(OPTION.fsc))
    end
else
    OPTION.fsc = false;  % default
end

% com
if isfield(OPTION, 'com')
    if ~isscalar(OPTION.com) || ~islogical(OPTION.com)
        error('EMC:com', 'OPTION.com should be a boolean, got %s', class(OPTION.com))
    end
else
    OPTION.com = false;  % default
end

% lowpass
if isfield(OPTION, 'lowpass')
    if ~isscalar(OPTION.lowpass) || ~isnumeric(OPTION.lowpass) || ...
       isinf(OPTION.lowpass) || ~(OPTION.lowpass >= 0)
        error('EMC:LOWPASS', 'OPTION.lowpass should be a nonnegative float|int')
    elseif OPTION.lowpass < PIXEL * 2
        OPTION.lowpass = PIXEL * 2;
    end
else
    OPTION.lowpass = max(14, PIXEL * 2);  % default
end

% threshold
if isfield(OPTION, 'threshold')
    if ~isscalar(OPTION.threshold) || ~isnumeric(OPTION.threshold)
        error('EMC:threshold', 'OPTION.threshold should be a numeric scalar, got %s of size: %s', ...
              class(OPTION.threshold), mat2str(size(OPTION.threshold)))
    end
else
    OPTION.threshold = 2.5;  % default
end

% hydration_scaling
if isfield(OPTION, 'hydration_scaling')
    if ~isscalar(OPTION.hydration_scaling) || ~isnumeric(OPTION.hydration_scaling) || ...
       isinf(OPTION.hydration_scaling) || ~(OPTION.hydration_scaling >= 0)
        error('EMC:hydration_scaling', 'OPTION.hydration_scaling should be a nonnegative float|int')
    end
else
    OPTION.hydration_scaling = (-2.8e-3) .* OPTION.lowpass.^2 + 0.14 .* OPTION.lowpass + 1.5;  % default
end

[flg.precision, flg.isOnGpu, flg.method] = EMC_getClass(IMAGE);

% precision
if strcmpi(flg.precision, 'single')
    tofloat = @(array) single(array);
elseif strcmpi(flg.precision, 'double')
    tofloat = @(array) double(array);
else
    error('EMC:precision', "IMAGE should be 'single' or 'double'")
end

end  % checkIN
