function [MASK, COM, FRACTION] = EMC_maskReference(IMAGE, PIXEL, OPTION)
%
% [MASK, COM, FRACTION] = EMC_maskReference(IMAGE, PIXEL, OPTION)
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
%                                       default = see line 266.
%
%     -> 'precision' (str):             Precision of the output MASK.
%                                       NOTE: This is changed before lowpass filter.
%                                       default = same as input IMAGE.
%
% Output:
%   MASK (numeric):                     Reference mask.
%
%   COM (vector):                       Center of mass of the 2d/3d reference mask.
%                                       If 'com'=false, return nan.
%
%   FRACTION (float):                   Estimated particle/background ratio in the MASK.
%                                       If 'fsc'=false, return nan.
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
%   EMC_gaussianKernel, EMC_setPrecision, EMC_convn, EMC_centerOfMass.
%

% Created:  18Jan2020, R2019a
% Version:  v.1.0.
%

%% MAIN
[IMAGE, SIZE, METHOD, OPTION] = checkIN(IMAGE, PIXEL, OPTION);

% Use this mask to suppress densities at the edges of IMAGE.
if strcmp(METHOD, 'gpu')
    minusEdges = EMC_resize(ones(SIZE, OPTION.precision, 'gpuArray'), nan, {'force_taper', true});
else
    minusEdges = EMC_resize(ones(SIZE, OPTION.precision), nan, {'force_taper', true});
end

% Get the IMAGE ready: median filter (remove hot pixels); force edge to go to zeros to suppress
% edge artifacts; lowpass filter
IMAGE = EMC_applyBandpass(...
            EMC_medianFilter(IMAGE) .* minusEdges, ...
            EMC_getBandpass(SIZE, PIXEL, nan, OPTION.lowpass, METHOD, {'precision', OPTION.precision; ...
                                                                       'lowpassRoll', 'extended'}), ...
            {});

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
currentMask = IMAGE > maxThreshold;  % seeds

% The dilatation kernel is one of the key part of the algorithm as it restricts the selection of
% region with lowest density to regions that are in close proximity to the already selected
% regions: connectivity-based expansion.
dilationKernel = EMC_gaussianKernel([1,3], 3, {'precision', OPTION.precision; 'method', METHOD});

% Connectivity-based expansion of the current mask.
for threshold = dilationThresholds .* maxThreshold
    currentMask = EMC_setPrecision(currentMask, OPTION.precision);
    currentKernel = dilationKernel;

    for i = 1:ceil(threshold.^2 ./ 3) - 1
        currentKernel = convn(currentKernel, currentKernel);
    end
    currentMask = (~currentMask .* EMC_convn(currentMask, currentKernel) > 0) + currentMask;

    % This part is crucial as it restricts the expansion to the close pixel higher than the current
    % threshold. This is where the 'connectivity' really happens.
    currentMask = (IMAGE .* currentMask > threshold);
end

% At this stage, the currentMask is considered as the particle volume (non-hydrated).
% Save this volume to compute the particle fraction.
if OPTION.fsc
    particleVolEstimate = sum(currentMask, 'all');
end

% Expand the currentMask (particle volume) by at least 10 Angstroms.
rad = max(3, floor(10 ./ PIXEL));
for i = 1:size(currentMask, 3)
    currentMask(:,:,i) = currentMask(:,:,i) + bwdist(currentMask(:,:,i)) < (rad / 2);
end
currentMask = EMC_setPrecision(currentMask, OPTION.precision);

% At this point, currentMask is probably a slight over estimate of the protein envelope, however,
% we want to leave room for regions which aren't yet aligned well and thus appear weak. Otherwise,
% the centering on just the strong density may prevent them from ever improving.

% For FSC masks, add an additional dilation. By expanding the mask to include more solvent,
% one can arbitrarily decrease the measured SNR. Even though this expansion is almost certainly
% containing surrounding solvent, estimate the signal reduction in the taper as the mask could
% cut through densities.
if OPTION.fsc
    smoothKernel = EMC_gaussianKernel([1, rad], rad/2, {'precision', OPTION.precision; 'method', METHOD});
    fscMask = EMC_convn(currentMask, convn(smoothKernel, smoothKernel));
    fscMask = fscMask ./ max(fscMask(:));
    maskVolume = sum(fscMask>0, 'all');

    powerReduction = sum(IMAGE.^2 .* fscMask>0, 'all') ./ sum(IMAGE.^2 .* fscMask, 'all');

    % Scale the particle volume; 'remove' its hydration volume.
    % TODO: so we assume particleVolEstimate contains some solvent then?
    particleVolEstimate = particleVolEstimate ./ OPTION.hydration_scaling;

    % Finally, estimate the fraction of the mask taken by the particle, by
    % comparing the particle volume (not-hydrated) to the mask volume (scalled down
    % to take into account the power reduction due to the taper).
    FRACTION = particleVolEstimate ./ (maskVolume .* powerReduction);

    fprintf(['FSC mask: Estimated particule volume : %d voxels\n', ...
             '               Estimated mask volume : %d voxels\n', ...
             '                     Power reduction : %2.3f\n', ...
             '                   Particle fraction : %2.3f\n'], ...
             particleVolEstimate, maskVolume, powerReduction, FRACTION);

    % Make sure no wrap-around artifacts.
    MASK = fscMask .* minusEdges;

    % To compute the COM, restrict to the region most likely to contain only the protein.
    if OPTION.com; COM = EMC_centerOfMass(MASK .* currentMask, OPTION.origin); else; COM = nan; end
else
    % Make sure no wrap-around artifacts.
    MASK = currentMask .* minusEdges;

    if OPTION.com; COM = EMC_centerOfMass(MASK, OPTION.origin); else; COM = nan; end
    FRACTION = nan;
end

end  % EMC_maskReference


function [IMAGE, SIZE, METHOD, OPTION] = checkIN(IMAGE, PIXEL, OPTION)

if ~isnumeric(IMAGE)
    error('EMC:IMAGE', 'IMAGE should be numeric, got %s', class(IMAGE))
elseif ~isreal(IMAGE)
    error('EMC:IMAGE', 'IMAGE should be real, got complex')
end
[~, SIZE, ~] = EMC_is3d(size(IMAGE));

if EMC_isOnGpu(IMAGE)
    METHOD = 'gpu';
else
    METHOD = 'cpu';
end

% PIXEL
if ~isscalar(PIXEL) || ~isnumeric(PIXEL) || isinf(PIXEL) || ~(PIXEL > 0)
    error('EMC:PIXEL', 'PIXEL should be a positive float|int')
end

OPTION = EMC_getOption(OPTION, {'origin', 'fsc', 'com', 'lowpass', ...
                                'threshold', 'hydration_scaling', 'precision'}, false);

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

% precision
if isfield(OPTION, 'precision')
    if ~(ischar(OPTION.precision) || isstring(OPTION.precision)) || ...
       ~strcmpi(OPTION.precision, 'single') && ~strcmpi(OPTION.precision, 'double')
        error('EMC:precision', "OPTION.precision should be 'single' or 'double'")
    end
else
    OPTION.precision = 'single';  % default
end

end  % checkIN
