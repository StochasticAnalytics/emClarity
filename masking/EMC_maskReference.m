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
%     -> 'fsc' (bool):                  Is the mask meant to be used for FSC calculation?
%                                       If so, use a softer dilation for FSC masking compared to the
%                                       more stringent dilation used for alignment where the focus is
%                                       on stronger density.
%                                       default = false
%
%     -> 'lowpass' (float):             Cutoff, in Angstrom, of the low pass filter to apply
%                                       to the IMAGE before computing the reference mask.
%                                       NOTE: If lower than Nyquist, adjust it to Nyquist.
%                                       default = max(14, PIXEL*2)
%
%     -> 'threshold' (float):           The first part of the algorithm is based on a connectivity expansion
%                                       scheme. The pixels/voxels higher than a given threshold will be used
%                                       as seeds for expansion...
%
%	  -> 'com' (bool):                  Compute the center of mass of the reference mask (MASK).
%                                       default = false
%
% 	  -> 'origin' (int):                Origin of the IMAGE.
%                                       0, 1 or 2; see EMC_coordVectors for more details.
%                                       NOTE: This option is only used to compute the center of mass and
%                                             is therefore ignored if 'com'=false.
%                                       NOTE: 'origin'=-1 is not allowed as the input IMAGE should be in
%                                             real space.
%                                       defaut = 1
%
%     -> 'fraction' (bool):            	Compute the estimated particle fraction.
%                                       default = false
%
%     -> 'hydrationScaling'(float):   	Estimate of the molecular volume at the hydratation radius
%                                       of the underlying atoms. The default value (calculated empirically)
%                                       depends on the mask resolution, therefore, the 'lowpass' parameter.
%                                       NOTE: This option is only used with 'fraction'=true.
%                                       TODO: When the map resolution is lower than the masking
%                                             resolution, this will again underestimate the scaling,
%                                             artificialy *de*pressing the FSC.
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
%                                       If 'fraction'=false, return nan.
%


%% MAIN
[IMAGE, SIZE, OPTION, flg, tofloat] = checkIN(IMAGE, PIXEL, OPTION);

% Use this mask to suppress densities at edge of IMAGE.
minusEdges = EMC_shapeMask('rectangle', SIZE, ceil(SIZE/2), flg.method, ...
                           {'precision', OPTION.precision; 'origin', OPTION.origin});

% Get the IMAGE ready: median filter; force edge to go to zeros to suppress
%                      edge artifacts; lowpass filter
if flg.is3d
    IMAGE = EMC_applyBandpass(EMC_setMethod(medfilt3(gather(IMAGE), [3,3,3]), flg.method) .* minusEdges, ...
                              EMC_bandpass(SIZE, PIXEL, nan, OPTION.lowpass, flg.method, ...
                                           {'precision', OPTION.precision}));
else
    IMAGE = EMC_applyBandpass(medfilt2(IMAGE, [3,3]) .* minusEdges, ...
                              EMC_bandpass(SIZE, PIXEL, nan, OPTION.lowpass, flg.method, ...
                                           {'precision', OPTION.precision}));
end

% Make sure no wrap-around artifacts.
IMAGE = IMAGE .* minusEdges;

% Decrease progressively the threshold to select regions with lower densities.
if OPTION.fsc
    dilationThresholds = [0.9, 0.85, 0.75, 0.7, 0.65, 0.5, 0.35, 0.2, 0.1];  % loose
else
    dilationThresholds = [1.0, 0.9];  % tight
end

% Seeds: select regions of the IMAGE with highest values while making sure there is at least on pixel.
maxThreshold = max(OPTION.threshold .* std(IMAGE(IMAGE > 0), 0, 'all'), max(IMAGE, [], 'all'));
currentMask = IMAGE > maxThreshold;  % seeds

% The dilatation kernel is one of the key part of the algorithm as it restricts the selection of
% region with lowest density to regions that are in close proximity to the already selected
% regions: connectivity-based expansion.
dilationKernel = EMC_gaussianKernel([1,3], 3, METHOD, {'precision', OPTION.precision; ...
                                                       'method', flg.method});

% Connectivity-based expansion of the current mask.
for threshold = dilationThresholds .* maxThreshold
    currentMask = tofloat(currentMask);
    currentKernel = dilationKernel;

    for i = 1:ceil(threshold.^2 ./ 3); currentKernel = convn(currentKernel, currentKernel); end
    currentMask = (~currentMask .* EMC_convn(currentMask, currentKernel) > 0) + currentMask;
    
    % This part is crucial as it restricts the expansion to the close pixel higher than the current
    % threshold. This is where the 'connectivity' really happens.
    currentMask = (IMAGE .* currentMask > threshold);
end

% At this stage, the currentMask is considered as the particle volume (non-hydrated).
% Same this volume to compute the particle fraction.
if OPTION.fraction
    particleVolEstimate = sum(currentMask, 'all');
end

% Expand the currentMask (particle volume) by at least 10 Angstroms.
rad = max(3, floor(10 ./ PIXEL));
for i = 1:size(currentMask,3)
    currentMask(:,:,i) = currentMask(:,:,i) + bwdist(currentMask(:,:,i)) < (rad / 2);
end

% For FSC masks, apply a gaussian blur to the mask.
if OPTION.fsc
    smoothKernel = EMC_gaussianKernel([1, rad], rad/2, {'precision', OPTION.precision; 'method', flg.method});
    currentMask = EMC_convn(tofloat(currentMask), convn(smoothKernel, smoothKernel));
    currentMask = currentMask ./ max(currentMask(:));
end

% Make sure no wrap-around artifacts.
MASK = currentMask .* minusEdges;

% Compute the center of mass of the final MASK.
% Q4B8: If fsc mask, why aren't you computing the com on the final MASK? line 458 of BH_mask3d is weird:
% binaryVol = (binaryVol - min(binaryVol(MASK > 0.01)) ).*MASK; with binaryVol being the binary mask and
% MASK being binaryVol with blurring and edgetaper.
if OPTION.com
    COM = EMC_centerOfMass(MASK, OPTION.origin);
else
    COM = nan;
end

% Estimate the particle fraction
if OPTION.fraction
    % The mask is not a logical mask, it has values between 0 and 1. As such,
    % estimate the signal reduction in these regions.
    powerReduction = sum(IMAGE.^2 .* currentMask>0, 'all') ./ ...
                     sum(IMAGE.^2 .* currentMask, 'all');

    maskVolume = sum(currentMask>0, 'all');

    % Scale the particle volume to add its hydration volume.
    particleVolEstimate = particleVolEstimate ./ OPTION.hydrationScaling;

    % Finally, estimate the fraction of the mask taken by the particle by
    % comparing the hydrated particle volume to the mask volume (scalled down
    % to take into account the power reduction due to the roll off)
    FRACTION = particleVolEstimate ./ (maskVolume .* powerReduction);

    fprintf(['Size: %s - Precision: %s - Method: %s\n', ...
             'Estimated particule volume : %d voxels\n', ...
             'Estimated mask volume      : %d voxels\n', ...
             'Power reduction            : %2.3f\n', ...
             'Particle fraction          : %2.3f\n'], ...
             mat2str(SIZE), OPTION.precision, flg.method, ...
             particleVolEstimate, maskVolume, powerReduction, FRACTION);
else
    FRACTION = nan;
end

end


function [IMAGE, SIZE, OPTION, flg, tofloat] = checkIN(IMAGE, PIXEL, OPTION)

if ~isnumeric(IMAGE)
    error('EMC:IMAGE', 'IMAGE should be numeric, got %s', class(IMAGE))
elseif ~isreal(IMAGE)
    error('EMC:IMAGE', 'IMAGE should be real, got complex')
end
[flg.is3d, SIZE, flg.ndim] = EMC_is3d(size(IMAGE));

if EMC_isOnGpu(IMAGE)
    flg.method = 'gpu';
else
    flg.method = 'cpu';
end

% PIXEL
if ~isscalar(PIXEL) || ~isnumeric(PIXEL) || isinf(PIXEL) || ~(PIXEL > 0)
    error('EMC:PIXEL', 'PIXEL should be a positive float|int')
end

OPTION = EMC_getOption(OPTION, {'origin', 'fsc', 'com', 'lowpass', ...
                                'threshold', 'fraction', 'hydrationScaling', 'precision'}, false);

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

% fraction
if isfield(OPTION, 'fraction')
    if  ~isscalar(OPTION.fraction) || ~islogical(OPTION.fraction)
        error('EMC:fraction', 'OPTION.fraction should be a boolean, got %s', class(OPTION.fraction))
    end
else
    OPTION.fraction = false;  % default
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

% hydrationScaling
if isfield(OPTION, 'hydrationScaling')
    if ~isscalar(OPTION.hydrationScaling) || ~isnumeric(OPTION.hydrationScaling) || ...
       isinf(OPTION.hydrationScaling) || ~(OPTION.hydrationScaling >= 0)
        error('EMC:hydrationScaling', 'OPTION.hydrationScaling should be a nonnegative float|int')
    end
else
    OPTION.hydrationScaling = (-2.8e-3) .* OPTION.lowpass.^2 + 0.14 .* OPTION.lowpass + 1.5;  % default
end

% precision
if isfield(OPTION, precision)
    if strcmpi(OPTION.precision, 'single')
        tofloat = @(array) single(array);
    elseif strcmpi(OPTION.precision, 'double')
        tofloat = @(array) double(array);
    else
        error('EMC:precision', "OPTION.precision should be 'single' or 'double'")
    end
else
    OPTION.precision = 'single';  % default
    tofloat = @(array) single(array);
end

end  % checkIN
