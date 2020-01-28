function [MASK, COM, PARTICLE_FRACTION] = EMC_referenceMask(IMAGE, PIXEL, OPTION)
%
% [MASK, COM, PARTICLE_FRACTION] = EMC_referenceMask(IMAGE, PIXEL, OPTION)
% Compute a real-space reference mask of a 3d/2d particle.
%
% IMAGE (numeric):                  2d/3d (real space) image to make a reference mask of.
%
% PIXEL (float):                    Pixel size (A/pix) of the IMAGE.
%
% METHOD (str):                     Device to use; 'gpu' or 'cpu'
%
% OPTION (cell | struct):           Optional parameters.
%                                   If cell: {field, value; ...}, note the ';' between parameters.
%                                   NOTE: Can be empty.
%                                   NOTE: Unknown fields will raise an error.
%
% 	-> 'fsc' (bool):                Is the mask meant to be used for FSC calculation?
%                                   If so, use a softer dilation for FSC masking compared to the more
%                                   stringent dilation used for alignment where the focus is
%                                   on stronger density.
%                                   default = false
%
% 	-> 'lowpass' (float):           Low pass filter cut off (in Angtrom) to apply to the IMAGE before
%                                   computing the reference mask.
%                                   default = 14
%                                   Q4B1: 14A is quite arbitrary? By default I was thinking to define the cut
%                                   as a fraction of nyquist, with maybe a minimun value... I don't have the
%                                   full picture here yet, so I don't know if it is a good idea.
%
%  	-> 'threshold' (float):         The first part of the algorithm is based on a connectivity expansion
%                                   scheme. The pixels/voxels higher than a given threshold will be used
%                                   as seeds for expansion...
%
%                                   Q4B2: Depending on the particle/background ratio, could the number of
%                                   initial pixel vary significantly?
%                                   I wanted to change it to something like, "select the 2% highest pixels
%                                   as seed" but without sorting the pixel values, I dont know how to do
%                                   this because the distribution can vary... Do we want to select a
%                                   constant percentage of pixel though? If the particle is 90% of the volume,
%                                   we might want to select more pixels than if the particle is 40% of the
%                                   volume...? Maybe the fact to setting the mean to 0 and having the
%                                   threshold relative to the std of the positive pixel is a way of doing
%                                   this. I guess you already thought about this?
%
%                                   Q4B3: Can we consider the noise to be gaussian in cryoEM?
%
%  	-> 'com' (bool):                Compute the center of mass of the reference mask (MASK).
%                                   default = false
%
%	-> 'origin' (int):              Origin of the IMAGE.
%                                   0, 1 or 2; see EMC_multi_gridVectors for more details.
%                                   NOTE: This option is only used to compute the center of mass and
%                                         is therefore ignored if 'com'=false.
%                                   NOTE: 'origin'=-1 is not allowed as the input IMAGE should be in
%                                         real space.
%                                   defaut = 1
%
%   -> 'particle_fraction' (bool):  Compute the estimated particle fraction.
%                                   default = false
%
% 	-> 'hydration_scaling'(float):  Estimate of the molecular volume at the hydratation radius
%                                   of the underlying atoms. The default value (calculated empirically)
%                                   depends on the mask resolution, therefore, the 'lowpass' parameter.
%                                   NOTE: This option is only used with 'particle_fraction'=true.
%                                   TODO: When the map resolution is lower than the masking
%                                         resolution, this will again underestimate the scaling,
%                                         artificialy *de*pressing the FSC.
%
%   -> 'precision' (str):           Precision of the output MASK.
%                                   NOTE: This is changed before lowpass filter.
%                                   default = same as input IMAGE.
%
%   -> 'verbose' (int):             Verbosity; 0: none, 1: default, 2:debug
%                                   default = 1
%
%---------
% RETURN:                           MASK: Reference mask.
%                                   COM:  Center of mass of the 2d/3d reference mask, if asked, else, nan.
%                                   PARTICLE_FRACTION: Estimated particle/background ratio in the MASK.
%
%---------
% EXAMPLE:                          [MASK, ~] = EMC_referenceMask(IMAGE, 1.4, {})
%                                   [MASK, ~] = EMC_referenceMask(IMAGE, 1.4, ...
%                                                                 {'low_pass', 5 ; 'com', true})

%% MAIN
[IMAGE, SIZE, OPTION, flg, ndim, tofloat] = checkIN(IMAGE, PIXEL, METHOD, OPTION);

% Use this mask to taper the edges of IMAGE.
taperEdge = EMC_shapeMask('rectangle', SIZE, ceil(SIZE/2), METHOD, {'precision', OPTION.precision});

% Lowpass filter only if cutoff < nyquist.
%Q4B4: What's the point of doing a lowpass with cutoff after Nyquist?
if flg.lowpass
    lowpass = EMC_bandpass(SIZE, PIXEL, nan, OPTION.lowpass, METHOD, {'precision', OPTION.precision});
else
    lowpass = nan;
end

% Median filter, lowpass and taper the edges to 0.
if (flg.is3d)
    % medfilt3 asks for a lot of GPU memory, so do it on the cpu and transfer it back
    % automatically when multiplying with taperEdge.
    IMAGE = EMC_FilterCenterNormalize(medfilt3(gather(IMAGE), [3,3,3]) .* taperEdge, lowpass, {});
else
    IMAGE = EMC_FilterCenterNormalize(medfilt2(IMAGE, [3,3]) .* taperEdge, lowpass, {});
end

% Apply taper to the edges of the IMAGE also after lowpass filter.
% Q4B5: Wouldn't it be enough to taper only after lowpass?
IMAGE = IMAGE .* taperEdge;

% Decrease progressively the threshold to select regions with lower densities.
if (OPTION.fsc)
    dilationThresholds = [0.9, 0.85, 0.75, 0.7, 0.65, 0.5, 0.35, 0.2, 0.1];  % loose
else
    dilationThresholds = [1.0, 0.9];  % tight
end

% The dilatation kernel is one of the key part of the algorithm as it restricts the selection of
% region with lowest density to regions that are in close proximity to the already selected
% regions: connectivity-based expansion.
dilationKernel = EMC_multi_gaussianKernel(3.*ones(1,ndim), 3.0, METHOD, {'precision', OPTION.precision});

% Seeds: select regions of the IMAGE with highest values.
maxThreshold = OPTION.threshold .* (std(IMAGE(IMAGE(:) > 0)));
maxAllowed = max(IMAGE(:));
if maxThreshold > maxAllowed; maxThreshold = maxAllowed; end  % make sure there is at last one pixel.
currentMask = IMAGE > maxThreshold;  % seeds

% Connectivity-based expansion of the current mask.
for threshold = dilationThresholds .* maxThreshold
    currentMask = tofloat(currentMask);

    % Q4B6: Why aren't you using a bigger kernel and do the convn only once?
    dilationIter = ceil(threshold.^2 ./ 3);
    for i = 1:dilationIter
        currentMask = (~currentMask .* convn(currentMask, dilationKernel, 'same') > 0) + currentMask;
    end
    % This part is crucial as it restricts the expansion to the close pixel higher than the current
    % threshold. This is where the 'connectivity' really happens.
    currentMask = (IMAGE .* currentMask > threshold);
end

% At this stage, the currentMask is considered as the particle volume (non-hydrated)
% when calculating the particle fraction.
if (OPTION.particle_fraction)
    particleVolEstimate = sum(currentMask(:));  % check if sum(currentMask(:) == 1) is faster.
end

% The connectivity-based expansion will limit the expansion to the pixels that are higher
% than the current threshold. As such, it will not necessary be very smooth (a negative pixel in
% the middle of a selected region will not be selected for example). This next part tries to prevent
% that by expending the mask by ~1.5 pixels around the selected region (only based on Euclidian
% distance and not on pixel value).
rad = max(3, floor(10 ./ PIXEL));
% Q4B7: Why do you want a large expansion if the pixel size is smaller?

for i = 1:size(currentMask,3)
    euclidianDistance = bwdist(currentMask(:,:,i));  % bwdist is only 2D, so loop through volume.
    currentMask(:,:,i) = currentMask(:,:,i) + euclidianDistance < (rad / 2);
end
clear euclidianDistance
% Q4B8: To be honest, in practive I don't see the difference between computing the euclidian distante and
% a gaussian blur. With proper kernel size and sigma, one can achieve the same result, right?

% For FSC masks, apply a gaussian blur to the mask.
if (OPTION.fsc)
    currentMask = tofloat(currentMask);
    smoothKernel = EMC_multi_gaussianKernel(zeros(1, ndim) + rad, rad / 2, METHOD, ...
                                            {'precision', OPTION.precision});
   	for i = 1:2  % set to one for tight mask demo
        currentMask = convn(currentMask, smoothKernel,'same');
    end
    currentMask = currentMask ./ max(currentMask(:));
end

% Apply a final taper to make sure the edges go to 0.
MASK = currentMask .* taperEdge;

% Compute the center of mass of the final MASK.
% Q4B8: If fsc mask, why aren't you computing the com on the final MASK? line 458 of BH_mask3d is weird:
% binaryVol = (binaryVol - min(binaryVol(MASK > 0.01)) ).*MASK; with binaryVol being the binary mask and
% MASK being binaryVol with blurring and edgetaper.
if (OPTION.com)
    COM = EMC_centerOfMass(MASK, nan, OPTION.origin);
else
    COM = nan;
end

% Estimate the particle fraction
if (OPTION.particle_fraction)
    % The mask is not a logical mask, it has values between 0 and 1. As such,
    % estimate the signal reduction in these regions. This is >1.
    powerReduction = sum(IMAGE(:).^2.*(currentMask(:)>0))./ ...
                     sum(IMAGE(:).^2.* currentMask(:));

    maskVolume = sum(currentMask(:)>0);

    % particleVolEstimate was calculted after the first expansion and is an estimation of
    % particle volume. Now, scale it to add its hydration surface.
    particleVolEstimate = particleVolEstimate ./ localParticleScaling;

    % Finally, compare it (hydrated particle) to the mask volume (scalled down to take into account roll off)
    % to estimate the fraction of the mask taken by the particle.
    PARTICLE_FRACTION = particleVolEstimate ./ (maskVolume .* powerReduction);

    % Q4B9: is my understanding OK? Also, what 1/PARTICLE_FRACTION represents?

    if OPTION.verbose > 0
        fprintf('Estimated partVol, %d voxels\nmaskVol %d voxels\npwrReduction %2.3f\npartFract %2.3f\n',...
                 particleVolEstimate, maskVolume, powerReduction,PARTICLE_FRACTION);
    end
else
    PARTICLE_FRACTION = nan;
end

end


function [COM] = EMC_centerOfMass(IMAGE, ORIGIN)
%
% Compute the center of mass of an 2d/3d IMAGE.
%
% IMAGE (numeric):          2d/3d IMAGE.
%
% ORIGIN (int):             Origin convention; -1, 0, 1 or 2;
%                           see EMC_multi_vectorCoordinates for more details.
%
%--------
% RETURN:                   COM: center of mass [x,y] or [x,y,z]
%


[vX, vY, vZ] = EMC_multi_vectorCoordinates(size(IMAGE), 'gpu', {'origin', ORIGIN});
array = IMAGE - min(IMAGE);  % only positive values (set min value to 0)

if (flg3d)
    COM = [sum(sum(sum(array.*vX))), ...
           sum(sum(sum(array.*vY))), ...
           sum(sum(sum(array.*vZ)))] ./ sum(array(:));
else
    COM = [sum(sum(array.*vX)), sum(sum(array.*vY))] ./ sum(array(:));
end

end  % end EMC_centerOfMass


function [IMAGE, SIZE, OPTION, flg, ndim, tofloat] = checkIN(IMAGE, PIXEL, METHOD, OPTION)
% checkIN of EMC_referenceMask

if ~isnumeric(IMAGE)
    error('IMAGE should be numeric, got %s', class(IMAGE))
end
[SIZE, flg.is3d, ndim] = EMC_is3d(size(IMAGE));

validateattributes(PIXEL, {'numeric'}, {'nonnegative', 'numel', 1}, 'checkIN', 'PIXEL');

OPTION = EMC_extract_option(OPTION, ...
            {'origin', 'fsc', 'com', 'lowpass', 'threshold', 'hydration_scaling', 'precision'}, false);
% origin
if isfield(OPTION, 'origin')
  	if ~(OPTION.origin == 0 || OPTION.origin == 1 || OPTION.origin == 2)
        error("origin should be 0, 1, or 2, got %d", ORIGIN)
  	end
else
  	OPTION.origin = 1;  % default
end

% fsc
if isfield(OPTION, 'fsc')
    if ~islogical(OPTION.fsc)
        error('fsc should be a boolean, got %s', class(OPTION.fsc))
    end
else
    OPTION.fsc = false;  % default
end

% com
if isfield(OPTION, 'com')
    if ~islogical(OPTION.com)
        error('com should be a boolean, got %s', class(OPTION.com))
    end
else
    OPTION.com = false;  % default
end

% particle_fraction
if isfield(OPTION, 'particle_fraction')
    if ~islogical(OPTION.particle_fraction)
        error('particle_fraction should be a boolean, got %s', class(OPTION.particle_fraction))
    end
else
    OPTION.particle_fraction = false;  % default
end

% lowpass
if isfield(OPTION, 'lowpass')
    if isnan(OPTION.lowpass)
        flg.lowpass = false;
    elseif isfloat(OPTION.lowpass) || isinteger(OPTION.lowpass)
        if OPTION.lowpass < PIXEL * 2
            flg.lowpass = false;
        else
            flg.lowpass = true;
        end
    else
        error('lowpass should be a float|int or nan')
    end
else
    OPTION.lowpass = 14;  % default
    if OPTION.lowpass < PIXEL * 2
        flg.lowpass = false;
    else
        flg.lowpass = true;
    end
end

% threshold
if isfield(OPTION, 'threshold')
    validateattributes(OPTION.threshold, {'numeric'}, {'nonnegative', 'numel', 1}, 'checkIN', 'threshold');
else
    OPTION.threshold = 2.5;  % default
end

% hydration_scaling
if isfield(OPTION, 'hydration_scaling')
    validateattributes(OPTION.hydration_scaling, {'numeric'}, {'nonnegative', 'numel', 1}, ...
                       'checkIN', 'hydration_scaling');
else
    OPTION.hydration_scaling = (-2.8e-3) .* OPTION.lowpass.^2 + 0.14 .* OPTION.lowpass + 1.5;  % default
end

% precision is checked by EMC_resize.
[precision, onGPU] = EMC_getPrecision(IMAGE);
if strcmpi(METHOD, 'gpu')
    flg.gpu = true;
    if ~onGPU
        IMAGE = gpuArray(IMAGE);
    end
elseif strcmpi(METHOD, 'cpu')
    flg.gpu= false;
    if onGPU
        IMAGE = gather(IMAGE);
    end
else
    error("METHOD should be 'gpu' or cpu''")
end

if isfield(OPTION, precision)
    if strcmpi(OPTION.precision, 'single')
        tofloat = @(array) single(array);
    elseif strcmpi(OPTION.precision, 'double')
        tofloat = @(array) double(array);
    else
        error("precision should be 'single' or 'double'")
    end
else
    OPTION.precision = 'single';  % default
    tofloat = @(array) single(array);
end

end  % checkIN
