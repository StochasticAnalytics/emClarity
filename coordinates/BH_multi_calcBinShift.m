function [ binSize, binShift ] = BH_multi_calcBinShift(coords, samplingRate, force_odd_dimension)
% Address fractional shifts on binning
%   Coordinates are stored relative to the lower left corner of the full
%   tilt/tomo. On binning a shift is needed to keep that origin in the same
%   place, so that any operations done on the binned data, when scaled back
%   to full sampling, leaves the coordinates unchanged. I think the easiest
%   way to do this is to shift the data on binning in 2d, and leave the
%   coordinates alone.

  % Expecting just the x,y,z for a tilt series and the binning. Also may
  % shift to have an odd dimension so that Imod origin is always the same.
  binSize = floor(coords./samplingRate);
  if (force_odd_dimension)
    binSize = binSize - (1-mod(binSize,2));
  end
  
  originFull = emc_get_origin_index(coords);
  originBin  = emc_get_origin_index(binSize);

  % 1 ++++++^+++^
  % 2 _ _ _ _ _ _
  % 3 ___ ___ ___ 
  % If the continuous specimen in on line 1 and the unbinned image is sampling that specimen as in line 2
  % Our goal is to have the feature that is on origin 1 (pixel 4 = ^) to be on origin 2 (pixel 2 = ^)

  % You know, typing this out makes me think it is unneeded.
  
  % This is the shift we need to apply to the binned image to make sure
  % that the origin is in the same place.
  binShift = -1.*(samplingRate.*originBin - originFull) ./ samplingRate;


end

