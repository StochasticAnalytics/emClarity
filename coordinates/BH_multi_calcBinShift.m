function [ binSize, binShift ] = BH_multi_calcBinShift(coords, samplingRate)
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
binSize = binSize - (1-mod(binSize,2));

originFull = floor(coords ./2) + 1;
originBin  = floor(binSize./2) + 1;
% This is the shift we need to apply to the binned image to make sure
% that the origin is in the same place.
binShift = -1.*(samplingRate.*originBin - originFull) ./ samplingRate;

end

