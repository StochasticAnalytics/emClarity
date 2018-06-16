function [ mask ] = BH_multi_randomizeTaper( mask )
%Randomize the taper region in a mask to reduce correlation
%   Detailed explanation goes here

taperRegion = (mask < 1 & mask > 0);


% Random value between [1,2]
randVector = (rand(size(taperRegion)) + 1).^0.25;

mask(taperRegion) = mask(taperRegion) .^ randVector(taperRegion);

end

