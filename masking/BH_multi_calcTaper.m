function [ taper ] = BH_multi_calcTaper(pixelFallOff)
% Calculate a cosine edge
% This doesn't include the 1 that precedes it or the 0 that would follow, 
% just the transition region.

taper = 0.5+0.5.*cos((((1:pixelFallOff)).*pi)./(length((1:pixelFallOff+1))));



end
