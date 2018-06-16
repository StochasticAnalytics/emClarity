function [ maskType, maskSize, maskRadius, maskCenter ] = BH_multi_maskCheck( pBH, PREFIX, pixelSize )
%Check that any masked sub-regions are in bounds with proper apodization given
%the offset.
%   

maskType   = pBH.(sprintf('%s_mType',PREFIX));
% Although the minimum roll-off is set for 6 pixels on each edge, use 14 rather 
% than 12 to make sure a shift in the origin doesn't preculde the volume from
% analysis.
maskSize = 2.*floor(pBH.('particleRadius') ./ pixelSize)+14;
% center can be floats.
maskRadius = round(pBH.(sprintf('%s_mRadius',PREFIX)) ./ pixelSize);

maskSize = max([maskSize;2.*maskRadius+14],[],1);
maskCenter = pBH.(sprintf('%s_mCenter',PREFIX)) ./ pixelSize;

if any(maskSize < maskRadius+abs(maskCenter(1:3))+6)
  error('The requested subRegion radius and origin offset are larger than allowed by the maskRadius.\n')
end



end

