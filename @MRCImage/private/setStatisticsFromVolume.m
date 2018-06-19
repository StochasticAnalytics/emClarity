% mRCImage = setStatisticsFromVolume(mRCImage)
%
% This file is part of PEET (Particle Estimation for Electron Tomography).
% Copyright 2000-2014 The Regents of the University of Colorado & BL3DEMC:
%           The Boulder Laboratory For 3D Electron Microscopy of Cells.
% See PEETCopyright.txt for more details.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  $Author: John Heumann $
%
%  $Date: 2014/01/13 20:00:38 $
%
%  $Revision: 6b413b88334c $
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
% I suppose it wasn't imagined this would be used for very large volumes like 
% tomograms, where setting as a double the full image would be very costly, p
% especially for tomos > 50GB -- set to loop over planes.

%% even looping is very slow, and at least for now, I certainly don't use the 
%% values in the header for anything, much less anything that need 2x precision.


function mRCImage = setStatisticsFromVolume(mRCImage)


if isreal(mRCImage.volume)

%  if ( numel(mRCImage.volume) < 2.25e9 )
%    mRCImage.header.minDensity = min(mRCImage.volume(:));
 %   mRCImage.header.maxDensity = max(mRCImage.volume(:));
    
%   mRCImage.header.meanDensity = mean(mRCImage.volume(:));
%   mRCImage.header.densityRMS = std(mRCImage.volume(:));
%  else
    mRCImage.header.minDensity = 0;%min(mRCImage.volume(:));
    mRCImage.header.maxDensity = 0;%max(mRCImage.volume(:));
    
    mRCImage.header.meanDensity = 0;%mean(mRCImage.volume(:));
    mRCImage.header.densityRMS = 0;%std(mRCImage.volume(:));
%  end
else % Use magnitudes for complex values
  error('Not set to handle saving complex images in this context!\n');
end

%function mRCImage = setStatisticsFromVolume(mRCImage)
%if isreal(mRCImage.volume)
%  mRCImage.header.minDensity = double(min(mRCImage.volume(:)));
%  mRCImage.header.maxDensity = double(max(mRCImage.volume(:)));
%  temp = double(mean(mRCImage.volume(:)));
%  mRCImage.header.meanDensity = temp;
%  mRCImage.header.densityRMS =                                         ...
%    sqrt(mean((double(mRCImage.volume(:)) - temp).^2));
%else % Use magnitudes for complex values
%  temp = (abs(double(mRCImage.volume)));
%  mRCImage.header.minDensity = min(temp(:));
%  mRCImage.header.maxDensity = max(temp(:));
%  mRCImage.header.meanDensity = mean(temp(:));
%  mRCImage.header.densityRMS =                                         ...
%    sqrt(mean((mRCImage.volume(:) - mRCImage.header.meanDensity).^2));
%end
