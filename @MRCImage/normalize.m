%normalize      Normalize each projection to unity integrated density
%
%   mRCImage = normalize(mRCImage)
%
%   mRCImage    The MRCImage object.
%
%
%   Bugs: none known
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

function mRCImage = normalize(mRCImage)

%  Load the volume if it is not already loaded
if ~ mRCImage.flgVolume
  mRCImage = loadVolume(mRCImage);
end


% Normalize each projection to unity integral

for iProj = 1:mRCImage.header.nZ
  totalSum = sum(sum(mRCImage.volume(:, :, iProj)));
  mRCImage.volume(:, :, iProj) = mRCImage.volume(:, :, iProj) ./ totalSum;
end

