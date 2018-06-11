%getMX          Return the number of x elements
%
%   mX = getMX(mRCImage)
%
%   mX          The number of X elements.
%
%   mRCImage    The MRCImage object.
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

function mX = getMX(mRCImage)
mX = mRCImage.header.mX;
