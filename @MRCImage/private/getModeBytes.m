%getModeBytes   Return the number of bytes per element for the specified mode
%
%   nBytes = getModeBytes(mRCImage)
%
%   nBytes      Get the number of bytes per data element in the MRCImage file.
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

function nBytes = getModeBytes(mRCImage)

switch  mRCImage.header.mode
 case 0
  nBytes = 1; % byte
 case 1
  nBytes = 2; % int16
 case 2
  nBytes = 4; % single
 case 3
  nBytes = 4; % complex int16
 case 4
  nBytes = 8; % complex single
 case 6
   nBytes = 2; % uint16
 case 12
   nBytes = 2; % half (16-bit float)
 otherwise
  nBytes = -1;
end
