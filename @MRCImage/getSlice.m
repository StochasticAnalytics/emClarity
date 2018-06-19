%getSlice       Get a slice of data 
%
%   slice = getImage(mRCImage, idxDomain, index)
%
%   slice       The selected image slice
%
%   mRCImage    The MRCImage object.
%
%   idxDomain   The index of the domain to hold constant for extracting the
%               slice: 1=I, 2=J, 3=K
%
%   index       The value of the constant domain.
%
%   Bugs: inefficient in that it loads to much data
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

function slice = getSlice(mRCImage, idxDomain, index)

dimensions = [mRCImage.header.nX mRCImage.header.nY mRCImage.header.nZ];


if index < 1 || index > dimensions(idxDomain)
  PEETError('Image index out of range!');
end

%FIXME implement a more efficient access than getVolume
iRange = [1 dimensions(1)];
jRange = [1 dimensions(2)];
kRange = [1 dimensions(3)];

switch idxDomain
 case 1,
  iRange = [index index];
 case 2,
  jRange = [index index];
 case 3,
  kRange = [index index];
 otherwise,
  PEETError('Domain index must be 1, 2, or 3!');
end

slice = getVolume(mRCImage, iRange, jRange, kRange);
