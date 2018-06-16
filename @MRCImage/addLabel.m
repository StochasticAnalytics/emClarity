%addLabel       Add a label (max of 10 labels allowed) to an MRCImage.
%               A label is a string of no more than 80 characters.
%
%   mRCImage = addLabel(mRCImage, label)
%
%   mRCImage    The MRCImage object.
%
%   label       The label to be added.
%
%   Add a label to the header section of an MRCImage object
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

function mRCImage = addLabel(mRCImage, label)
if mRCImage.header.nLabels < 10
  mRCImage.header.nLabels =  mRCImage.header.nLabels + 1;
  fullLabel = [label blanks(80 - length(label))];
  mRCImage.header.labels(mRCImage.header.nLabels, :) = fullLabel;

end

