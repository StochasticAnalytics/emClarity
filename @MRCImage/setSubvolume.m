% setSubvolume      Set a subvolume of an MRCImage
%
%   mRCImage = setSubolume(mRCImage, vol, center, suppressStatistics)
%
%   mRCImage    The opened MRCImage object.
%
%   vol         The new subvolume
%
%   vol         The volume to be inserted.
%
%   center      OPTIONAL: The indices at which to center the inserted 
%               subvolume. Default positions subvolume at origin of the
%               original volume.
%
%   suppressStatistics OPTIONAL: if true, suppress computation of image
%                      statistics for speed. Default = false. NOTE: valid 
%                      statistics *MUST* be computed before using the volume.
%               
%   setVolume sets a subvolume of the MRCImage object. The subvolume must
%   lie entirely within the existing volume and be of the same type.
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

function mRCImage = setSubvolume(mRCImage, vol, center, suppressStatistics)

% Replace a subvolume (up to the entire volume) of an existing volume.
% Expansion or change of mode is not allowed.
szVol = size(vol);
szMRC = getDimensions(mRCImage);
% Handle setting a 2D slice in a 3D volume
if length(szVol) < length(szMRC)
  n = length(szMRC) - length(szVol);
  szVol = [szVol ones(1, n)];
end

if nargin < 4
  suppressStatistics = false;
end

if nargin < 3 
  idxMin = [1 1 1];
  idxMax = szVol;
else
  idxMax = center + ceil(szVol ./ 2) - 1;
  idxMin = center - floor(szVol ./ 2);
end

if any(idxMin < 1)
  fprintf('Minimum index ');
  fprintf('%f', idxMin);
  fprintf('\n');
  PEETError('Minimum index out of range!');
end

if any(idxMax > szMRC)
  fprintf('Maximum index ');
  fprintf('%f', idxMin);
  fprintf('\n');
  PEETError('Maximum index out of range!');
end

% If the volume is already loaded return the selected indices
if mRCImage.flgVolume
  mRCImage.volume(idxMin(1):idxMax(1), ...
    idxMin(2):idxMax(2), idxMin(3):idxMax(3)) = vol;
else
  PEETError(['MRCImage.setVolume is not yet implemented for unloaded ' ...
    'volumes!']);
end

% Update the header to reflect the new min, max, mean and rms voxel values
if ~suppressStatistics
  mRCImage = setStatisticsFromVolume(mRCImage);
end

