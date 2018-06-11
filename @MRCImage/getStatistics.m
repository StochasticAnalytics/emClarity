%getStatistics  Return the selected statistic(s) of the volume
%
%   stat = getStatistics(mRCImage, statistic, domain)
%
%   stat        The requested statistic(s)
%
%   mRCImage    The mRCImage object to analyze.
%
%   statistic   The statistic to calculate:
%               'min', 'max', 'mean', 'median'
%
%   domain      OPTIONAL: The domain over which the statistic will be
%               calculated: ('z')
%               'x', 'y', 'z', 'global'
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

function stat = getStatistics(mRCImage, statistic, domain)
this=mRCImage;

%  Permute/stack the data if according to the domain selection
if nargin < 3
  domain = 'z';
end

switch lower(domain)
 case 'x'
    data = reshape(permute(this.volume, [2 3 1]),                      ...
                   this.header.nY * this.header.nZ,                    ...
                   this.header.nX);
 case 'y'
    data = reshape(permute(this.volume, [1 3 2]),                      ...
                   this.header.nX * this.header.nZ, ...
                   this.header.nY);
 case 'z'
   data = reshape(this.volume, this.header.nX * this.header.nY,        ...
                  this.header.nZ);
 case 'global'
  data = this.volume(:);
 
 otherwise
  PEETError(['Invalid domain selector: ' domain]);
end


switch lower(statistic)
 case 'min'
  stat = min(data);
 
 case 'max'
   stat = max(data);
 
 case 'mean',
   stat = mean(data);
 
 case 'std',
   stat = std(data);
 
 case 'median',
    stat = median(data);
 
 otherwise
  PEETError(['Unimplemented statistic: ' statistic]);
end
