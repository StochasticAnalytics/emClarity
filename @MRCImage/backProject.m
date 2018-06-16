%backProject    Compute the filtered back projection of the MRCImage object 
%
%   backProj = backProject(mRCImage, projectionAngles, filter, interp, D)
%
%   result      Output description
%
%   parm        Input description [units: MKS]
%
%   Optional    OPTIONAL: This parameter is optional (default: value)
%
%
%   TEMPLATE Describe function, it's methods and results.
%
%   Calls: none
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
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [backProj, filtproj]= ...
    backProject(mRCImage, theta, xPoints, ySlices, zPoints, fc, ft)

theta = theta * (pi / 180);
nSamples = getNX(mRCImage);
nProjections = getNZ(mRCImage);

% if xPoints is empty cacluate the full volume given the range of z and the
% the number of samples
%
% Assumptions: the range of tilts is greater than:
%  atan(max(zPoints)/max(xPoints))

% Calculate the FFT size necessary to prevent spatial aliasing and the
% ramp filter
nFFT = 2 ^ nextpow2(2 * nSamples);
[H nRamp] = calcRamp(fc, ft, nFFT);


if isempty(xPoints)
  projDomainMax = (nSamples - nRamp - 1) / 2;
  xMax = floor(min(sqrt((projDomainMax)^2 - max(zPoints)^2), ...
             sqrt((projDomainMax)^2 - min(zPoints)^2))) - 1;
  xPoints = -xMax:xMax;
end

if isempty(ySlices)
  nSlices = getNY(mRCImage);
  ySlices = [1 nSlices];
else
  nSlices = length(ySlices);
end

% Calculate the interpolation indices for each projection angle
% using a linear interpolation
[idxLow, wLow, idxHigh, wHigh] = ...
    interpWeights(nSamples - nRamp, theta, xPoints, zPoints);

% Loop over each slice computing the filtered back projection
nX = length(xPoints);
nZ = length(zPoints);

volume = zeros(nX, nZ, nSlices);

projMin = double(getStatistics(mRCImage, 'mean'));
for iSlice = 1:nSlices
  disp(int2str(iSlice))
  % Extract the current slice
  idxSlice = ySlices(iSlice);
  proj = double(squeeze(getVolume(mRCImage, ...
                                  [1 nSamples], ...
                                  [idxSlice idxSlice], ...
                                  [1 nProjections])));
%  for iProj = 1:nProjections
%    proj(:, iProj ) = proj(:, iProj) - projMin(iProj);
%  end
%  proj = proj + min(proj(:)) + 1;
%  proj = log(proj);
  
  filtproj = filterProjection(proj, H, nFFT, nSamples, nRamp);
  
  volume(:, :, iSlice) = backProjection(filtproj, idxLow, wLow, idxHigh, wHigh);
end

% Construct the MRCImage containing the volume
backProj = volume;
return

function filtproj = filterProjection(proj, H, nFFT, nSamples, nRamp)
PROJ = fft(proj, nFFT);
nProj = size(PROJ,2);
filtproj = zeros(nSamples - nRamp, nProj);
for i = 1:nProj
  %% TODO get the length of the ramp filter to remove the transients
  temp = real(ifft(PROJ(:, i) .* H));
  filtproj(:, i) = temp(nRamp+1:nSamples);
end

return

function rec = backProjection(filtproj, idxLow, wLow, idxHigh, wHigh)
nProjections = size(filtproj, 2);
rec = zeros(size(idxLow, 1), size(idxLow, 2));

for iTheta = 1:nProjections
  proj = filtproj(:, iTheta);
  rec = rec ...
      + wLow(:, :, iTheta) .* proj(idxLow(:, :, iTheta)) ...
      + wHigh(:, :, iTheta) .* proj(idxHigh(:, :, iTheta));
end
return

%%
%%  Calculate the linear back projection interpolate weights and indices
%%  TODO: probably need to zero out the weights and set the indices to one
%%  for any X-Z coordinate out of the domain of the projection
%%
function [idxLow, wLow, idxHigh, wHigh] = ...
    interpWeights(nSamples, theta, xPoints, zPoints)

nXPoints = length(xPoints);
nZPoints = length(zPoints);
[xPoints zPoints] = ndgrid(xPoints, zPoints);

% t contains the coordiate along the projection for each point in X and Z
% and for each theta
% TODO: can (should) this be vectorized
nTheta = length(theta);
t = zeros(nXPoints, nZPoints, nTheta);
for iTheta = 1:nTheta
  t(:, :, iTheta) = xPoints .* cos(theta(iTheta)) ...
      + zPoints .* sin(theta(iTheta)) + nSamples / 2;
end
% Generate the linear interpolation weights and indices for each
% The mapping from position to index requires the addition of 1 since
% MATLAB uses base 1 indexing
idxLow = floor(t + 1);
idxHigh = idxLow + 1;
wHigh = (t - idxLow);
wLow = 1 - wHigh;


%%
%%  Calculate the ramp filter
%%
%%  fc   The cutoff frequency of the filter (cycles/sample)
%%
%%  ft   The transition width between the pass and stop bands (cycles/sample)
%%
function [H, n] = calcRamp(fc, ft, nPoints)
[n wn beta ftype] = kaiserord([2*fc 2*(fc + ft)], [1 0], [0.05 0.05]);
fprintf('Kaiser filter order: %d\n', n);
hlp = fir1(n, wn, ftype, kaiser(n+1, beta));
HLP = fft(hlp', nPoints);
H = HLP .* ([0:nPoints/2-1 nPoints/2-1:-1:0] ./ (nPoints / 2))';

