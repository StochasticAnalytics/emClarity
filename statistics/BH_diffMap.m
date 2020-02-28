function [ diffMap, normMap ] = BH_diffMap( refMap, ...
                                            particle, particleCTF, ...
                                            flgNorm, pixelSize,...
                                            radialGrid, padVal )
%Scale a higher SNR map down do approximate the power of a noisy map
%   Maps are assumed to be the same size and also be masked in real space. If
%   the particle ctf is just 1, and flgNorm than this should behave like Niko's
%   diffmap, if not ctf is 1 and flgNorm is 0 than this is just an ordinary
%   difference.
%
%   The ctf can be a binary wedge mask or a full 3d-CTF which pre-multiplies the
%   refMap FFT so that the normalization factor has a similar sampling compared
%   to the particle.
normMap = '';
% Next check to see if a radial grid is provided and that it is the right size.
if isnumeric(radialGrid)
  if (size(radialGrid) ~= size(particleCTF))
    error('radialGrid and particleCTF are different sizes');
  end
elseif (flgNorm)
  [radialGrid,~,~,~,~,~] = BH_multi_gridCoordinates(size(refMap),'Cartesian',...
                                                    'GPU',{'none'},1,0,1);
  radialGrid = radialGrid ./ pixelSize;
end


% Assuming the maps are of the same particle, are aligned, and of the same size.
% Whether the ctf is a binary wedge, or a full 3d-CTF pre multiply the refMap so
% that the average power calculated is accurate.

% Don't take abs to re-use later.
if isreal(refMap(1))
  refMap = fftn(BH_padZeros3d(refMap,padVal(1,:),padVal(2,:),'GPU','single')).*particleCTF;
else
  refMap = refMap.*particleCTF; 
end

if isreal(particle)
  particle= fftn(BH_padZeros3d(particle,padVal(1,:),padVal(2,:),'GPU','single'));
else
  particle= BH_padZeros3d(particle,padVal(1,:),padVal(2,:),'GPU','single',0,1);
end


% Just in case they aren't centered, set mean to zero
refMap(1) = 0;
particle(1) = 0;

refMap = refMap ./ (sqrt(sum(sum(sum(abs(refMap(particleCTF > 1e-2)).^2))))./numel(refMap));
particle = particle ./ (sqrt(sum(sum(sum(abs(particle(particleCTF > 1e-2)).^2))))./numel(particle));

if (flgNorm)
  binDiv = 40;
  bin = floor(size(refMap,1)/binDiv);
  inc = 0.5 / (bin*pixelSize);
  shellsFreq = zeros(bin, 1, 'gpuArray');
  shellsFSC = zeros(bin, 1, 'gpuArray');

  for q = 1:bin

    iMask = gpuArray((q-1)*inc <= radialGrid & radialGrid < (q)*inc);
    shellsFreq(q, 1) = inc.*(q-1/2);

%    % The have the same number of pixels so no need to consider in average
%    shellsFSC(q, 1) = real(sum(double(abs(refMap(iMask)))))./...
%                      real(sum(double(abs(particle(iMask)))));

    % The have the same number of pixels so no need to consider in average
    shellsFSC(q, 1) = real(sum(double(abs(refMap(iMask)))))./...
                      real(sum(double(abs(particle(iMask)))));


  end
  clear iMask 

  % cFit to oversample the scaling factor
  fitScale = fit(gather(shellsFreq(:,1)),gather(shellsFSC(:,1)),'cubicSpline');
  normMap = refMap ./ reshape(fitScale(radialGrid),size(radialGrid));
  clear fitScale radialGrid
  diffMap = real(ifftn(normMap - particle));
else
%  refMap = refMap ./ sum(abs(refMap(:)).^2); % FIXME 
%  particle = particle ./ sum(abs(particle(:)).^2);

  normMap = '';
  diffMap = real(ifftn(refMap - particle));
end


% % figure, imshow3D(gather(real(ifftn(particle))));
% % figure, imshow3D(gather(real(ifftn(refMap))));
% % if isnumeric(normMap); figure, imshow3D(gather(real(ifftn(normMap)))); end
% % figure, imshow3D(gather(real((diffMap))));
% % 
% % 
% % error('sdf')

% [ PEAK_COORD ] =  BH_multi_xcf_Translational( particle, conj(refMap), ...
%                                                         '', [1,1,1])
% % % SAVE_IMG(MRCImage(gather(real(ifftn(refMap)))),'refMap.mrc',pixelSize);
% % % if ~isempty(normMap)
% % %   SAVE_IMG(MRCImage(gather(real(ifftn(normMap)))),'normMap.mrc',pixelSize);
% % % end
% % % SAVE_IMG(MRCImage(gather(real(ifftn(particle)))),'particle.mrc',pixelSize);
% % % SAVE_IMG(MRCImage(gather(diffMap)),'diffMap.mrc',pixelSize);
% % % error('sdf')

% This second trimming is a waste of time in pca
if any(padVal)
  diffMap = BH_padZeros3d(diffMap,-1.*padVal(1,:),-1.*padVal(2,:),'GPU','single');
  if isnumeric(normMap)
    normMap = BH_padZeros3d(normMap,-1.*padVal(1,:),-1.*padVal(2,:),'GPU','single');
  end
end

clear refMap particle


end

