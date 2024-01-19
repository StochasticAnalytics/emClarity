function [nParProcesses, iterList] = BH_multi_parallelJobs(nTomograms, ...
  nGPUs, ...
  calcSize,flgAvg)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Rough Scaling of processes by mem available. Assuming all gpus are
% equivalent, which is often true, but this should be improved. The current
% approach below seems to be okay for 11-12 Gb cards ... try just scaling
% proportionally for now to get to 8 or 16
gDev = parallel.gpu.GPUDevice.getDevice(1);
pInfo = parcluster();
totMem = gDev.TotalMemory;
scaleMem = 1;
if totMem > 7.9e9 && totMem < 10.8e9
  % Currently not working on 8Gb cards, but this could probably happen
  scaleMem = 0.62;
elseif totMem > 10.9e9 && totMem < 12.2e9
  % 1080Ti or Tesla or Titan Xp
  scaleMem = 1.0;
elseif totMem > 12.8e9 && totMem < 17.2e9
  % TitanV V/p100 % seems to crash when it shouldn't, figure out later.
  scaleMem = 1.3;
else
  scaleMem = 2.2;
end
fprintf('found totalMem on GPU1 %3.3e, nWorkers %d, so scaling nProcs by %2.2f\n',totMem,pInfo.NumWorkers,scaleMem);
% need to add restrainst on parpool size, and also try to balance based on number of tomograms.
% even better o explicitly depend on memory available ( this will allow smaller cards to be used )

% nParProcesses =  nGPUs .*  (floor(scaleMem.*384./calcSize).^2+ceil(scaleMem.*288./calcSize));
nParProcesses =  ceil(scaleMem .* nGPUs .*  (floor(384./calcSize).^2+ceil(256./calcSize)))
fprintf('scaleMem %f nGpus %f calcSize %f nParProc %f\n',scaleMem,nGPUs,calcSize,nParProcesses);
%nParProcesses = min(2*nTomograms, nGPUs .*  (floor(scaleMem.*384./calcSize).^2+ceil(scaleMem.*288./calcSize)));
nParProcesses = min(nParProcesses, pInfo.NumWorkers);
fprintf('nParProcesses in parallelJobs %d\n',nParProcesses);
if ( flgAvg )
  if flgAvg == -1
    % Linear interpolation
    nParProcesses = min(nParProcesses,3*nGPUs);
  else
    % Fourier interp, requires more mem. Unoptimized
    nParProcesses = min(flgAvg,nParProcesses);
  end
end

fprintf('nParProcesses in parallelJobs1 %d\n',nParProcesses);



maxLen = length(1:nParProcesses:nTomograms);
newLen = maxLen;
while newLen == maxLen
  newLen = length(1:nParProcesses-2:nTomograms);
  if newLen == maxLen
    nParProcesses = nParProcesses - 1;
  else
    break
  end
end


nParProcesses = nParProcesses - mod(nParProcesses - nGPUs,2);
if (nParProcesses < nGPUs*2+1)
  nParProcesses = 2*nGPUs+1;
end
nParProcesses = min(nParProcesses,nTomograms);
if (flgAvg > 0)
  nParProcesses = min(nParProcesses,flgAvg);
end

nIters = min(nParProcesses,nTomograms);

% fprintf('Using %d workers in %d batches\n',nParProcesses,(1+nTomograms)./nParProcesses);
% Divide the tomograms up over each gpu
iterList = cell(nIters,1);
for iParProc = 1:nIters
  iterList{iParProc} = iParProc:nIters:nTomograms;
end

end

