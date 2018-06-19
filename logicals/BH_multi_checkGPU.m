function [ useGPU ] = BH_multi_checkGPU( gpuIDX )
%Check to see if the gpu is available and assign it to be used.
%   No explicit mem checks, aside from sending to the gpu with the most
%   available if the requested gpu doesn't exist.

% Check if a specific gpu is requested, and that it is available.

nGpus = gpuDeviceCount;
% if the number is okay just use it.
if any(ismember(1:nGpus, gpuIDX))
  useGPU = gpuIDX;
  fprintf('Using gpuIDX: %d\n', gpuIDX);
else
  % otherwise send to the gpu with the most available memory, which isn't
  % necessarily going to be enough. Later try to add explicit checks.
  mostMem = zeros(nGpus,2);
  for iGpu = 1:nGpus
    gpuDev = gpuDevice(iGpu);
    mostMem(iGpu,:) = [iGpu, gpuDev.AvailableMemory/gpuDev.TotalMemory];
  end
  [val,ind] = max(mostMem(:,2));
  fprintf('There are %d gpus, you selected %d which is not valid.\n', nGpus,gpuIDX);
  fprintf('Sending to gpu with the most available memory %d with %0.2f%%\n', ind, 100.*val);
  useGPU = ind;
end


end

