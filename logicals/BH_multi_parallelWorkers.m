function [ nWorkers ] = BH_multi_parallelWorkers(nWorkers_wanted)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

pInfo = parcluster();
if (pInfo.NumWorkers < nWorkers_wanted)
  fprintf('\nnWorkers requested but only %d are visible, reducing\n', ...
    nWorkers_wanted,pInfo.NumWorkers);
  nWorkers = pInfo.NumWorkers;
else
  nWorkers = nWorkers_wanted;
end

