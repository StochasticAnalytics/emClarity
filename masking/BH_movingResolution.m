function [ rmsMasks ] = BH_movingResolution(INPUT_VOL, nREGIONS, RADIUS)
%Generate masks based on spatially variable variance for reference
%generation.
%   Detailed explanation goes here

% TODO experiment with kmeans options
kDistMeasure='cityblock';
kReplicates=36;
kMaxIterations=100000;

% TODO see if this makes sense
minClassMembership = 0.02 * numel(INPUT_VOL); % fraction of total elememts


% TODO should check on GPU

% TODO fix moving Avg RMS to use gaussian not square.
localVar = gather(BH_movingRMS(INPUT_VOL - BH_movingAverage( ...
                                      INPUT_VOL,RADIUS.*[1,1,1]), RADIUS.*[1,1,1]));


% TODO should check on existance of parpool and use the correct nWorkers
[class, ~, ~] = kmeans(gather( localVar(:)), nREGIONS, ...
                              'replicates', kReplicates, ...
                              'Distance', kDistMeasure, ...
                              'MaxIter', kMaxIterations, ... % Default was 100
                              'Options', statset('UseParallel', 1) );
                            
% Check that there are no "trivial" classes
classKeep = [];
for iClass = 1:nREGIONS
  if sum(class == iClass) >= minClassMembership
    classKeep = [classKeep iClass];
  end
end

nToKeep = length(classKeep);
fprintf('\nKeeping %d of %d requested local regions\n', nToKeep,nREGIONS);
clear localVar

rmsMasks = cell(nToKeep,1);
% TODO is this stdDev reasonable?
smoothingKernel = BH_multi_gaussian3d(-1.*size(INPUT_VOL),4);
for iClass = 1:nToKeep
  rmsMasks{iClass} = zeros(size(INPUT_VOL),'single','gpuArray');
  rmsMasks{iClass}(class == classKeep(iClass)) = 1;
%   rmsMasks{iClass} = gather(real(ifftn(fftn(rmsMasks{iClass}).*smoothingKernel)));
%   rmsMasks{iClass} = rmsMasks{iClass} ./ max(rmsMasks{iClass}(:));
end



end

