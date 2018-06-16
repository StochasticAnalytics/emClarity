function [] = BH_benchmark( fileNameOut, fastScratchDisk,nWorkers)
% Time a series of operations to compare different architectures under the
% following variable conditions
%
%   Size - <=12 x 128
%          <=6 x 256^3 
%          <=3 x 384^3 (0.2265 Gb/each) x 3 times mem
%
%   Num processes - MPS has been changing, particularly on volta
%
%   In vRAM, In main memory, On SSD
%     --> I'm not sure the disks available are equivalent, may need to
%     scale this based on the models
%
%   Operations - 
%
%     Simulate average by: a multiply (mask)
%                          an Interoplation
%                          stats in real space (center and norm)
%
%
%     Simulate alignment by: a multiply (mask0
%                            FFT --> bandpass, norm
%                            conj(multiply)
%                            IFFT
%                            A maxVal op 
%                            This 2x
%
%
%
nSeconds = 3;

SIZE = 256;
nWorkers = str2num(nWorkers);
if nWorkers < 0
  nWorkers = abs(nWorkers);
  flgLoadEveryOtherTrial = 1
else
  flgLoadEveryOtherTrial = 0
end
finalContainer = cell(3,1);

for memLevel = [0,1,2]

  try
    parpool(nWorkers);
  catch
    delete(gcp('nocreate'))
    parpool(nWorkers);
  end

  % For monitoring get existing pids on system, start get new (to kill at the
  % end
  [~,existingPIDS] = system('pgrep nvidia-smi');
  existingPIDS = str2num(existingPIDS)
  system(sprintf('nvidia-smi dmon -d %d -i 1 -s pctum > %s_timingLog_%d_%d.txt &',nSeconds,fileNameOut,memLevel,nWorkers));
  
  [~,updatedPIDS] = system('pgrep nvidia-smi');
  updatedPIDS = str2num(updatedPIDS)

  % Get the newPID
  if isempty(existingPIDS)
    dmonPID = updatedPIDS
  else
    dmonPID = updatedPIDS(~ismember(updatedPIDS,existingPIDS))
  end
  
  clear updatedPIDS existingPIDS

  
  nTrials = 10;
  parContainer = cell(nWorkers,1);
  parfor iProc = 1:nWorkers
    
    g = gpuDevice(1);

    % Create the volume to operate on
    if memLevel == 2
      % Create and write to disk
      volumeData = sprintf('%s/bench_%d_%d.mrc',fastScratchDisk,iProc,SIZE);
      volOp = randn(SIZE.*[1,1,1],'single');
      SAVE_IMG(MRCImage(volOp),volumeData);
      volOp = [];
    elseif memLevel == 1
      volumeData = 0;
      volMem = randn(SIZE.*[1,1,1],'single');
    elseif memLevel == 0
      volumeData = 0;
      volOp = randn(SIZE.*[1,1,1],'single','gpuArray');    
    else
      error('memLevel 0,1,2 gpu,main,disk, not %d',memLevel);
    end

    % Run the "averaging" benchmark
    aliTiming = 0;
    readTiming = 0;
    for iTrial = 1:nTrials    

      tic
      volOp = [];
      if ( volumeData )
        volOp = gpuArray(getVolume(MRCImage(volumeData)));
      elseif memLevel == 1
        volOp = gpuArray(volMem);
      else
        volOp = randn(SIZE.*[1,1,1],'single','gpuArray');
      end
      readTiming = readTiming + toc;

      tic
      bp = BH_bandpass3d(SIZE.*[1,1,1]+2,0,60,4,'GPU',1);
      volOp = volOp .* (volOp > 1);
      volOp = BH_bandLimitCenterNormalize(volOp,bp,'',[1,1,1;1,1,1],'singleTaper');
      bp = [];
      for iAng = 1:50
        % for the case where data i/o becomes an issue.
        if ( flgLoadEveryOtherTrial && (1-mod(iAng,2)))
          volOp = [];
          if ( volumeData )
            volOp = gpuArray(getVolume(MRCImage(volumeData)));
          elseif memLevel == 1
            volOp = gpuArray(volMem);
          else
            volOp = randn(SIZE.*[1,1,1],'single','gpuArray');            
          end 
        end
        
        volTmp = BH_resample3d(volOp,randn(1,3).*100,randn(1,3).*5,'Bah','GPU','inv');      
        volTmp = real(ifftn(fftn(volOp).*conj(fftn(volTmp))));
        max(volTmp(:));
        volTmp = [];
      end
      volOp = volOp(2:end-1,2:end-1,2:end-1);
      aliTiming = aliTiming + toc;
    end
    
    
    

   parContainer{iProc} = [readTiming,aliTiming]./nTrials;
  end % end of parfor loop
  
  avgTiming = [0,0];
  for iProc = 1:nWorkers
    avgTiming = avgTiming + parContainer{iProc};
  end
  % Average time per worker, then normalize to number of workers
  finalContainer{memLevel+1} = avgTiming./(nWorkers^2);
  
  [failToKill]=system(sprintf('kill %d',dmonPID));
  
%   if ( failToKill )
%     finalContainer{memLevel+1}{2} = 'did not kill nvidia-smi dmon';
%   else
%     finalContainer{memLevel+1}{2} = 'successfully killed nvidia-smi dmon';
%   end
end

save(sprintf('%s_%d.mat',fileNameOut,nWorkers),'finalContainer');

end % end of benchmark function


