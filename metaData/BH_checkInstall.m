function [  ] = BH_checkInstall( runPath )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


% system(sprintf('%s',getenv('BH_CHECKINSTALL')));
system(sprintf('%s',runPath));

fOUT = fopen('emClarity_checkInstall.txt','a');


[status, returnVal] = system('which chimera');

if (status)
  fprintf(fOUT,'\nNo chimera installation detected on system path.\n\n');
else
  [~,version] = system('chimera --version');
  fprintf(fOUT,'\nChimera found on path at %s\n',returnVal);
  fprintf(fOUT,'%s\n\n',version);
end

[status, returnVal] = system('which imod');


if (status)
  fprintf(fOUT,'\nNo imod installation detected on system path.\n\n');
else
  [~,version] = system('imod -h | head -n 1');
  fprintf(fOUT,'\nimod found on path at %s\n',returnVal);
  fprintf(fOUT,'%s\n\n',version);
end

nGPUs = gpuDeviceCount;

fprintf(fOUT,'Found %d gpus on the system\n', nGPUs);

for iGPU = 1:nGPUs
  fprintf(fOUT,'\n\n##########\ngpu = %d\n##########\n\n',iGPU);
  g = gpuDevice(iGPU);
  fprintf(fOUT,['Name:\t%s\nComputeCapability:\t%s\nDriverVersion:\t%f\n',...
    'ToolkitVersion:\t%f\nTotalMemory:\t%e\nAvailableMemory:\t%e\n',...
    'MultiProcessorCount:\t%d\nClockRate:\t%f\nComputeMode:\t%s\n\n'],...
    g.Name,g.ComputeCapability,g.DriverVersion,g.ToolkitVersion,...
    g.TotalMemory,g.AvailableMemory,g.MultiprocessorCount,...
    g.ClockRateKHz,g.ComputeMode);
end

fclose(fOUT);

end

