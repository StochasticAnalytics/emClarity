function [ ] = BH_runCtfFind2(PARAMETER_FILE, stackName, tltName, ctfParams)
%Fit the ctf to a background subtracted PS using ctffind4
%   CTF params
%     PixelSize (Ang) 
%     KeV 
%     CS (mm)
%     Amplitude Contrast

system('mkdir -p fixedStacks/ctf/forCtfFind');
pBH = BH_parseParameterFile(PARAMETER_FILE);

n_threads = pBH.nCpuCores;

[~,tltBase,~] = fileparts(tltName);

% rng('shuffle');
% randPrfx = sprintf('%s_%d',tltBase,randi(1e6,[1,1]));

randPrfx = sprintf('%s',tltBase);

ctfFindPath = getenv('EMC_CTFFIND');

tmpCache= pBH.('fastScratchDisk');

if strcmpi(tmpCache, 'ram') 
  if isempty(getenv('EMC_CACHE_MEM'))
    fprintf('Did not find a variable for EMC_CACHE_MEM\nSkipping ram\n');
    tmpCache= '';
  else
    tmpCache=getenv('MCR_CACHE_ROOT');
    fprintf('Using the tmp EMC cache in ram at %s\n',tmpCache);
  end
end

try
  size_of_amplitude_spectrum = pBH.('size_of_amplitude_spectrum');
catch
  size_of_amplitude_spectrum = 512;
end

tiltAngles = load(tltName);
% split the stack up
fprintf('%s\n',ctfFindPath)

% mirror directory in memory for larger temp files
system(sprintf('mkdir -p %s/forCtfFind',tmpCache));

stack_header = getHeader(MRCImage(stackName,0));

d1 = stack_header.nX;
d2 = stack_header.nY;
d3 = stack_header.nZ;

% Temp override for dev
d3 = 1;
slice_name = cell(d3,1);

for iPrj = 1:d3
  slice_name{iPrj} = sprintf('%s/forCtfFind/%s_%d.mrc',tmpCache,randPrfx,iPrj);
%   system(sprintf('newstack -fromone -secs %d %s %s > /dev/null',iPrj, stackName, slice_name{iPrj}));
  system(sprintf('newstack -fromone -secs %d %s %s',iPrj, stackName, slice_name{iPrj}));
end


system(sprintf('mv fixedStacks/ctf/%s fixedStacks/ctf/%s_orig',tltName,tltName));

tmpTLT = load(sprintf('%s_orig',tltName));
meanDefocus = mean(tmpTLT(:,15))*-1.0*10^10;
fprintf('Searching around an estimated mean defocus of %3.6f Angstrom\n');

% write the run script, this should link to a distributed version with
% special name, but for testing use the beta.

%         **   Welcome to Ctffind   **
% 
%             Version : 4.1.14
%            Compiled : Dec  9 2021
%     Library Version : 2.0.0-alpha-22-9131cb2-dirty
%         From Branch : add_gpu_methods_for_unblur
%                Mode : Interactive
% 
% Input image file name [rot.mrc]                    : 
% Output diagnostic image file name
% [diagnostic_output.mrc]                            : 
% Pixel size [3.24]                                  : 
% Acceleration voltage [200]                         : 
% Spherical aberration [2.70]                        : 
% Amplitude contrast [0.1]                           : 
% Size of amplitude spectrum to compute [512]        : 
% Minimum resolution [30.0]                          : 
% Maximum resolution [12]                            : 
% Minimum defocus [15000]                            : 
% Maximum defocus [60000]                            : 
% Defocus search step [100.0]                        : 
% Do you know what astigmatism is present? [No]      : 
% Slower, more exhaustive search? [No]               : 
% Use a restraint on astigmatism? [yes]              : 
% Expected (tolerated) astigmatism [200.0]           : 
% Find additional phase shift? [No]                  : 
% Determine sample tilt? [yes]                       : 
% Do you want to set expert options? [yes]           : 
% Resample micrograph if pixel size too small? [Yes] : 
% Do you already know the defocus? [No]              : 
% Desired number of parallel threads [7]             : 

scriptName = sprintf('.%s.sh',randPrfx);
fID = fopen(scriptName,'w');
fprintf(fID,'#!/bin/bash\n\n');
for iPrj = 1:d3 % I want to fit to lower resolution at higher tilts
  tltIDX = find(tiltAngles(:,1) == iPrj);

  % put in a line to limit number of cores, or use the threaded version
  fprintf(fID,'\n%s << eof &',ctfFindPath);
  fprintf(fID,'\n%s/forCtfFind/%s_%d.mrc\n',tmpCache, randPrfx,iPrj);
  fprintf(fID,'fixedStacks/ctf/forCtfFind/%s_diagnostic_%d.mrc\n',randPrfx,iPrj);
  fprintf(fID,'%f\n%f\n%f\n%f\n%d\n%f\n%f\n%d\n%d\n%d\n',ctfParams(1:4), ...
                                                         size_of_amplitude_spectrum,30,2.5*ctfParams(1)./sqrt(cosd(tiltAngles(tltIDX,4))),...
                                                         0.75*meanDefocus,...
                                                         1.25*meanDefocus,...
                                                         50.0);
  fprintf(fID,'no\nno\nyes\n200.0\nno\nyes\nyes\nyes\nno\n%d\neof\n\n', n_threads);
end
fprintf(fID,'wait\n');
fclose(fID);

system(sprintf('chmod a=wrx %s',scriptName));

[runFail] = system(sprintf('./%s',scriptName));

if (runFail)
  system(sprintf('cp ./%s tmpFail',scriptName));
  system(sprintf('mv tmpFail ./%s',scriptName));
  [runFail] = system(sprintf('./%s',scriptName));
  if (runFail)
    error('Tried to run %s twice and failed\n',scriptName);
  end
end

error('run')
% will this wait for return?

baseName = sprintf('fixedStacks/ctf/forCtfFind/%s_diagnostic_',randPrfx);
tmpName  = sprintf('fixedStacks/ctf/forCtfFind/%s_tmp',randPrfx);

system(sprintf('newstack %s?.mrc %s??.mrc %s_full.st',baseName,baseName,baseName));
system(sprintf('rm %s?.mrc %s??.mrc',baseName,baseName));
system(sprintf('rm -f %s',tmpName));

for iPrj = 1:d3
  
	system(sprintf('tail -n -1 %s%d.txt | awk  ''{print (($2-$3)/2)*10^-10,3.141592/180*$4,-1*(($2+$3)/2)*10^-10 }'' >> %s', baseName,iPrj,tmpName));                               
                                         
end

% TODO ground truth to confirm orientation of astigmatism



system(sprintf('awk ''FNR==NR{a[FNR]=$1;b[FNR]=$2;c[FNR]=$3 ;next}{ print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,a[$1],b[$1],$14,c[$1],$16,$17,$18,$19,$20,$21,$22,$23}'' %s fixedStacks/ctf/%s_orig > fixedStacks/ctf/%s',tmpName,tltName,tltName));
