function [ ] = BH_runCtfFind(stackName, tltName, ctfParams, tiltAngles)
%Fit the ctf to a background subtracted PS using ctffind4
%   CTF params
%     PixelSize (Ang)
%     KeV
%     CS (mm)
%     Amplitude Contrast

system('mkdir -p fixedStacks/ctf/forCtfFind');

rng('shuffle');
randPrfx = sprintf('%s_%d',tltName,randi(1e6,[1,1]));

ctfFindPath = getenv('EMC_CTFFIND');

fprintf('%s\n',ctfFindPath);% split the stack up
fullStack = getVolume(MRCImage(stackName));
[d1,d2,d3] = size(fullStack); % FIXME d1 assumed to equal d2 Add check in saving

for iPrj = 1:d3
  SAVE_IMG(MRCImage(fullStack(:,:,iPrj)),sprintf('fixedStacks/ctf/forCtfFind/%s_%d.mrc',randPrfx,iPrj));
end

% % Check to make sure this hasn't alread been done
% if ~exist(sprintf('fixedStacks/ctf/%s_orig',tltName), 'file')
system(sprintf('mv fixedStacks/ctf/%s fixedStacks/ctf/%s_orig',tltName,tltName));
% end

tmpTLT = load(sprintf('fixedStacks/ctf/%s_orig',tltName));
meanDefocus = mean(abs(tmpTLT(:,15)))*10^10;
fprintf('Searching around an estimated mean defocus of %3.6f Angstrom\n');

% write the run script, this should link to a distributed version with
% special name, but for testing use the beta.
scriptName = sprintf('.%s.sh',randPrfx);
fID = fopen(scriptName,'w');

fprintf(fID,'#!/bin/bash\n\n');
for iPrj = 1:d3 % I want to fit to lower resolution at higher tilts
  tltIDX = find(tiltAngles(:,1) == iPrj);
  
  % put in a line to limit number of cores, or use the threaded version
  fprintf(fID,'\n%s --amplitude-spectrum-input << eof &',ctfFindPath);
  fprintf(fID,'\nfixedStacks/ctf/forCtfFind/%s_%d.mrc\n',randPrfx,iPrj);
  fprintf(fID,'fixedStacks/ctf/forCtfFind/%s_diagnostic_%d.mrc\n',randPrfx,iPrj);
  fprintf(fID,'%f\n%f\n%f\n%f\n%d\n%f\n%f\n%d\n%d\n%d\n',ctfParams(1:4), ...
    d1,30,3*ctfParams(1)./cosd(tiltAngles(tltIDX,4)),...
    0.75*meanDefocus,...
    1.25*meanDefocus,...
    25.0);
  fprintf(fID,'no\nno\nyes\n500.0\nno\nno\nno\neof\n\n');
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

% will this wait for return?

baseName = sprintf('fixedStacks/ctf/forCtfFind/%s_diagnostic_',randPrfx);
tmpName  = sprintf('fixedStacks/ctf/forCtfFind/%s_tmp',randPrfx);

system(sprintf('newstack %s?.mrc %s??.mrc %s_full.st',baseName,baseName,baseName));
system(sprintf('rm %s?.mrc %s??.mrc',baseName,baseName));
system(sprintf('rm -f %s',tmpName));

for iPrj = 1:d3
  
  system(sprintf('tail -n -1 %s%d.txt | awk  ''{print (($2-$3)/2)*10^-10,3.1415926535/180.0*$4,-1*(($2+$3)/2)*10^-10 }'' >> %s', baseName,iPrj,tmpName));
  
end

% TODO ground truth to confirm orientation of astigmatism



system(sprintf('awk ''FNR==NR{a[FNR]=$1;b[FNR]=$2;c[FNR]=$3 ;next}{ print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,a[$1],b[$1],$14,c[$1],$16,$17,$18,$19,$20,$21,$22,$23}'' %s fixedStacks/ctf/%s_orig > fixedStacks/ctf/%s',tmpName,tltName,tltName));
