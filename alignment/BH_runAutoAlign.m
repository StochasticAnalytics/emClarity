function [ ] = BH_runAutoAlign(runPath,findBeadsPath,stackIN,tiltAngles,pixelSize,imgRotation,varargin)
%Run auto tiltseries alignment. For now, basically a script
%  sadf

  RESOLUTION_CUTOFF=18;
  LOW_RES_CUTOFF=800;
  PATCH_SIZE=3600 ;
  MAX_SAMPLING_RATE=3.0;
  MIN_SAMPLING_RATE=10.0;
  ITERATIONS_PER_BIN=3;
  CLEAN_UP_RESULTS=false;
  MAG_OPTION=5;
  tiltAngleOffset=0.0;
  
  imgRotation = str2num(imgRotation);
  pixelSize = str2num(pixelSize);
  
if nargin > 5
  TILT_OPTION = varargin{1}
else
  TILT_OPTION = '0'
end

inputStack = single(getVolume(MRCImage(stackIN)));


[~,baseName,ext] = fileparts(stackIN);
fixedName = sprintf('fixedStacks/%s.fixed.rot',baseName);

[~,tiltName,tiltExt] = fileparts(tiltAngles);
startDir = pwd;
wrkDir = sprintf('emC_autoAlign_%s', baseName);
system(sprintf('mkdir -p %s',wrkDir));
system(sprintf('mkdir -p fixedStacks'));

cd('fixedStacks');
system(sprintf('ln -s ../%s %s.fixed',stackIN,baseName));
cd('../');



if ~strcmp(tiltExt,'.rawtlt')
  system(sprintf('ln -s %s %s.rawtlt',tiltAngles,tiltName));
end

binHigh=ceil(MIN_SAMPLING_RATE ./ pixelSize);
if MAX_SAMPLING_RATE > 4
  binLow = MAX_SAMPLING_RATE;
else
  binLow = ceil(MAX_SAMPLING_RATE ./ pixelSize);
end
binInc = -1*ceil((binHigh- binLow)./3);

rotMat = [cosd(imgRotation),-1*sind(imgRotation),sind(imgRotation),cosd(imgRotation)];

[nX,nY,nZ] = size(inputStack);
fprintf('Preprocessing tilt-series\n');

gradientAliasFilter = BH_bandpass3d([nX,nY,1],1e-6,LOW_RES_CUTOFF,RESOLUTION_CUTOFF,'GPU',pixelSize);
%                      gradientAliasFilter = {BH_bandpass3d(1.*[nX,nY,1],0,0,0,'GPU','nyquistHigh'),...
%                        BH_bandpass3d([nX,nY,1],1e-6,LOW_RES_CUTOFF,RESOLUTION_CUTOFF,'GPU',pixelSize)};
medianFilter = 3;
for iPrj = 1:nZ
%   tmpPrj = BH_preProcessStack(gpuArray(inputStack(:,:,iPrj)),gradientAliasFilter,medianFilter);
  tmpPrj = real(ifftn(fftn(gpuArray(inputStack(:,:,iPrj))).*gradientAliasFilter));
  tmpPrj = BH_resample2d(tmpPrj,rotMat,[0,0],'Bah','GPU','forward',1,size(tmpPrj));
  inputStack(:,:,iPrj) = gather(tmpPrj);
end
fprintf('finished reprocessing tilt-series\n');
SAVE_IMG(MRCImage(gather(inputStack)),fixedName,pixelSize);

clear tmpPrj inputStack


cd(wrkDir)

rotFile = fopen('preRotXf.xf','w');
fprintf(rotFile,'%6.3f %6.3f %6.3f %6.3f 0.0 0.0\n',rotMat);
fclose(rotFile);



                                     
system('pwd')
fprintf('Running %s\n',runPath);
system(sprintf('%s %s %f %f %d %d %d %d %d %d %s %d > ./.autoAliLog_%s.txt',...
                                                   runPath, ...
                                                   baseName, ...
                                                   pixelSize, ...
                                                   imgRotation, ...
                                                   binHigh, ...
                                                   binLow, ...
                                                   binInc,...
                                                   nX,nY,nZ,ext,...
                                                   TILT_OPTION,baseName));
    
cd(sprintf('%s',startDir));

system(sprintf('rm %s',fixedName));

if strcmpi(TILT_OPTION,'0')
  % If not fitting tilt angles we need a copy of them with .tlt
   system(sprintf('cp %s fixedStacks/%s.tlt',tiltAngles,tiltName));
end

cd fixedStacks
system(sprintf('%s %s %d %d %d %d', findBeadsPath, baseName,...
                                       nX,nY,3000,...
                                       ceil(1.05*100/pixelSize)));
cd ..
end

