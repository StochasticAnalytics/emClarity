function [ ] = BH_runAutoAlign(PARAMETER_FILE, runPath,findBeadsPath,stackIN,tiltAngles,imgRotation)
%Run auto tiltseries alignment. For now, basically a script
%  sadf

% TODO add options for experimenting.
pBH = BH_parseParameterFile(PARAMETER_FILE);

pixelSize = pBH.('PIXEL_SIZE').*10^10;
imgRotation = str2double(imgRotation);

try 
  RESOLUTION_CUTOFF = pBH.('autoAli_max_resolution');
catch
  RESOLUTION_CUTOFF=18;
end
try
  MAX_SAMPLING_RATE = pBH.('autoAli_max_sampling_rate');
catch
  MAX_SAMPLING_RATE = 4.0;
end
try
  PATCH_SIZE_FACTOR = pBH.('autoAli_patch_size_factor');
catch
  PATCH_SIZE_FACTOR = 4;
end

try
  REFINE_ON_BEADS = pBH.('autoAli_refine_on_beads');
catch
  REFINE_ON_BEADS = false;
end

  LOW_RES_CUTOFF=800;
  MIN_SAMPLING_RATE=10.0;
  ITERATIONS_PER_BIN=3;
  CLEAN_UP_RESULTS=false;
  MAG_OPTION=5;
  tiltAngleOffset=0.0;
  TILT_OPTION = 0;

  
  
if nargin > 6
  MAX_SAMPLING_RATE = str2double(varargin{1})
else
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
system(sprintf('ln -sf %s/%s %s.fixed',startDir,stackIN,baseName));
cd('../');



if ~strcmp(tiltExt,'.rawtlt')
  system(sprintf('ln -sf %s %s.rawtlt',tiltAngles,tiltName));
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

% gradientAliasFilter = BH_bandpass3d([nX,nY,1],1e-6,LOW_RES_CUTOFF,RESOLUTION_CUTOFF,'GPU',pixelSize);
                     gradientAliasFilter = {BH_bandpass3d(1.*[nX,nY,1],0,0,0,'GPU','nyquistHigh'),...
                       BH_bandpass3d([nX,nY,1],1e-6,LOW_RES_CUTOFF,RESOLUTION_CUTOFF,'GPU',pixelSize)};
medianFilter = 3;
for iPrj = 1:nZ
%   tmpPrj = BH_preProcessStack(gpuArray(inputStack(:,:,iPrj)),gradientAliasFilter,medianFilter);
  tmpPrj = real(ifftn(fftn(gpuArray(inputStack(:,:,iPrj))).*gradientAliasFilter{1}));
  tmpPrj = medfilt2(tmpPrj,medianFilter.*[1,1]);
  tmpPrj = real(ifftn(fftn(tmpPrj).*gradientAliasFilter{2}));
  tmpPrj = BH_resample2d(tmpPrj,rotMat,[0,0],'Bah','GPU','forward',1,size(tmpPrj));
  inputStack(:,:,iPrj) = gather(tmpPrj);
end
fprintf('finished preprocessing tilt-series\n');
SAVE_IMG(inputStack,fixedName,pixelSize);

clear tmpPrj inputStack


cd(wrkDir)

rotFile = fopen('preRotXf.xf','w');
fprintf(rotFile,'%6.3f %6.3f %6.3f %6.3f 0.0 0.0\n',rotMat);
fclose(rotFile);


                                             
                                     
system('pwd')
fprintf('Running %s\n',runPath);
system(sprintf('%s %s %f %f %d %d %d %d %d %d %s %d > ./emC_autoAliLog_%s.txt',...
                                                   runPath, ...
                                                   baseName, ...
                                                   pixelSize, ...
                                                   imgRotation, ...
                                                   binHigh, ...
                                                   binLow, ...
                                                   binInc,...
                                                   nX,nY,nZ,ext,...
                                                   PATCH_SIZE_FACTOR,baseName));
    
cd(sprintf('%s',startDir));

system(sprintf('rm %s',fixedName));

if strcmpi(TILT_OPTION,'0')
  % If not fitting tilt angles we need a copy of them with .tlt
   system(sprintf('cp %s fixedStacks/%s.tlt',tiltAngles,tiltName));
end

if (REFINE_ON_BEADS)
  cd(sprintf('%s',wrkDir'));
  extList = {'fixed','tlt','xf','local'}; % stack is skipped in second round. leave as number 1
  for iExt = 1:length(extList)
    system(sprintf('ln -sf ../fixedStacks/%s.%s %s.%s', ...
                    baseName,extList{iExt},baseName,extList{iExt}));
  end  
  
  % Stopping for now at a bin5, this should be dynamic along with a handful
  % of other options.
  max_binning = 5;
  BH_refine_on_beads(baseName,nX,nY,3000,pixelSize,1.05.*100);
  
  for iExt = 2:length(extList)
    system(sprintf('mv ../fixedStacks/%s.%s ../fixedStacks/%s.%s_patchTracking', ...
                    baseName,extList{iExt},baseName,extList{iExt}));
    system(sprintf('cp %s_%d.%s ../fixedStacks/%s.%s', ...
                    baseName,max_binning,extList{iExt},baseName,extList{iExt}));                  
  end   
  
  
  system(sprintf('imodtrans -i ../fixedStacks/%s.fixed %s_%d_fit.fid ../fixedStacks/%s.erase',...
                 baseName,baseName,max_binning,baseName));
  cd ..
  
else
  cd fixedStacks
  system(sprintf('%s %s %d %d %d %d', findBeadsPath, baseName,...
                                         nX,nY,3000,...
                                         ceil(1.05*100/pixelSize)));
  cd ..
end

end

