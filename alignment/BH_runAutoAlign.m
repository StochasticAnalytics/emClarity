function [ ] = BH_runAutoAlign(PARAMETER_FILE, runPath,findBeadsPath,stackIN,tiltAngles,imgRotation,varargin)
%Run auto tiltseries alignment. For now, basically a script
%  sadf

% TODO add options for experimenting.
pBH = BH_parseParameterFile(PARAMETER_FILE);

skip_tilts = 0;
if nargin > 6
  skip_tilts = str2num(varargin{1});  
end

pixelSize = pBH.('PIXEL_SIZE').*10^10;
imgRotation = str2double(imgRotation);

try 
  RESOLUTION_CUTOFF = pBH.('autoAli_max_resolution');
catch
  RESOLUTION_CUTOFF=18;
end

% Min and max sampling rate in Ang/Pix (for patch tracking)
try
  MIN_SAMPLING_RATE = pBH.('autoAli_min_sampling_rate');
catch
  MIN_SAMPLING_RATE = 10.0;
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

% Check this first to allow only patch tracking even if there are beads
try
  REFINE_ON_BEADS = pBH.('autoAli_refine_on_beads');
catch
  REFINE_ON_BEADS = false;
end

% Check this first to allow only patch tracking even if there are beads
try
  BORDER_SIZE_PIXELS = pBH.('autoAli_patch_tracking_border');
catch
  BORDER_SIZE_PIXELS = 64;
end

% Check this first to allow only patch tracking even if there are beads
try
  N_ITERS_NO_ROT = pBH.('autoAli_n_iters_no_rotation');
catch
  N_ITERS_NO_ROT = 3;
end

try
  PATCH_OVERLAP = pBH.('autoAli_patch_overlap');
catch
  PATCH_OVERLAP = 0.5;
end

try
  ITERATIONS_PER_BIN = pBH.('autoAli_iterations_per_bin');
catch
  ITERATIONS_PER_BIN = 3;
end

% FIXME this should probably be specified in Ang
try
  FIRST_ITER_SHIFT_LIMIT_PIXELS  =  ceil(pBH.('autoAli_max_shift_in_angstroms')./pixelSize);   
catch
  FIRST_ITER_SHIFT_LIMIT_PIXELS = ceil(40 ./ pixelSize);
end

try
  DIVIDE_SHIFT_LIMIT_BY =  pBH.('autoAli_max_shift_factor');
catch
  DIVIDE_SHIFT_LIMIT_BY = 1;
  % int(max_shift / (iter^DIVI...)) + 1
end
                                                   
% Now get the bead diameter, if it is zeros override the default to refine
% on beads after patch tracking.
beadDiameter = pBH.('beadDiameter') * 10^10;
if beadDiameter == 0
  REFINE_ON_BEADS = false;
end



  LOW_RES_CUTOFF=800;
  CLEAN_UP_RESULTS=false;
  MAG_OPTION=5;
  tiltAngleOffset=0.0;
  TILT_OPTION = 0;

inputMRC = MRCImage(stackIN,0);
inputStack = single(getVolume(inputMRC));

skip_tilts_logical = [];
if (skip_tilts)
  skip_tilts_logical = ~ismember(1:size(inputStack,3),skip_tilts);
  inputStack = inputStack(:,:,skip_tilts_logical);
end

[~,baseName,ext] = fileparts(stackIN);
fixedName = sprintf('fixedStacks/%s.fixed.rot',baseName);

ext = sprintf('%s.rot',ext);

[~,tiltName,tiltExt] = fileparts(tiltAngles);
startDir = pwd;
wrkDir = sprintf('emC_autoAlign_%s', baseName);
system(sprintf('mkdir -p %s',wrkDir));
system(sprintf('mkdir -p fixedStacks'));

cd('fixedStacks');

if (skip_tilts)
  iHeader = getHeader(inputMRC);
  iPixelHeader = [iHeader.cellDimensionX/iHeader.nX, ...
                  iHeader.cellDimensionY/iHeader.nY, ...
                  iHeader.cellDimensionZ/iHeader.nZ];
                
  iOriginHeader= [iHeader.xOrigin , ...
                  iHeader.yOrigin , ...
                  iHeader.zOrigin ];
  SAVE_IMG(inputStack,sprintf('%s.fixed',baseName),iPixelHeader,iOriginHeader);
  
  f = load(sprintf('../%s',tiltAngles));
  if length(tiltAngles) ~= sum(skip_tilts_logical)
    if exist(sprintf('../%s.orig',tiltAngles),'file')
      % We must be re_rerunning, so copy the orig back
      system(sprintf('cp ../%s.orig ../%s',tiltAngles,tiltAngles));
      f = load(sprintf('../%s',tiltAngles));
    else
      % Create a backup 
      system(sprintf('cp ../%s ../%s.orig',tiltAngles,tiltAngles));
    end
    f = f(skip_tilts_logical);
    fout = fopen(sprintf('../%s',tiltAngles),'w');
    fprintf(fout,'%s\n',f');
    fclose(fout);
  end
  clear f
      
else
  system(sprintf('ln -sf %s/%s %s.fixed',startDir,stackIN,baseName));
end

cd('../');


if strcmp(tiltExt,'.rawtlt')
  system(sprintf('ln -sf %s %s.rawtlt',tiltAngles,tiltName));
end

binHigh=ceil(MIN_SAMPLING_RATE ./ pixelSize);
if MAX_SAMPLING_RATE > 4
  binLow = MAX_SAMPLING_RATE;
else
  binLow = ceil(MAX_SAMPLING_RATE ./ pixelSize);
end
binInc = -1*ceil((binHigh- binLow)./3);


[nX,nY,nZ] = size(inputStack);


% Check for rotations close to 90 or 270 that would reduce the total
% useable area, and if found switch the x/y dimension. This could be more
% precise
maxAngle = atand(nY./nX); % will be positive
switch_axes = false;

if ( abs(abs(imgRotation) - 180) > maxAngle )
  fprintf('Your image rotation will result in a loss of data. Switching X/Y axes\n')
  
  switch_axes = true;
  
  rotStack = zeros(nY,nX,nZ,'single');
%   
  ny = nY;
  nY = nX;
  nX = ny;
  tmpFile = sprintf('%s/%s_tmp.st',getenv('MCR_CACHE_ROOT'),baseName);
  for iPrj = 1:nZ
    
    system(sprintf('newstack -fromone -secs %d -rotate 90 fixedStacks/%s.fixed %s >/dev/null',iPrj,baseName,tmpFile));
    rotStack(:,:,iPrj) = getVolume(MRCImage(sprintf('%s',tmpFile)));

    system(sprintf('rm %s',tmpFile));
  end
  
  inputStack = rotStack; clear rotStack
  imgRotation = imgRotation + 90;
end


rotMat = [cosd(imgRotation),-1*sind(imgRotation),sind(imgRotation),cosd(imgRotation)];

fprintf('Preprocessing tilt-series\n');

%gradientAliasFilter = BH_bandpass3d([nX,nY,1],1e-6,LOW_RES_CUTOFF,RESOLUTION_CUTOFF,'GPU',pixelSize);
                      gradientAliasFilter = {BH_bandpass3d(1.*[nX,nY,1],0,0,0,'GPU','nyquistHigh'),...
                        BH_bandpass3d([nX,nY,1],1e-6,LOW_RES_CUTOFF,RESOLUTION_CUTOFF,'GPU',pixelSize)};
if pixelSize < 2                     
  medianFilter = 5;
else
  medianFilter = 3;
end

for iPrj = 1:nZ
%   tmpPrj = BH_preProcessStack(gpuArray(inputStack(:,:,iPrj)),gradientAliasFilter,medianFilter);
  tmpPrj = real(ifftn(fftn(gpuArray(inputStack(:,:,iPrj))).*gradientAliasFilter{1}));
  tmpPrj = medfilt2(tmpPrj,medianFilter.*[1,1]);
  tmpPrj = real(ifftn(fftn(tmpPrj).*gradientAliasFilter{2}));
  tmpPrj = BH_resample2d(tmpPrj,rotMat,[0,0],'Bah','GPU','forward',1,size(tmpPrj));
  inputStack(:,:,iPrj) = gather(tmpPrj);
end

SAVE_IMG(inputStack,fixedName,pixelSize);
fprintf('finished preprocessing tilt-series\n');

clear tmpPrj inputStack


cd(wrkDir)

rotFile = fopen('preRotXf.xf','w');
fprintf(rotFile,'%6.3f %6.3f %6.3f %6.3f 0.0 0.0\n',rotMat);
fclose(rotFile);


                                            
                                     
system('pwd')
fprintf('Running %s\n',runPath);
system(sprintf('%s %s %f %f %d %d %d %d %d %d %s %d %d %d %f %f %f %d %d %d > ./emC_autoAliLog_%s.txt',...
                                                   runPath, ...
                                                   baseName, ...
                                                   pixelSize, ...
                                                   imgRotation, ...
                                                   binHigh, ...
                                                   binLow, ...
                                                   binInc,...
                                                   nX,nY,nZ,ext,...
                                                   PATCH_SIZE_FACTOR,...
                                                   N_ITERS_NO_ROT,...
                                                   BORDER_SIZE_PIXELS,...
                                                   PATCH_OVERLAP,...
                                                   RESOLUTION_CUTOFF,...
                                                   LOW_RES_CUTOFF,...
                                                   ITERATIONS_PER_BIN,...
                                                   FIRST_ITER_SHIFT_LIMIT_PIXELS,...
                                                   DIVIDE_SHIFT_LIMIT_BY,...
                                                   baseName));
    
cd(sprintf('%s',startDir));

system(sprintf('rm %s',fixedName));


if (switch_axes)
  % backup the original
  system(sprintf('mv fixedStacks/%s.fixed fixedStacks/%s.fixed_nonSwapped',baseName,baseName));
  % rotate by 90
  system(sprintf('newstack -rotate 90 fixedStacks/%s.fixed_nonSwapped fixedStacks/%s.fixed',baseName,baseName));
end

if strcmpi(TILT_OPTION,'0')
  % If not fitting tilt angles we need a copy of them with .tlt
   system(sprintf('cp %s fixedStacks/%s.tlt',tiltAngles,tiltName));
end

to_few_beads = false;

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
  [ to_few_beads ] = BH_refine_on_beads(baseName,nX,nY,3000,pixelSize,1.05.*100);
  
  if (to_few_beads)
    fprintf('\nWARNING: to few beads found. Using iterative patch tracking results\n');
  else
    for iExt = 2:length(extList)
      system(sprintf('mv ../fixedStacks/%s.%s ../fixedStacks/%s.%s_patchTracking', ...
                      baseName,extList{iExt},baseName,extList{iExt}));
      system(sprintf('cp %s_%d.%s ../fixedStacks/%s.%s', ...
                      baseName,max_binning,extList{iExt},baseName,extList{iExt}));                  
    end   
  
  
    system(sprintf('imodtrans -i ../fixedStacks/%s.fixed %s_%d_fit.fid ../fixedStacks/%s.erase',...
                 baseName,baseName,max_binning,baseName));
               
    system(sprintf('newstack -xf ../fixedStacks/%s.xf -bin 12../fixedStacks/%s.fixed ../fixedStacks/%s_bin12.ali',baseName,baseName,baseName));
  end
  
  cd ..
end

if (to_few_beads || ~REFINE_ON_BEADS)
  cd fixedStacks
  system(sprintf('%s %s %d %d %d %d', findBeadsPath, baseName,...
                                         nX,nY,3000,...
                                         ceil(1.05*100/pixelSize)));
  cd ..
end

end

