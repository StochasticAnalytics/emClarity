function [ ] = BH_runAutoAlign(PARAMETER_FILE, runPath,findBeadsPath,stackIN,tiltAngles,imgRotation,varargin)
%Run auto tiltseries alignment. For now, basically a script
%  sadf

% TODO add options for experimenting.
emc = BH_parseParameterFile(PARAMETER_FILE);

skip_tilts = 0;
if nargin > 6
  skip_tilts = EMC_str2double(varargin{1});
end

imgRotation = EMC_str2double(imgRotation);

try
  RESOLUTION_CUTOFF = emc.('autoAli_max_resolution');
catch
  RESOLUTION_CUTOFF=18;
end

% Min and max sampling rate in Ang/Pix (for patch tracking)
try
  MIN_SAMPLING_RATE = emc.('autoAli_min_sampling_rate');
catch
  MIN_SAMPLING_RATE = 10.0;
end
try
  MAX_SAMPLING_RATE = emc.('autoAli_max_sampling_rate');
catch
  MAX_SAMPLING_RATE = 4.0;
end

try
  PATCH_SIZE_FACTOR = emc.('autoAli_patch_size_factor');
catch
  PATCH_SIZE_FACTOR = 4;
end

% Check this first to allow only patch tracking even if there are beads
try
  REFINE_ON_BEADS = emc.('autoAli_refine_on_beads');
catch
  REFINE_ON_BEADS = false;
end

% Check this first to allow only patch tracking even if there are beads
try
  BORDER_SIZE_PIXELS = emc.('autoAli_patch_tracking_border');
catch
  BORDER_SIZE_PIXELS = 64;
end

% Check this first to allow only patch tracking even if there are beads
try
  N_ITERS_NO_ROT = emc.('autoAli_n_iters_no_rotation');
catch
  N_ITERS_NO_ROT = 3;
end

try
  PATCH_OVERLAP = emc.('autoAli_patch_overlap');
catch
  PATCH_OVERLAP = 0.5;
end

try
  ITERATIONS_PER_BIN = emc.('autoAli_iterations_per_bin');
catch
  ITERATIONS_PER_BIN = 3;
end

% FIXME this should probably be specified in Ang
try
  FIRST_ITER_SHIFT_LIMIT_PIXELS  =  ceil(emc.('autoAli_max_shift_in_angstroms')./emc.pixel_size_angstroms);
catch
  FIRST_ITER_SHIFT_LIMIT_PIXELS = ceil(40 ./ emc.pixel_size_angstroms);
end

try
  DIVIDE_SHIFT_LIMIT_BY =  emc.('autoAli_max_shift_factor');
catch
  DIVIDE_SHIFT_LIMIT_BY = 1;
  % int(max_shift / (iter^DIVI...)) + 1
end

% Now get the bead diameter, if it is zeros override the default to refine
% on beads after patch tracking.
beadDiameter = emc.('beadDiameter') * 10^10;
if beadDiameter == 0
  REFINE_ON_BEADS = false;
end



LOW_RES_CUTOFF=800;
CLEAN_UP_RESULTS=false;
MAG_OPTION=5;
tiltAngleOffset=0.0;
TILT_OPTION = 0;

fprintf("Stack in is %s\n",stackIN);
inputMRC = MRCImage(stackIN,0);
inputStack = single(getVolume(inputMRC));

skip_tilts_logical = [];
if (skip_tilts)
  skip_tilts_logical = ~ismember(1:size(inputStack,3),skip_tilts);
  inputStack = inputStack(:,:,skip_tilts_logical);
else
  skip_tilts_logical = true(size(inputStack,3),1);
end

[~,baseName,ext] = fileparts(stackIN);
fixedName = sprintf('fixedStacks/%s.fixed.preprocessed',baseName);

ext = sprintf('%s.preprocessed',ext);

startDir = pwd;
wrkDir = sprintf('emC_autoAlign_%s', baseName);
system(sprintf('mkdir -p %s',wrkDir));
system(sprintf('mkdir -p fixedStacks'));

cd('fixedStacks');

% Get this info here, as it is possibly used in the skip_tilts branch or in
% the rotate to avoid information loss branch
iHeader = getHeader(inputMRC);
iPixelHeader = [iHeader.cellDimensionX/iHeader.nX, ...
  iHeader.cellDimensionY/iHeader.nY, ...
  iHeader.cellDimensionZ/iHeader.nZ];

iOriginHeader= [iHeader.xOrigin , ...
  iHeader.yOrigin , ...
  iHeader.zOrigin ];

f = load(sprintf('../%s',tiltAngles));
f = f(skip_tilts_logical);
fout = fopen(sprintf('%s.rawtlt',baseName),'w');
fprintf(fout,'%f\n',f');
fclose(fout);
clear f

cd('../');


binHigh=ceil(MIN_SAMPLING_RATE ./ emc.pixel_size_angstroms);
if MAX_SAMPLING_RATE > 4
  binLow = MAX_SAMPLING_RATE;
else
  binLow = ceil(MAX_SAMPLING_RATE ./ emc.pixel_size_angstroms);
end
binInc = -1*ceil((binHigh- binLow)./3);


[nX,nY,nZ] = size(inputStack);


% % Check for rotations close to 90 or 270 that would reduce the total
% % useable area, and if found switch the x/y dimension. This could be more
% % precise
% maxAngle = atand(nY./nX)
% ; % will be positive
% switch_axes = false;
% abs(abs(imgRotation) - 180)

a = ones(nX,nY,'single','gpuArray');
p = BH_multi_padVal([nX,nY],max([nX,nY]).*[2,2]);
pad = BH_padZeros3d(a,'fwd',p,'GPU','single');

b = BH_resample2d(pad,[imgRotation,0,0],[0,0],'Bah','GPU','inv',1,size(pad));
s = pad+b;
score_1 = sum(sum(s==2))./sum(b(:));

pad = rot90(pad);
b = BH_resample2d(pad,[90-imgRotation,0,0],[0,0],'Bah','GPU','forward',1,size(pad));
s = pad+b;
score_2 = sum(sum(s==2))./sum(b(:));

if score_2 > score_1
  switch_axes = true;
else
  switch_axes = false;
end


if (switch_axes)
  
  % if ( abs(abs(imgRotation) - 180) > maxAngle )
  fprintf('Your image rotation will result in a loss of data. Switching X/Y axes\n')
  %   switch_axes = true;
  
  rotStack = zeros(nY,nX,nZ,'single');
  %
  ny = nY;
  nY = nX;
  nX = ny;
  tmpFile = sprintf('%s/%s_tmp.st',getenv('MCR_CACHE_ROOT'),baseName);
  for iPrj = 1:nZ
    % Once we've done this, we want to work as if this is how the stack
    % came off the scope.
    system(sprintf('newstack -fromone -secs %d -rotate 90 fixedStacks/%s.fixed %s >/dev/null',iPrj,baseName,tmpFile));
    rotStack(:,:,iPrj) = getVolume(MRCImage(sprintf('%s',tmpFile)));
    
    system(sprintf('rm %s',tmpFile));
  end
  inputStack = rotStack; clear rotStack
  imgRotation = imgRotation + 90;
  SAVE_IMG(inputStack,sprintf('fixedStacks/%s.fixed',baseName),iPixelHeader,iOriginHeader);
elseif ( skip_tilts)
  % Originally saved in the skip_tilts block, but that is redundant if we
  % save in the switch_axes block in the new implementation.
  SAVE_IMG(inputStack,sprintf('fixedStacks/%s.fixed',baseName),iPixelHeader,iOriginHeader);
else
  % No modifications, so just link to the original stack
  cd('fixedStacks');
  system(sprintf('ln -sf ../%s %s.fixed',stackIN,baseName));
  cd('../')
end


rotMat = [cosd(imgRotation),-1*sind(imgRotation),sind(imgRotation),cosd(imgRotation)];

fprintf('Preprocessing tilt-series\n');

%gradientAliasFilter = BH_bandpass3d([nX,nY,1],1e-6,LOW_RES_CUTOFF,RESOLUTION_CUTOFF,'GPU',emc.pixel_size_angstroms);
gradientAliasFilter = {BH_bandpass3d(1.*[nX,nY,1],0,0,0,'GPU','nyquistHigh'),...
  BH_bandpass3d([nX,nY,1],1e-6,LOW_RES_CUTOFF,RESOLUTION_CUTOFF,'GPU',emc.pixel_size_angstroms)};
if emc.pixel_size_angstroms < 2
  medianFilter = 5;
else
  medianFilter = 3;
end

for iPrj = 1:nZ
  %   tmpPrj = BH_preProcessStack(gpuArray(inputStack(:,:,iPrj)),gradientAliasFilter,medianFilter);
  tmpPrj = real(ifftn(fftn(gpuArray(inputStack(:,:,iPrj))).*gradientAliasFilter{1}));
  tmpPrj = medfilt2(tmpPrj,medianFilter.*[1,1]);
  tmpPrj = real(ifftn(fftn(tmpPrj).*gradientAliasFilter{2}));
  %   tmpPrj = BH_resample2d(tmpPrj,rotMat,[0,0],'Bah','GPU','inv',1,size(tmpPrj));
  inputStack(:,:,iPrj) = gather(tmpPrj);
end

SAVE_IMG(inputStack,fixedName,emc.pixel_size_angstroms);
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
  emc.pixel_size_angstroms, ...
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

% 2021-May-08 BAH, not needed b/c fixed/name.fixed should remain rotated by
% 90
% if (switch_axes)
%   % backup the original
%   system(sprintf('mv fixedStacks/%s.fixed fixedStacks/%s.fixed_nonSwapped',baseName,baseName));
%   % rotate by 90
%   system(sprintf('newstack -rotate 90 fixedStacks/%s.fixed_nonSwapped fixedStacks/%s.fixed',baseName,baseName));
% end

if strcmpi(TILT_OPTION,'0')
  % If not fitting tilt angles we need a copy of them with .tlt
  system(sprintf('cp fixedStacks/%s.rawtlt fixedStacks/%s.tlt',baseName,baseName));
end

to_few_beads = false;

if (REFINE_ON_BEADS)
  cd(sprintf('%s',wrkDir'));
  extList = {'tlt','xf','local'}; % stack is skipped in second round. leave as number 1
  for iExt = 1:length(extList)
    system(sprintf('ln -sf ../fixedStacks/%s.%s %s.%s', ...
      baseName,extList{iExt},baseName,extList{iExt}));
  end
  
  system(sprintf('ln -sf ../fixedStacks/%s.%s.preprocessed %s.%s', ...
    baseName,'fixed',baseName,'fixed'));
  
  
  % Stopping for now at a bin5, this should be dynamic along with a handful
  % of other options.
  min_sampling_rate = 5;
  [ to_few_beads ] = BH_refine_on_beads(baseName,nX,nY,3000,emc.pixel_size_angstroms,1.05.*100, min_sampling_rate);
  
  if (to_few_beads)
    fprintf('\nWARNING: to few beads found. Using iterative patch tracking results\n');
  else
    for iExt = 1:length(extList)
      system(sprintf('mv ../fixedStacks/%s.%s ../fixedStacks/%s.%s_patchTracking', ...
        baseName,extList{iExt},baseName,extList{iExt}));
      system(sprintf('cp %s_%d.%s ../fixedStacks/%s.%s', ...
        baseName,min_sampling_rate,extList{iExt},baseName,extList{iExt}));
    end
    
    
    system(sprintf('imodtrans -i ../fixedStacks/%s.fixed %s_%d_fit.fid ../fixedStacks/%s.erase',...
      baseName,baseName,min_sampling_rate,baseName));
    
    system(sprintf('newstack -xf ../fixedStacks/%s.xf -bin 12 ../fixedStacks/%s.fixed ../fixedStacks/%s_bin12.ali',baseName,baseName,baseName));
  end
  
  cd ..
end

if (to_few_beads || ~REFINE_ON_BEADS)
  cd fixedStacks
  system(sprintf('%s %s %d %d %d %d', findBeadsPath, baseName,...
    nX,nY,3000,...
    ceil(1.05*100/emc.pixel_size_angstroms)));
  cd ..
end

system(sprintf('rm %s',fixedName));

end

