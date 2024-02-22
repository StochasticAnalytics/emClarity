function [ ] = BH_geometryInitialize( PARAMETER_FILE, varargin )
%Take a set of coordinates and optionally angles and set up geometry.
%   For now just BH, but also to work with Relion set up.

% Option for testing euler angles will remove.
global bh_global_do_profile;
if isempty(bh_global_do_profile)
  bh_global_do_profile = false;
end

convertEulers = 0;
direction = '';

if (bh_global_do_profile)
  profile on;
end

% Initialize the subTomoMeta structure
subTomoMeta = struct();
mapBackIter = 0;
doImport=false;
angleSgn=1;
convention='Bah';
direction='fwd';
emc = BH_parseParameterFile(PARAMETER_FILE);
if nargin > 1
  if length(varargin) == 1
    mapBackIter = EMC_str2double(varargin{1});
  else
    
    eulerSet = {'zyx','zxz','zyz'};
    defMatSet= {'IMOD','Bah','SPIDER'};
    doImport = true;
    convention = varargin{1};
    direction = varargin{2};
    angleSgn = 1;
    
    conventionMessage = ['\nThere are four ways to use each euler convention :\n ',...
      'zyz fwd, zyz inv, -1*(zyz) fwd, -1*(zyz) inv \n', ...
      'the default is positive, but adding the addidional neg argument will produce the last two\n',...
      '\n\n Assuming the zyz,zxz,zyx you selected is correct, one of these should work.\n', ...
      'PLEASE let me know if you have confirmed a convention for a particluar software and ',...
      'how you made extracted the angles and I will add a software specific flag. Thank you!\n'];
    
    switch convention
      case 'Protomo'
        fprintf('Using the Protomo convention\n\tZYZ, fwd, neg\n');
        convention = 'Bah',
        direction = 'forward';
        angleSgn = -1;
      case 'Dynamo'
        fprintf('%s',conventionMessage);
        error('For Dynamo please try zyz fwd first');
      case 'Relion'
        fprintf('%s',conventionMessage);
        error('For Releion please try zyz inv first');
      case 'Eman2'
        fprintf('%s',conventionMessage);
        error('For Eman2 please try zyz fwd first');
      case 'Peet'
        fprintf('%s',conventionMessage);
        error('For Peet please try zxz fwd first');
      case 'Slicer'
        fprintf('%s',conventionMessage);
        error('For Imod slicer please try zyz fwd first');
      otherwise
        if ismember(convention,eulerSet)
          convention = defMatSet{find(ismember(eulerSet,convention))};
          if (strcmpi(direction,'fwd'))
            direction = 'forward';
          elseif ~(strcmpi(direction,'inv'))
            error('angular convention direction must be fwd or inv');
          end
        else
          error('Not found')
        end
    end
    
    if length(varargin) == 3
      if strcmpi(varargin{3},'neg')
        angleSgn = -1;
      else
        error('If supplying a third angular argument, it must be neg to negate the angles in your csv');
      end
    end
    
    
  end
end
nGPUs = emc.('nGPUs');
% This will have to do until a better approach based on PSF of positions
% in projection space linked together can be used to define groups that
% don't have co-mingled resolution. When splitOnTOmos, always run in serial, so
% that we have the best chance of diTomoributing the defocus variateion/tomoqualtiy.
splitOnTomos = emc.('fscGoldSplitOnTomos');
if (splitOnTomos)
  nGPUs = 1;
  fprintf('override nGPUs to just 1 for initial step to evenly split crowded tomos because fscGoldSplitOnTomos is true')
end

% Resolution lower than this is not gold standard, and will also be mixed
% in the references to keep orientations from diverging. Should be > 2.25 x
% your expected resolution. 1.5x to be totally uncorrelated, and another
% 1.5x because aligning at ~40A in template matching will (with good data)
% often lead to and FSC of ~ 26A initially.

% 40 may be too conservative, especially for low defocus tomos
try
  lowResCut = emc.('lowResCut');
catch
  lowResCut = 40;
end

maxGoldStandard = lowResCut;%emc.('Tmp_bandpassFilter')(3);
dupSampling = emc.('Tmp_samplingRate')
dupRadius = 2;
dupTolerance = (2.*dupRadius)+1;
dupMask = zeros(dupTolerance, dupTolerance, dupTolerance, 'single');
dupMask = dupMask + 1;

% If we are in the working directory, following template matching, there should
% be the director convmap, holding convolution maps, model files etc.

if ~isdir('convmap')
  error('Did not find directory named <convmap>')
end

% Each tomogram will have its own file suffix .path which holds the basename,
% path, file extension, tilt file. Get these and concat into a cell.

if (doImport)
  getPath = dir('convmap/*.csv');
else
  getPath = dir('convmap/*.mod');
end

nTomogramsTotal = length(getPath);
fileInfo = cell(nTomogramsTotal,4);
n_tomos_added = 1;
getCoords = dir('recon/*.coords');
nStacks = length(getCoords);

% There is no real need to have a separate data here. The main difference
% is these have all tomos from one tilt-series and instead of the origin on
% Y, they liTomo the slices. Should convert to include both values and also
% include the size of the tilt-series which could then be removed from the
% TLT geometry. I think it important to not duplicate the information as
% this could lead to bugs if one is changed and the other not. The only
% other concern is then linking each tomogram to the parent tilt-series.
for iStack = 1:nStacks
  
  [ recGeom, tiltName, nTomosPossible, tilt_geometry ] = BH_multi_recGeom( sprintf('recon/%s',getCoords(iStack).name), mapBackIter);
  % Initialize
  
  if (doImport)
    iPath = dir(sprintf('convmap/%s_*.csv',tiltName));
  else
    iPath = dir(sprintf('convmap/%s_*.mod',tiltName));
  end
  
  nTomos = length(iPath);
  if nTomos > nTomosPossible
    error('The number of model files in convmap/*.mod is greater than the number in the recon/*.coords\n');
  elseif nTomos < nTomosPossible
    fprintf('\n\n\tThere are fewer model files in convmap for %s than are described in your coords file.\n',tiltName);
    fprintf('This is okay, but take note this is what you intended.\n');
  end
  subTomoMeta.('mapBackGeometry').(tiltName).('nTomos') = nTomos;
  subTomoMeta.('mapBackGeometry').(tiltName).('tomoCprRePrjSize') = 512;

  for iTomo = 1:nTomos
    
    if (doImport)
      modName = strsplit(iPath(iTomo).name,'.csv');
      modName = strsplit(modName{1},'_');
      tomoIdx = EMC_str2double(modName{2});
    else
      modName = strsplit(iPath(iTomo).name,'_');
      tomoIdx = EMC_str2double(modName{end-1});
    end

    % We are storing this info to make it available when checking for duplicates
    fileInfo{n_tomos_added,1} = tiltName;
    fileInfo{n_tomos_added,2} = sprintf('%s_%d', tiltName, tomoIdx);
    fileInfo{n_tomos_added,3} = sprintf('%s_%d_bin%d',tiltName, tomoIdx, dupSampling);
    tomoName = sprintf('%s_%d',tiltName,tomoIdx);
    
    subTomoMeta.('tiltGeometry').(fileInfo{n_tomos_added,2}) = tilt_geometry;
    n_tomos_added = n_tomos_added + 1;
    
    % Store a reference to the parent tilt-series for every tomogram
    subTomoMeta.('mapBackGeometry').('tomoName').(tomoName).('tiltName') = tiltName;
    % Store the tomoIdx for every tomogram, currently used to refer back to recGEom, but I'm going to put this into a struct
    subTomoMeta.('mapBackGeometry').('tomoName').(tomoName).('tomoIdx') = tomoIdx;
    subTomoMeta.('mapBackGeometry').('tomoCoords').(tomoName).('is_active') = true;
    subTomoMeta.('mapBackGeometry').('tomoCoords').(tomoName).('y_i') = recGeom{tomoIdx}.y_i;
    subTomoMeta.('mapBackGeometry').('tomoCoords').(tomoName).('y_f') = recGeom{tomoIdx}.y_f;
    subTomoMeta.('mapBackGeometry').('tomoCoords').(tomoName).('NX')  = recGeom{tomoIdx}.NX;
    subTomoMeta.('mapBackGeometry').('tomoCoords').(tomoName).('NY')  = recGeom{tomoIdx}.NY;
    subTomoMeta.('mapBackGeometry').('tomoCoords').(tomoName).('NZ')  = recGeom{tomoIdx}.NZ;
    subTomoMeta.('mapBackGeometry').('tomoCoords').(tomoName).('dX_specimen_to_tomo') = recGeom{tomoIdx}.dX_specimen_to_tomo;
    subTomoMeta.('mapBackGeometry').('tomoCoords').(tomoName).('dY_specimen_to_tomo') = recGeom{tomoIdx}.dY_specimen_to_tomo;
    subTomoMeta.('mapBackGeometry').('tomoCoords').(tomoName).('dZ_specimen_to_tomo') = recGeom{tomoIdx}.dZ_specimen_to_tomo;
  end 
end % end of loop over stacks

% For now, just assuming all of the maps are in the same place and have the same
% suffix - generalize later.

if nGPUs > nTomogramsTotal
  nGPUs = nTomogramsTotal;
end


iterLiTomo = cell(nGPUs,1);
for iGPU = 1:nGPUs
  iterLiTomo{iGPU} = iGPU:nGPUs:nTomogramsTotal;
  iterLiTomo{iGPU}
end

try
  EMC_parpool(nGPUs)
catch
  delete(gcp('nocreate'))
  EMC_parpool(nGPUs)
end

parResults = cell(nGPUs,1);

dupInTheLoop = dupSampling
parfor iGPU = 1:nGPUs
  % for iGPU = 1:nGPUs % revert
  
  D = gpuDevice(iGPU);
  
  tomoResults = struct();
  for iTomo = iterLiTomo{iGPU}
    
    if (doImport)
      mapName = fileInfo{iTomo,2};
    else
      mapName = fileInfo{iTomo,3};
    end
    
    % Load in the template matching geometry for the tomogram, and the model file
    % which may (or may not) have been edited.
    if isempty(fileInfo{iTomo,2})
      fprintf("\n\tWarning, missing tomoname for fileInfo{%d,2}\n",iTomo);
      pause(1)
      continue
    end
    tomoIdx = subTomoMeta.mapBackGeometry.tomoName.(fileInfo{iTomo,2}).tomoIdx;
    tiltName   = subTomoMeta.mapBackGeometry.tomoName.(fileInfo{iTomo,2}).tiltName;
    try
      tmpSearchGeom = importdata(sprintf('convmap/%s.csv',mapName));
    catch
      error('Did not find the file convmap/%s.csv\n\nMake sure the binning in your parameter file is correct.',mapName)
    end
    
    if (doImport)
      tmpCSV = tmpSearchGeom;
      % TODO replace the 26 with the global for the number of lines FIXME
      tmpSearchGeom = zeros([size(tmpCSV,1),26],'single');
      tmpSearchGeom(:,[11:16]) = tmpCSV;
      
      for iAng = 1:size(tmpSearchGeom,1)
        tmpSearchGeom(iAng,17:25) = reshape( ...
          BH_defineMatrix( ...
          angleSgn .* tmpSearchGeom(iAng,14:16),...
          convention,direction),1,9);
      end
      
      % I'm assuming that the proper scaling was done, let the user know
      fprintf('\nImporting coordinates, assuming to be scaled properly to match full reconstruction size\n');
    end
    
    % Leave in for now, but check with Yunjie to remove for new import
    % style
    % New check for all -1 and then convert to Protomo
    if ( abs(sum(tmpSearchGeom(1,17:25)) + 9) < 1e-3 )
      fprintf('\nconverting euler angles from Protomo trf convention\n');
      for iAng = 1:size(tmpSearchGeom,1)
        
        tmpSearchGeom(iAng,17:25) = reshape( ...
          BH_defineMatrix( ...
          -1.* flip(tmpSearchGeom(iAng,14:16)),...
          'Bah','inv'),1,9);
        
        tmpSearchGeom(iAng,11:13) = tmpSearchGeom(iAng,11:13)+1;
      end
    end
    if (convertEulers)
      
      fprintf('converting euler angles to spider/relion convention\n');
      for iAng = 1:size(tmpSearchGeom,1)
        tmpSearchGeom(iAng,17:25) = reshape( ...
          BH_defineMatrix( ...
          convertEulers.* ...
          tmpSearchGeom(iAng,14:16),...
          'SPIDER',direction),1,9);
        
        tmpSearchGeom(iAng,11:13) = tmpSearchGeom(iAng,11:13)+1;
      end
    end
    % Note the model coordinates are already scaled down by the sampling rate used
    % in the template search, which is the same sampling used here for the
    % comparison.
    
    % Putting in a try/catch for synthetic data without template matching step
    try
      tiltName = sprintf('convmap/%s_convmap.mrc',mapName);
      iHeader = getHeader(MRCImage(tiltName));
      % first convert the imod model file to a temporary text file
      tmpFile = sprintf('tmp_%d.txt',iGPU);
      system(sprintf('model2point convmap/%s.mod %s', mapName, tmpFile))
      modGeom = load(tmpFile);
      system(sprintf('rm %s', tmpFile));
      flgLookForPoints = 1;
      
      % leave IDX in main memory because it is just for reference.
      sx = floor(iHeader.nX);
      sy = floor(iHeader.nY);
      sz = floor(iHeader.nZ);
      
    catch
      fprintf('did not find convmap, assuming synthetic data, with no deleted model points\n');
      flgLookForPoints = 0;
    end
    
    % Make sure nothing has gone wrong in translating the convmap to the
    % full size
    if (flgLookForPoints)
      
      positionMatrix = zeros(sx, sy, sz, 'single', 'gpuArray');
      positionIDX    = zeros(sx, sy, sz, 'uint32');
      
      % Make a volume with ones in the position of the centers of the tomos.
      for iSubTomo = 1:size(tmpSearchGeom,1)
        
        subTomoOrigin = fix(tmpSearchGeom(iSubTomo,11:13)./dupInTheLoop);
        if any(subTomoOrigin < 1 + dupRadius) || any([sx,sy,sz] < subTomoOrigin + dupRadius)
          tmpSearchGeom(iSubTomo,26:26:26*emc.nPeaks) = -9999;
        else
          positionMatrix(subTomoOrigin(1),subTomoOrigin(2),subTomoOrigin(3)) = 1;
          positionIDX(subTomoOrigin(1),subTomoOrigin(2),subTomoOrigin(3)) = ...
            tmpSearchGeom(iSubTomo, 4);
        end
      end % loop building position matrix
      
      for iSubTomo = 1:size(modGeom,1)
        
        subTomoOrigin = fix(modGeom(iSubTomo,:));
        if all(subTomoOrigin > 1) && all([sx,sy,sz] > subTomoOrigin)
          
          positionMatrix(subTomoOrigin(1),subTomoOrigin(2),subTomoOrigin(3)) = ...
            positionMatrix(subTomoOrigin(1),subTomoOrigin(2),subTomoOrigin(3)) + 1;
        end % sometimes a point gets moved out of bounds.
      end
      
      % Convolve positionmatrix with duplicate mask. Just in case there are
      % rounding errors at any point use convolution to check a neighborhood of
      % +/- 1 pixel.
      
      overlapMatrix = convn(positionMatrix, gpuArray(dupMask), 'same');
      
      idxLiTomo = positionIDX((overlapMatrix > 1));

      tomoResults.(fileInfo{iTomo,2}) = tmpSearchGeom(ismember(tmpSearchGeom(:,4), idxLiTomo),:);
    else
      
      tomoResults.(fileInfo{iTomo,2}) = tmpSearchGeom;
    end
    
    % Fix the retained points
  end
  parResults{iGPU} = tomoResults;
  
end

% By my convention, every particle should have a unique id.

nIDX = 1;
for iGPU = 1:nGPUs
  for iTomo = iterLiTomo{iGPU}
    mapName = fileInfo{iTomo,2};
    tmpGeom = parResults{iGPU}.(mapName);
    
    % Sort so that CTFs can be left in main mem, and only pulled when needed and only
    % once per round of alignment.
    tmpGeom = sortrows(tmpGeom,9);
    
    
    for iSubTomo = 1:size(tmpGeom,1)
      tmpGeom(iSubTomo, 4:26:26*emc.nPeaks) = nIDX;
      
      nIDX = nIDX + 1;
    end
        
    subTomoMeta.('cycle000').('geometry').(mapName) = tmpGeom;
  end
end

fprintf('nSubTomos initial = %d\n', (nIDX-1));

subTomoMeta.('nSubTomoInitial') = nIDX-1;

preFscSplit = gather(subTomoMeta);

% Randomly divide the data into half sets.
[ subTomoMeta ] = BH_fscSplit( preFscSplit, splitOnTomos, emc.nPeaks);
subTomoMeta.('currentCycle') = 0;
subTomoMeta.('currentTomoCPR') = mapBackIter;
subTomoMeta.('currentResForDefocusError') = lowResCut;
subTomoMeta.('maxGoldStandard') = maxGoldStandard;
save(sprintf('%s.mat', emc.('subTomoMeta')), 'subTomoMeta');


if (bh_global_do_profile)
  profile off;
  profsave;
end

end % end of initialize function

