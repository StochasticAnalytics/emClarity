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
% that we have the best chance of distributing the defocus variateion/tomoqualtiy.
splitOnTomos = emc.('fscGoldSplitOnTomos');
if (splitOnTomos)
  nGPUs = 1;
  fprintf('override nGPUs to just 1 for initial step to evenly split crowded tomos because fscGoldSplitOnTomos is true')
end
nOrientations=1;%nOrientations = emc.('pseudoMLnumber');
nCTFgroups = 9;

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

checkDir = exist('convmap');
if checkDir ~= 7 
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

getCoords = dir('recon/*.coords')
getCoords(1).name
nStacks = length(getCoords)

% There is no real need to have a separate data here. The main difference
% is these have all tomos from one tilt-series and instead of the origin on
% Y, they list the slices. Should convert to include both values and also
% include the size of the tilt-series which could then be removed from the
% TLT geometry. I think it important to not duplicate the information as
% this could lead to bugs if one is changed and the other not. The only
% other concern is then linking each tomogram to the parent tilt-series. 
for iStack = 1:nStacks
 


sprintf('recon/%s',getCoords(iStack).name)


[ recGeom, tiltName, nTomosPossible] = BH_multi_recGeom( sprintf('recon/%s',getCoords(iStack).name) );
  % Initialize
 
                                                       
  subTomoMeta.('mapBackGeometry').(tiltName).('coords') = zeros(nTomosPossible,6);
  if (doImport)
    iPath = dir(sprintf('convmap/%s_*.csv',tiltName))
  else
    iPath = dir(sprintf('convmap/%s_*.mod',tiltName));
  end
  
  nTomos = length(iPath);
  fprintf('nTomos %d by mods, nTomosPossible %d, by coords\n',nTomos,nTomosPossible);
  if nTomos > nTomosPossible
    error('The number of model files in convmap/*.mod is greater than the number in the recon/*.coords\n');
  elseif nTomos < nTomosPossible
    fprintf('\n\nThere are fewer model files in convmap for %s than are described in your coords file.\n',tiltName);
    fprintf('This is okay, but take note this is what you intended.\n');
  end
  subTomoMeta.('mapBackGeometry').(sprintf('%s',tiltName)).('nTomos') = nTomos;
  subTomoMeta.('mapBackGeometry').(sprintf('%s',tiltName)).('tomoCprRePrjSize') = 512;

  for iSt = 1:nTomos

    if (doImport)
      modName = strsplit(iPath(iSt).name,'.csv');
      modName = strsplit(modName{1},'_');
      tomoNumber = EMC_str2double(modName{2});
    else
      modName = strsplit(iPath(iSt).name,'_');
      tomoNumber = EMC_str2double(modName{end-1});
    end
    

    subTomoMeta.('mapBackGeometry').(tiltName).('coords')(tomoNumber,:) = ...
                                                        recGeom(tomoNumber,:);
    subTomoMeta.('mapBackGeometry').('tomoName').(...
                sprintf('%s_%d',tiltName,tomoNumber)).('tiltName') = tiltName;
    subTomoMeta.('mapBackGeometry').('tomoName').(...
                sprintf('%s_%d',tiltName,tomoNumber)).('tomoNumber') = tomoNumber;                
  end
      
    
                                                         
end


for iTomo = 1:nTomogramsTotal
  iTomo

  sprintf('convmap/%s',getPath(iTomo).name)
  modName = strsplit(getPath(iTomo).name,'_')
  if (doImport)
    tiltName = modName{1}
    tomoNumber = strsplit(modName{2},'.csv');
    tomoNumber = EMC_str2double(tomoNumber{1})
  else
    tiltName = strjoin(modName(1:end-2),'_');
    tomoNumber = EMC_str2double(modName{end-1})
  end

  fileInfo{iTomo,1} = tiltName;
  fileInfo{iTomo,2} = sprintf('%s_%d',tiltName,tomoNumber);
  fileInfo{iTomo,3} = sprintf('%s_%d_bin%d',tiltName,tomoNumber,dupSampling);
  fileInfo{iTomo,4} = sprintf('fixedStacks/ctf/%s_ali%d_ctf.tlt',tiltName,mapBackIter+1);
  
   subTomoMeta.('tiltGeometry').(fileInfo{iTomo,2}) = load(fileInfo{iTomo,4});
   
    recCoords = importdata(sprintf('./recon/%s_recon.coords',tiltName))
    recCoords = recCoords.data
    
    % The reconstruction could be defined based on the aliStacks or the
    % fixedStacks. the dimensions 
    
    recGeom = [recCoords(2 + (tomoNumber-1)*6), ... % NX
               recCoords(4 + (tomoNumber-1)*6) - recCoords(3 + (tomoNumber-1)*6) + 1, ... % NY
               recCoords(5 + (tomoNumber-1)*6), ... % NZ
            -1*recCoords(6 + (tomoNumber-1)*6), ... OX (negative shift X in imod reconstruction command);
        floor((recCoords(4 + (tomoNumber-1)*6) + recCoords(3 + (tomoNumber-1)*6) - 1)/2 - (subTomoMeta.('tiltGeometry').(fileInfo{iTomo,2})(1,21))/2),... % oY -- need the tilt series size
               recCoords(7 + (tomoNumber-1)*6)];  %OZ (negative shift Z in imod reconstruction command -- but rotated during reconstruction so the -1 is implicit);
               
    subTomoMeta.('reconGeometry').(fileInfo{iTomo,2}) = ...
                                                   [recGeom(1:3);recGeom(4:6)];  
%   try
%     recGeom = load(sprintf('./recon/%s_recon.txt',fileInfo{iTomo,2}));
%     subTomoMeta.('reconGeometry').(fileInfo{iTomo,2}) = ...
%                                                    [recGeom(1:3);recGeom(4:6)];
%                                        
%   catch
%     fileInfo{iTomo,2}
%     error('error loading ./recon/%s_recon.txt',fileInfo{iTomo,2});
%   end
  

  % Check to make sure no out of bounds conditions were created in X Y
  rXrY = subTomoMeta.mapBackGeometry.(tiltName).coords;
  
    if rXrY(tomoNumber,2) < -75 
      error(['Out of bounds condition for %s yMin at %f,'...
             'please change recon.txt recon.coords'], ...
              fileInfo{iTomo,2},rXrY(tomoNumber,2))
    elseif rXrY(tomoNumber,2) < 1
      subTomoMeta.mapBackGeometry.(tiltName).coords(tomoNumber,2) = 1;
    end
    yMax = subTomoMeta.('tiltGeometry').(fileInfo{iTomo,2})(1,21);
    if rXrY(tomoNumber,3) > yMax + 75
      error(['Out of bounds condition for %s yMax at %f (max %d),'...
             'please change recon.txt recon.coords'], ...
              fileInfo{iTomo,2},rXrY(tomoNumber,3),yMax)
    elseif rXrY(tomoNumber,3) > yMax
      subTomoMeta.mapBackGeometry.(tiltName).coords(tomoNumber,3) = yMax;
    end
      
      
  subTomoMeta.('ctfGroupSize').(fileInfo{iTomo,2}) = [nCTFgroups,0];
  iX = subTomoMeta.('mapBackGeometry').(fileInfo{iTomo,1}).('coords')(tomoNumber,1);
  subTomoMeta.('ctfGroupSize').(fileInfo{iTomo,2})(2) = floor(iX./...
                                                              nCTFgroups);
% % %   subTomoMeta.('mapExt').(fileInfo{iTomo,2})  = fileInfo{iTomo,3};
% % %   subTomoMeta.('mapPath').(fileInfo{iTomo,1}) = fileInfo{iTomo,2};

  
end

% For now, just assuming all of the maps are in the same place and have the same
% suffix - generalize later.



if nGPUs > nTomogramsTotal
  nGPUs = nTomogramsTotal
end


iterList = cell(nGPUs,1);
for iGPU = 1:nGPUs
  iterList{iGPU} = iGPU:nGPUs:nTomogramsTotal;
  iterList{iGPU}
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
% % % for iGPU = 1:nGPUs

  D = gpuDevice(iGPU);
  
  tomoResults = struct();
  for iTomo = iterList{iGPU}

    if (doImport)
      mapName = fileInfo{iTomo,2};
    else
      mapName = fileInfo{iTomo,3};
    end

    % Load in the template matching geometry for the tomogram, and the model file 
    % which may (or may not) have been edited.
    tomoNumber = subTomoMeta.mapBackGeometry.tomoName.(fileInfo{iTomo,2}).tomoNumber;
    tiltName   = subTomoMeta.mapBackGeometry.tomoName.(fileInfo{iTomo,2}).tiltName;
    coords = subTomoMeta.mapBackGeometry.(tiltName).coords(tomoNumber,1:4);
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
      
          angleSgn
          convention
          direction
          tmpSearchGeom(iAng,17:25) = reshape( ...
                                        BH_defineMatrix( ...
                                          angleSgn .* tmpSearchGeom(iAng,14:16),...
                                          convention,direction),1,9);

        end  
        
        % I'm assuming that the proper scaling was done, let the user know
        fprintf('\nImporting coordinates, assuming to be scaled properly to match full reconstruction size\n');
      else
       % [ ~, binShiftTomo ] = BH_multi_calcBinShift( coords, dupInTheLoop);
       
       % tmpSearchGeom(:,11:13) = tmpSearchGeom(:,11:13) + repmat(binShiftTomo,size(tmpSearchGeom,1),1);
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
      sprintf('model2point convmap/%s.mod %s', mapName, tmpFile)
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


%      if (flgLookForPoints) && any(abs([sx,sy,sz].*dupInTheLoop - subTomoMeta.('reconGeometry').(fileInfo{iTomo,2})(1,1:3)) > 2.*dupInTheLoop)
%        fprintf('convmap/%s_convmap.mrc\n',mapName);
%        error('The binned (bin%d) convmap [%d %d %d] and recon size [%d %d %d] are > %f diff\n',dupInTheLoop,sx,sy,sz,subTomoMeta.('reconGeometry').(fileInfo{iTomo,2})(1,1:3),dupInTheLoop)
%      end
  

    % Assuming that 
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
     
      idxList = positionIDX((overlapMatrix > 1));
      size(overlapMatrix)
      size(positionIDX)
      

      tomoResults.(fileInfo{iTomo,2}) = tmpSearchGeom(ismember(tmpSearchGeom(:,4), idxList),:)
      sum(idxList(:))
      sum(ismember(tmpSearchGeom(:,4), idxList))
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
  for iTomo = iterList{iGPU}
    mapName = fileInfo{iTomo,2};
    tmpGeom = parResults{iGPU}.(mapName);
    
    tmpGeom(:,9:26:26*emc.nPeaks) = repmat(ceil(tmpGeom(:,11)./ ...
                        subTomoMeta.('ctfGroupSize').(mapName)(2)),1,emc.nPeaks);

    % Sort so that CTFs can be left in main mem, and only pulled when needed and only
    % once per round of alignment.
    tmpGeom = sortrows(tmpGeom,9);
  

    for iSubTomo = 1:size(tmpGeom,1)
      tmpGeom(iSubTomo, 4:26:26*emc.nPeaks) = nIDX;
      
      nIDX = nIDX +1; 
    end
    
    
    %tmpGeom(:,9) = ceil(tmpGeom(:,11)./ ...
    %                    subTomoMeta.('ctfGroupSize').(mapName)(2));
                      
    % Using my template matching, there should never be a tomo so close to
    % the edge for this to be problem, but when working with coordinates
    % from relion, I've noticed out of bounds conditions. Check explicitly
    % here.                                           
    tmpGeom( tmpGeom(:,9)> nCTFgroups, 9 ) = nCTFgroups;
    
    subTomoMeta.('cycle000').('geometry').(mapName) = tmpGeom;    
  end
end

fprintf('nSubTomos initial = %d\n', (nIDX-1));

subTomoMeta.('nSubTomoInitial') = nIDX-1;
% subTomoMeta.('nOrientationsPerParticle') = nOrientations;

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

