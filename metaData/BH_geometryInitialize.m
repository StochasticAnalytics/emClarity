function [ ] = BH_geometryInitialize( PARAMETER_FILE, varargin )
%Take a set of coordinates and optionally angles and set up geometry.
%   For now just BH, but also to work with Relion set up.

% Option for testing euler angles will remove.

convertEulers = 0;
direction = '';
% if nargin == 2

%   switch SPI
%     case 'fwd'
%       direction = 'forward'
%       convertEulers = 1
%     case 'nFwd'
%       direction = 'forward'
%       convertEulers = -1
%     case 'inv'
%       direction = 'inv'
%       convertEulers = 1
%     case 'nInv'
%       direction = 'inv'
%       convertEulers = -1
%     otherwise
%       error('SPI fwd nFwd inv nInv NOT %s\n',SPI)
%   end
% end
    
% Initialize the subTomoMeta structure
subTomoMeta = struct();

pBH = BH_parseParameterFile(PARAMETER_FILE);
if nargin == 2
  mapBackIter = str2double(varargin{1});
else
  mapBackIter = 0;
end
nGPUs = pBH.('nGPUs');
% This will have to do until a better approach based on PSF of positions
% in projection space linked together can be used to define groups that
% don't have co-mingled resolution. When splitOnTOmos, always run in serial, so
% that we have the best chance of distributing the defocus variateion/tomoqualtiy.
splitOnTomos = pBH.('fscGoldSplitOnTomos');
if (splitOnTomos)
  nGPUs = 1;
  fprintf('override nGPUs to just 1 for initial step to evenly split crowded tomos because fscGoldSplitOnTomos is true')
end
nOrientations=1;%nOrientations = pBH.('pseudoMLnumber');
nCTFgroups = 9;
try
  nPeaks = pBH.('nPeaks');
catch
  nPeaks = 1;
end

% Resolution lower than this is not gold standard, and will also be mixed
% in the references to keep orientations from diverging. Should be > 2.25 x
% your expected resolution. 1.5x to be totally uncorrelated, and another
% 1.5x because aligning at ~40A in template matching will (with good data)
% often lead to and FSC of ~ 26A initially. 

% 40 may be too conservative, especially for low defocus tomos
try
  lowResCut = pBH.('lowResCut');
catch
  lowResCut = 40;
end

maxGoldStandard = lowResCut;%pBH.('Tmp_bandpassFilter')(3);
dupSampling = pBH.('Tmp_samplingRate')
dupRadius = 2;
dupTolerance = (2.*dupRadius)+1;
dupMask = zeros(dupTolerance, dupTolerance, dupTolerance, 'single');
dupMask = dupMask + 1;

% If we are in the working directory, following template matching, there should
% be the director convmap, holding convolution maps, model files etc.

checkDir = dir('convmap');
if length(checkDir) < 1 || ~checkDir(1).isdir
  error('Did not find directory named <convmap>')
end

% Each tomogram will have its own file suffix .path which holds the basename,
% path, file extension, tilt file. Get these and concat into a cell.

getPath = dir('convmap/*.mod');

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


[ recGeom, tiltName, nTomosPossible] = BH_multi_recGeom( sprintf('recon/%s',getCoords(iStack).name) )
  % Initialize
 
                                                       
  subTomoMeta.('mapBackGeometry').(tiltName).('coords') = zeros(nTomosPossible,6);
  iPath = dir(sprintf('convmap/%s_*.mod',tiltName))
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

    modName = strsplit(iPath(iSt).name,'_');
    tomoNumber = str2double(modName{end-1});

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
  tiltName = strjoin(modName(1:end-2),'_');
  tomoNumber = str2double(modName{end-1})


  fileInfo{iTomo,1} = tiltName;
  fileInfo{iTomo,2} = sprintf('%s_%d',tiltName,tomoNumber);
  fileInfo{iTomo,3} = sprintf('%s_%d_bin%d',tiltName,tomoNumber,dupSampling);
  fileInfo{iTomo,4} = sprintf('fixedStacks/ctf/%s_ali%d_ctf.tlt',tiltName,mapBackIter+1);
  
  try
    recGeom = load(sprintf('./recon/%s_recon.txt',fileInfo{iTomo,2}));
    subTomoMeta.('reconGeometry').(fileInfo{iTomo,2}) = ...
                                                   [recGeom(1:3);recGeom(4:6)];
                                       
  catch
    fileInfo{iTomo,2}
    error('error loading ./recon/%s_recon.txt',fileInfo{iTomo,2});
  end
  
  subTomoMeta.('tiltGeometry').(fileInfo{iTomo,2}) = load(fileInfo{iTomo,4});
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






iterList = cell(nGPUs,1);
for iGPU = 1:nGPUs
  iterList{iGPU} = iGPU:nGPUs:nTomogramsTotal;
  iterList{iGPU}
end

try
  parpool(nGPUs)
catch
  delete(gcp)
  parpool(nGPUs)
end

parResults = cell(nGPUs,1);

dupInTheLoop = dupSampling
parfor iGPU = 1:nGPUs

  D = gpuDevice(iGPU);
  
  tomoResults = struct();
  for iTomo = iterList{iGPU}

    mapName = fileInfo{iTomo,3};
% % %     mapPath = fileInfo{iTomo,2};
% % %     mapExt  = fileInfo{iTomo,3};
    % Load in the template matching geometry for the tomogram, and the model file 
    % which may (or may not) have been edited.
    try
      tmpSearchGeom = importdata(sprintf('convmap/%s.csv',mapName));
    catch
      error('Did not find the file convmap/%s.csv\n\nMake sure the binning in your parameter file is correct.',mapName)
    end
    tmpSearchGeom(1,:)
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
      sprintf('model2point -float convmap/%s.mod %s', mapName, tmpFile)
      system(sprintf('model2point -float convmap/%s.mod %s', mapName, tmpFile))
      modGeom = load(tmpFile);
      system(sprintf('rm %s', tmpFile));
        flgLookForPoints = 1;

            % leave IDX in main memory because it is just for reference.
      sx = floor(iHeader.nX);
      sy = floor(iHeader.nY);
      sz = floor(iHeader.nZ);
      
    catch
      tiltName = sprintf('recon/%s.rec',mapName);
      sx= recGeom(1,1) ; sy = recGeom(1,2); sz = recGeom(1,3);
      fprintf('did not find convmap, assuming synthetic data, with no deleted model points\n');
      dupSampling = 1;
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
          tmpSearchGeom(iSubTomo,26:26:26*nPeaks) = -9999;
        else
        positionMatrix(subTomoOrigin(1),subTomoOrigin(2),subTomoOrigin(3)) = 1;
        positionIDX(subTomoOrigin(1),subTomoOrigin(2),subTomoOrigin(3)) = ...
                                                       tmpSearchGeom(iSubTomo, 4);

        end

      end % loop building position matrix

      % Add in the points from the modified list, which will increase the value at
      % any retained positions to 2.

      for iSubTomo = 1:size(modGeom,1)
        
        subTomoOrigin = fix(modGeom(iSubTomo,:))
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
  
  end
  parResults{iGPU} = tomoResults;
    
    
end

% By my convention, every particle should have a unique id.

nIDX = 1;
for iGPU = 1:nGPUs
  for iTomo = iterList{iGPU}
    mapName = fileInfo{iTomo,2};
    tmpGeom = parResults{iGPU}.(mapName);
    
    tmpGeom(:,9:26:26*nPeaks) = repmat(ceil(tmpGeom(:,11)./ ...
                        subTomoMeta.('ctfGroupSize').(mapName)(2)),1,nPeaks);

    % Sort so that CTFs can be left in main mem, and only pulled when needed and only
    % once per round of alignment.
    tmpGeom = sortrows(tmpGeom,9);
  
    for iSubTomo = 1:size(tmpGeom,1)
      tmpGeom(iSubTomo, 4:26:26*nPeaks) = nIDX;
      nIDX = nIDX +1; 
    end
    

    %tmpGeom(:,9) = ceil(tmpGeom(:,11)./ ...
    %                    subTomoMeta.('ctfGroupSize').(mapName)(2));
                      
% % %     % Using my template matching, there should never be a tomo so close to
% % %     % the edge for this to be problem, but when working with coordinates
% % %     % from relion, I've noticed out of bounds conditions. Check explicitly
% % %     % here.                                           
% % %     tmpGeom( tmpGeom(:,9) == nCTFgroups + 1, 9 ) = nCTFgroups;
    
    subTomoMeta.('cycle000').('geometry').(mapName) = tmpGeom;    
  end
end

fprintf('nSubTomos initial = %d\n', (nIDX-1));

subTomoMeta.('nSubTomoInitial') = nIDX-1;
% subTomoMeta.('nOrientationsPerParticle') = nOrientations;

preFscSplit = gather(subTomoMeta);

% Randomly divide the data into half sets.
[ subTomoMeta ] = BH_fscSplit( preFscSplit, splitOnTomos, nPeaks);
subTomoMeta.('currentCycle') = 0;
subTomoMeta.('currentTomoCPR') = mapBackIter;
subTomoMeta.('currentResForDefocusError') = lowResCut;
subTomoMeta.('maxGoldStandard') = maxGoldStandard;
save(sprintf('%s.mat', pBH.('subTomoMeta')), 'subTomoMeta');


end % end of initialize function

