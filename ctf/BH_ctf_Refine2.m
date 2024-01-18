function [ cccStorage, maxAst, maxAng, astigAngSearch] = BH_ctf_Refine2(PARAMETER_FILE, STACK_PRFX)
% Script to test refinement of ctf estimate by scaling tiles from tilted images
% to change their nominal magnification so that the defocus matches that of the
% mean 

% Load in the tomo and tilt info
pBH = BH_parseParameterFile(PARAMETER_FILE);
try
  load(sprintf('%s.mat', pBH.('subTomoMeta')), 'subTomoMeta');
  mapBackIter = subTomoMeta.currentTomoCPR; clear subTomoMeta
catch
  mapBackIter = 0;
end

try
  testNoRefine = pBH.('force_no_defocus_stretch');
catch
  testNoRefine = false;
end
if (testNoRefine)
  fprintf('\nWarning force_no_defocus_stretch is only for testing!\n\n');
end

% As long as the material coming into view at high tilt is at the same
% plane and does not have wildly different image stats, using it could
% improve the thon rings on tilted data. Zero will produce the "normal"
% process, 1 will use the same area as the min tilt, 
try 
  fraction_of_extra_tilt_data = pBH.('fraction_of_extra_tilt_data')
catch
  fraction_of_extra_tilt_data = 0.25
end
% set the search ranges - should change ctf_est to save the parameters used so
% this can be loaded automatically.

nWorkers = min(pBH.('nCpuCores'),7*pBH.('nGPUs'));
nWorkers = BH_multi_parallelWorkers(nWorkers)


gpuIDX = BH_multi_checkGPU(-1);
gDev = gpuDevice(gpuIDX);


reScaleRealSpace = 0;


PRJ_STACK = {sprintf('aliStacks/%s_ali%d.fixed',STACK_PRFX,mapBackIter+1)};
[pathName,fileName,extension] = fileparts(PRJ_STACK{1});
if isempty(pathName)
  pathName = '.';
end

PIXEL_SIZE = pBH.('PIXEL_SIZE');


Cs = pBH.('Cs');
VOLTAGE = pBH.('VOLTAGE');
AMPCONT = pBH.('AMPCONT');

ctfParams = [PIXEL_SIZE*10^10,VOLTAGE./1000,Cs.*1000,AMPCONT];


% Sanity check
if (PIXEL_SIZE > 20e-10 || PIXEL_SIZE < 0)
  error('pixel size should be [0,20e-10]');
elseif (Cs > 5*10^-3 || Cs < 0)
  error('Cs should be[1e-3,10e-3]');
elseif(VOLTAGE > 1000e3 || VOLTAGE < 20e3)
  error ('VOLTAGE should be [20e3,1000e3]');
elseif (AMPCONT < 0.025 || AMPCONT > 0.25)
  error('AMPCONT should be [0.025,0.25]');
else
  WAVELENGTH = 10^-12*1226.39/sqrt(VOLTAGE + 0.97845*10^-6*VOLTAGE^2) ;
end

if (Cs == 0)
  fprintf('You set Cs to zero, over-riding to 5 micron\n');
  Cs = 5e-6;
end

% Assuming that the first CTF zero is always less than this value 
FIXED_FIRSTZERO =  PIXEL_SIZE / 40*10^-10 ;

% Size to padTile to should be even, large, and preferably a power of 2
try
  paddedSize = pBH.('paddedSize');
catch
  paddedSize = 768;
end

% Tile size & overlap
tileOverlap = 4;
tileSize = floor(680e-10 / PIXEL_SIZE);
tileSize = tileSize + mod(tileSize,2);
fprintf('Using a tile size of %d',tileSize);
overlap = floor(tileSize ./ tileOverlap);


inc = (0.5 - FIXED_FIRSTZERO) / (paddedSize/2);
freqVector = [inc+FIXED_FIRSTZERO:inc:0.5 ];


tlt = {sprintf('fixedStacks/ctf/%s_ali%d_ctf.tlt',STACK_PRFX,mapBackIter+1)};
%PRJ_OUT = {fileName};

nStacks = length(tlt);
stacksFound = [];
INPUT_CELL = cell(nStacks,6);


for iStack = 1:nStacks
  if exist(tlt{iStack}, 'file') &&  exist(PRJ_STACK{iStack}, 'file')
    INPUT_CELL{iStack,1} = load(tlt{iStack});
    INPUT_CELL{iStack,2} = PRJ_STACK{iStack};
    [pathName,fileName,extension] = fileparts(PRJ_STACK{iStack});
    if isempty(pathName)
      pathName = '.';
    end
    INPUT_CELL{iStack,3} = pathName;
    INPUT_CELL{iStack,4} = fileName;
    INPUT_CELL{iStack,5} = extension;
    %INPUT_CELL{iStack,6} = PRJ_OUT;

  else
    fprintf('ignoring %s, because the file is not found.\n', tlt{iStack})
    
  end
end  

for iStack = 1%stacksFound

  STACK = single(getVolume(MRCImage(INPUT_CELL{iStack,2})));
  % The pixel size should be previously set correctly, but if it is not, then we
  % must maintain whatever is there in case beads are to be erased. The model
  % used for this process depends on the pixel size in the header when it was
  % created in IMod alignment.
  [~,iPixelHeader] = system(sprintf('header -pixel %s',INPUT_CELL{iStack,2}));
  iPixelHeader = EMC_str2double(iPixelHeader);
  [d1,d2,d3] = size(STACK)
    
    
    
  TLT = INPUT_CELL{iStack,1};
  pathName = INPUT_CELL{iStack,3};
  fileName = INPUT_CELL{iStack,4};
  extension = INPUT_CELL{iStack,5};

      
  SIZEOUT = [d1,d2];

  [radialForCTF,phi,~,~,~,~] = ...
    BH_multi_gridCoordinates([paddedSize,paddedSize,1],'Cylindrical','GPU',{'none'},1,1,0);


  radialForCTF = {radialForCTF./PIXEL_SIZE,1,phi}  ; 
  
  clear phi                                                        

  clear sumVector radialAvg
  sumVector(length(freqVector)) = gpuArray(double(0));
  radialAvg(length(freqVector)) = gpuArray(double(0));

  tic

  psTile = zeros([paddedSize,paddedSize,d3],'single');

  flgReplaceStack = 0;
  for iPrj = 1:d3
    iProjection = gpuArray(STACK(:,:,TLT(iPrj,1)));
    iProjection = iProjection - ...
                        BH_movingAverage(iProjection,[tileSize,tileSize]);
    iProjection = iProjection ./ ...
                        BH_movingRMS(iProjection,[tileSize,tileSize]);
    % Taking a cue from Alexis
    maxPixelSizeWanted = 2.0e-10;
    if TLT(iPrj,16) < maxPixelSizeWanted
        %fprintf(ftmp,'Resampling pixel size\n');
        %  Resample to 2Ang/pix
        padSq = BH_multi_padVal(size(iProjection),max(size(iProjection)).*[1,1]);

        iProjection = BH_padZeros3d(iProjection,padSq(1,:),padSq(2,:),'GPU','singleTaper');
        sizeIN = size(iProjection,1);
        % Replace with BH_fftShift if this works
        iProjection = fftshift(fftn(iProjection));
        trimVal = BH_multi_padVal(size(iProjection), floor(size(iProjection).*(TLT(iPrj,16)./maxPixelSizeWanted)));
        iProjection = real(ifftn(ifftshift(BH_padZeros3d(iProjection,trimVal(1,:),trimVal(2,:),'GPU','single'))));
        sizeOUT = size(iProjection);
        if iPrj == 1
          flgReplaceStack = 1;
          newSTACK = zeros([sizeOUT,d3],'single');
        end
        newSTACK(:,:,TLT(iPrj,1)) = gather(iProjection);
        clear iProjection
        % Actual new pixel size
        pixelSize = sizeIN./sizeOUT(1).*TLT(iPrj,16);
        
        % Update the CTF params with the new pixelSize
        ctfParams(1) = pixelSize.*10^10;
        
        %fprintf(ftmp,'%d %d %d %d %d %d\n',trimVal);
        %fprintf(ftmp,'pixelOld %3.3e, pixelNew %3.3e\n',TLT(iPrj,16),pixelSize);

        
      else
        pixelSize = TLT(iPrj,16);
      end
    end  % iPrj 1:d3

    if ( flgReplaceStack )
      STACK = newSTACK ; clear newSTACK;
    end
    [d1,d2,d3] = size(STACK);
    [Xnew, Ynew, ~, x1,y1, ~] = BH_multi_gridCoordinates([tileSize,d2], ...
    'Cartesian','GPU', ...
    {'none'},0,1,0);

    [X, Y, ~,~,~, ~] = BH_multi_gridCoordinates([tileSize,tileSize], ...
    'Cartesian','GPU', ...
    {'none'},0,1,0);

    coordShift = (-1).^(X+Y);
    clear X Y

    try 
      ppool = EMC_parpool(nWorkers);
    catch 
        delete(gcp('nocreate'));
        ppool = EMC_parpool(nWorkers);
    end

    for iPrj = 1:d3

      pFuture(iPrj) = parfeval(ppool,@runAvgTiles,2, TLT, paddedSize, tileSize, ...
                                                  d1,d2, iPrj, overlap, ...
                                                STACK(:,:,TLT(iPrj,1)), ...
                                                1, ...
                                                1, ...
                                                x1, y1, Xnew, Ynew,coordShift,  ...
                                                reScaleRealSpace,pixelSize,fraction_of_extra_tilt_data,testNoRefine);   
                                            

    end
  
    for iWorker = 1:d3
      fprintf('Refining defocus on prj %d/ %d\n',iWorker,d3);
      [iPrj, ctfCorr,pixelSize] = fetchNext(pFuture);

      psTile(:,:,TLT(iPrj,1)) = ctfCorr;
    end

    pixelSize = pixelSize*10^10;
    SAVE_IMG(MRCImage(gather(psTile)),sprintf('fixedStacks/ctf/%s-PS.mrc',fileName),pixelSize);
    bpLog = fftshift(BH_bandpass3d([size(psTile(:,:,1)),1],0,0,2.2.*pixelSize,'GPU',pixelSize));
    bpLog = bpLog > 0.99;
    bp  = fftshift(BH_bandpass3d([size(psTile(:,:,1)),1],0.25,20,2.*pixelSize,'GPU',pixelSize));
    bp2 = fftshift(BH_bandpass3d([size(psTile(:,:,1)),1],1e-6,400,2.*pixelSize,'GPU',pixelSize));
  
    for iPrj = 1:d3
      iTile = gpuArray(psTile(:,:,iPrj));
      iTile = iTile.*bp.*bp2;
      iTile(~bpLog) = mean(iTile(bpLog));
      psTile(:,:,iPrj) = gather(iTile);
      
    end

    SAVE_IMG(MRCImage(gather(psTile)),sprintf('fixedStacks/ctf/%s-PS2.mrc',fileName),pixelSize);
    delete(ppool);
    delete(gcp('nocreate'))
  
    % exit an fit the PS using CTFFIND4
      BH_runCtfFind(sprintf('fixedStacks/ctf/%s-PS2.mrc',fileName), ...
                    sprintf('%s_ctf.tlt',fileName), ctfParams,TLT)

 
  end
end


function [psTile,pixelSize] = runAvgTiles(TLT, paddedSize, tileSize, d1,d2, iPrj, overlap, ...
                          iProjection, evalMask, ...
                          ddZ, x1, y1, Xnew, Ynew, coordShift,  ...
                          reScaleRealSpace,pixelSize,fraction_of_extra_tilt_data,testNoRefine)

    DFo = TLT(iPrj,15);

    padTileOver = 256;
    tmpTile = zeros(paddedSize.*[1,1]+2*padTileOver,'single','gpuArray');
    
    tiltOrigin = ceil((size(iProjection,1)+1)./2);
    
    oXprj = ceil((size(iProjection,1)+1)./2);
    % Don't worry about extending the edges for thickness
    half_width = (size(iProjection,1)/2);

    maxEval = (fraction_of_extra_tilt_data + ...
               cosd(TLT(iPrj,4)).*(1-fraction_of_extra_tilt_data)) .* half_width;

    iEvalMask = floor(oXprj-maxEval):ceil(oXprj+maxEval);
    
    % Since I'm enforcing Y-tilt axis, then this could be dramatically sped up
    % by resampling strips along the sampling

    for iOuter = 1+tileSize/2:overlap:d1-tileSize/2
      randSize = randi(floor(overlap/2),1);
      if (randi(2,1) == 2)
        randSize = -1*randSize;
      end
      i = iOuter + randSize;
      if (i < tileSize/2 || i > d1-tileSize/2)
        continue;
      end
      
      % Slightly randomize the step size to avoid a Moire like effect that
      % presents particulary strongly with a continuous carbon layer.
      
      iDeltaZ = (i - tiltOrigin)*pixelSize*-1.*tand(TLT(iPrj,4));
      if any(ismember(i-tileSize/2+1:i+tileSize/2,iEvalMask)) %evalMask(i,paddedSize/2+1)

        mag =  (1+iDeltaZ./DFo).^0.5;
          
        estSize = 2048;
        ctf1 = BH_ctfCalc(pixelSize,TLT(iPrj,17),TLT(iPrj,18),DFo,estSize,TLT(iPrj,19),-1,1);
        ctf2 = BH_ctfCalc(pixelSize,TLT(iPrj,17),TLT(iPrj,18),iDeltaZ+DFo,estSize,TLT(iPrj,19),-1,1);
        ctf1 = ctf1(1:estSize/2);
        ctf2 = ctf2(1:estSize/2);
        firstZero = find(ctf1 > 0, 1, 'first');
        % secondZero= find(ctf1(firstZero:end) < 0 , 1, 'first') + firstZero - 1;

        %fprintf(ftmp,'firstZero %d %2.2f\n',firstZero,estSize/firstZero*pixelSize);
        % This range will depend on the size of the field of view. For now, setting manually for Florian's HIV
        % data, but will derive a formula to make sure the search is appropriate. Here we expect at most ~ 300 nm
        % deltaZ, the strongest difference is at the lowest defocus which is ~ 1500 nm, which gives an estimated mag
        % ~ 1.095
        defRange = mag-.1:.001:mag+.1;
        nDef = length(defRange);
        scoreDef = zeros(nDef,1,'gpuArray'); 
        for iDef = 1:nDef
          ci = interpn([1:estSize/2]',ctf2(1:estSize/2),[1:estSize/2]'./defRange(iDef),'linear',0);
          % Larger scalings will have zeros rather than extroplation, so%
  %            % don't let this influence the score.
  %           lastZero = find(abs(ci) > 0 , 1, 'last');
          %fprintf(ftmp,'%d %d %d %d',size(ci),size(ctf1)); 
          scoreDef(iDef) = sum(ci(firstZero:end).*ctf1(firstZero:end))./sqrt(sum(ci(firstZero:end).^2).*sum(ctf1(firstZero:end).^2));
          
        end
        [~,maxCoord] = max(scoreDef);
        mag = defRange(maxCoord);

        defRange = mag-.01:.0001:mag+.01;
        nDef = length(defRange);
        scoreDef = zeros(nDef,1,'gpuArray');
        for iDef = 1:nDef
          ci = interpn([1:estSize/2]',ctf2(1:estSize/2),[1:estSize/2]'./defRange(iDef),'linear',0);
          scoreDef(iDef) = sum(ci(firstZero:end).*ctf1(firstZero:end))./sqrt(sum(ci(firstZero:end).^2).*sum(ctf1(firstZero:end).^2));
        end
        [~,maxCoord] = max(scoreDef);
        mag = defRange(maxCoord);
         
        if (testNoRefine)
          mag = 1;
        end

        scaledStrip = iProjection(i-tileSize/2+1:i+tileSize/2,:); 

        for j = 1+tileSize/2:overlap:d2-tileSize/2    
        
          iTile = scaledStrip(:,j-tileSize/2+1:j+tileSize/2);%.*coordShift;
          if (reScaleRealSpace)
            scaledSize = paddedSize;
          else
            % Slightly randomize scaling
            if (randi(2,1) == 2)
              scaledSize = ceil(paddedSize .* mag) + randi(2,1) -1;
            else
              scaledSize = floor(paddedSize .* mag)+ randi(2,1) -1;
            end
            %scaledSize = floor(paddedSize ./ mag);
          end
                  
          [oX,oY] = size(tmpTile);
          oX = ceil((oX+1)./2);
          oY = ceil((oY+1)./2);

          iPadVal = BH_multi_padVal(size(iTile),[scaledSize,scaledSize]);


          oupSize = [floor(scaledSize./2),ceil(scaledSize./2); ...
                      floor(scaledSize./2),ceil(scaledSize./2)];
    
          % Get rid of th fftshift
          tmpTile(oX-oupSize(1,1):oX+oupSize(1,2)-1, ...
                  oY-oupSize(2,1):oY+oupSize(2,2)-1) = ...         
          tmpTile(oX-oupSize(1,1):oX+oupSize(1,2)-1, ...
                  oY-oupSize(2,1):oY+oupSize(2,2)-1) + ...
                  fftshift(abs(fftn(BH_padZeros3d(iTile,iPadVal(1,:),iPadVal(2,:), ...
                                          'GPU','singleTaper', mean(iTile(:))))));
                
        end % loop over j
      end % if over eval mask
    end % over tiles
  
    
    psTile = gather(BH_padZeros3d(tmpTile, [-1,-1].* ...
                                   padTileOver,[-1,-1].*padTileOver,...
                                   'GPU','single'));
    clear tmpTile iProjection ddZ evalMask Xnew Ynew x1 y1 
end
  
function [minRes] = calcMinResolution(TLT, radialForCTF,Cs,WAVELENGTH,AMPCONT)

    meanDef = mean(TLT(:,15));
    meanAst = mean(TLT(:,12));
    meanAng = mean(TLT(:,13));

    df1 = meanDef - meanAst;
    df2 = meanDef + meanAst;

    [ Hqz ] = BH_ctfCalc(radialForCTF,Cs,WAVELENGTH, ...
                         [df1,df2,meanAng],size(radialForCTF{1}), ...
                          AMPCONT,-1.0); 
    rV =  Hqz(1+size(Hqz, 1)/2,1+size(Hqz, 1)/2:end);
    freqVector  = radialForCTF{1}(1+size(Hqz, 1)/2,1+size(Hqz, 1)/2:end);
    firstZero = find(rV > 0, 1,'first');

    minRes = 1/freqVector(firstZero) * 10^10;

end
