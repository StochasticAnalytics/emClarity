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

flgAstigmatism= 1;
flgGroupProjections =1;% pBH.('flgGroupProjections');
calcAvg = 1;%pBH.('flgCalcAvg');
outputForCTFFIND = 1;
reScaleRealSpace = 0;
normFactor=0;%pBH.('normalizationFactor');


PRJ_STACK = {sprintf('aliStacks/%s_ali%d.fixed',STACK_PRFX,mapBackIter+1)};
[pathName,fileName,extension] = fileparts(PRJ_STACK{1})
if isempty(pathName)
  pathName = '.'
end

PIXEL_SIZE = pBH.('PIXEL_SIZE');


Cs = pBH.('Cs');
VOLTAGE = pBH.('VOLTAGE');
AMPCONT = pBH.('AMPCONT')

ctfParams = [PIXEL_SIZE*10^10,VOLTAGE./1000,Cs.*1000,AMPCONT]



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

if Cs == 0
  fprintf('You set Cs to zero, over-riding to 5 micron\n');
  Cs = 5e-6;
end

% Assuming that the first CTF zero is always less than this value 
FIXED_FIRSTZERO =  PIXEL_SIZE / 40*10^-10 ;

% highCutoff = PIXEL_SIZE/pBH.('defCutOff');
highCutoff = 1/pBH.('defCutOff');


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
overlap = floor(tileSize ./ tileOverlap)



% Starting at +/- 750nm
deltaZTolerance = 750e-9 / PIXEL_SIZE;
% Use to check for proper gradient.
zShift = 0;

inc = (0.5 - FIXED_FIRSTZERO) / (paddedSize/2);
freqVector = [inc+FIXED_FIRSTZERO:inc:0.5 ];


tlt = {sprintf('fixedStacks/ctf/%s_ali%d_ctf.tlt',STACK_PRFX,mapBackIter+1)};
%PRJ_OUT = {fileName};

nStacks = length(tlt);
 stacksFound = []
INPUT_CELL = cell(nStacks,6);


for i = 1:nStacks
  if exist(tlt{i}, 'file') &&  exist(PRJ_STACK{i}, 'file')
    INPUT_CELL{i,1} = load(tlt{i});
    INPUT_CELL{i,2} = PRJ_STACK{i};
    [pathName,fileName,extension] = fileparts(PRJ_STACK{i});
    if isempty(pathName)
      pathName = '.';
    end
    INPUT_CELL{i,3} = pathName;
    INPUT_CELL{i,4} = fileName;
    INPUT_CELL{i,5} = extension;
    %INPUT_CELL{i,6} = PRJ_OUT;
 
  else
    fprintf('ignoring %s, because the file is not found.\n', tlt{i})
    
  end
end  

for iStack = 1%stacksFound
  iStack
  
  STACK = single(getVolume(MRCImage(INPUT_CELL{iStack,2})));
  % The pixel size should be previously set correctly, but if it is not, then we
  % must maintain whatever is there in case beads are to be erased. The model
  % used for this process depends on the pixel size in the header when it was
  % created in IMod alignment.
  [~,iPixelHeader] = system(sprintf('header -pixel %s',INPUT_CELL{iStack,2}));
  iPixelHeader = str2num(iPixelHeader);
  [d1,d2,d3] = size(STACK)
  
  
  
  TLT = INPUT_CELL{iStack,1};
  pathName = INPUT_CELL{iStack,3}
  fileName = INPUT_CELL{iStack,4}
  extension = INPUT_CELL{iStack,5} 

  
      
  SIZEOUT = [d1,d2];

 



  [radialForCTF,phi,~,~,~,~] = ...
     BH_multi_gridCoordinates([paddedSize,paddedSize,1],'Cylindrical','GPU',{'none'},1,1,0);


  radialForCTF = {radialForCTF./PIXEL_SIZE,1,phi}  ; 
  
  clear phi

 
if (calcAvg)                                                            
%     [exposureFilter] = BH_exposureFilter(paddedSize.*[1,1], TLT,'cpu',1, 1);


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
    end   
  
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
    % with > 250,000 tiles, it doesn't make sense to call resample2d for just a
    % simple rescaling.



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
    
    for i = 1:d3
      fprintf('Refining defocus on prj %d/ %d\n',i,d3);
      [iPrj, ctfCorr,pixelSize] = fetchNext(pFuture);

      psTile(:,:,TLT(iPrj,1)) = ctfCorr;


    end
    %%pixelSize = PIXEL_SIZE*10^10;
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
else

    psTile = gpuArray(getVolume(MRCImage(sprintf('fixedStacks/ctf/%s-PS.mrc',fileName))));
  end
  
  if (calcAvg)
    delete(gcp('nocreate'))
  end  
 
 %%%%%%%%%%%%%%%%%%%%%%%%%
 if ( outputForCTFFIND )
   % exit an fit the PS using CTFFIND4
    BH_runCtfFind(sprintf('fixedStacks/ctf/%s-PS2.mrc',fileName), ...
                  sprintf('%s_ctf.tlt',fileName), ctfParams,TLT)
   return
 end
 %%%%%%%%%%%%%%%%%%%%%%%%%
  defInc = cell(3,1); astInc = cell(3,1);
  defRange = cell(3,1); astRange = cell(3,1);
  defSearch = cell(3,d3); astSearch = cell(3,1);
  
  % Find a close value for symmetric defocus. Do this on a rotationally averaged
  % image so that the value is centered between the astigmatic extremes rather
  % than potentially sitting at one end or the other.
  % Search around this value to get a ballpark on astigmatism.
  % Use this astigmatic value to refine symmetric results, and this to refine
  % astigmatic results.
  % Lather rinse and repeat 1x.
  
  
  defInc{1} = 25*10^-9;
  defRange{1} = 1000 *10^-9;
  
  defInc{2} = 10*10^-9;
  defRange{2} = 100*10^-9;

  defInc{3} = 5 *10^-9;
  defRange{3} = 50*10^-9;

  maxAstig = 200*10^-9;
  astigStep =10*10^-9;
  coarseAngStep = (pi/180)*10;
  
  astigDefSearch{1} = 0:astigStep:maxAstig;
%   astigAngSearch{1} = -pi/2:coarseAngStep:pi/2;
  astigAngSearch{1} = 0:coarseAngStep:pi;
  
  astigDefSearch{2} = -5*astigStep:astigStep/2:astigStep*5;
  astigAngSearch{2} = -2*coarseAngStep:coarseAngStep/5:coarseAngStep*2;
  
  astigDefSearch{3} = -3*astigStep:astigStep/4:astigStep*3;
  astigAngSearch{3} = -coarseAngStep/2:coarseAngStep/20:coarseAngStep/2;

  astigCCC = cell(3,d3);
  for iSearch = 1:3
    for iPrj = 1:d3
    astigCCC{iSearch}{iPrj} = zeros(length(astigDefSearch{iSearch})* ...
                     length(astigAngSearch{iSearch}),3, 'gpuArray'); 
    end
  end
  
  rotationalAvg = 1;
  if (rotationalAvg)
    rotBgSubPS = zeros(size(psTile),'single','gpuArray');
  end
  
  bgSubPS = zeros(size(psTile),'single','gpuArray');
  minRes = calcMinResolution(TLT, radialForCTF, Cs,WAVELENGTH,AMPCONT);
  fprintf('\nMin Resolution fit is %3.3f Angstrom.\n',minRes);
  nBgPix = floor(paddedSize.*PIXEL_SIZE.*10^10*sqrt(2)/minRes);
  nBgPix = nBgPix + mod(nBgPix,2)
  if (normFactor)
    nRMSpix = floor(nBgPix/normFactor) + mod(nBgPix/normFactor,2);
  end
  cccResults = zeros(d3,1);
% 

  [rot1, rot2, ~, r1,r2, ~] = BH_multi_gridCoordinates(paddedSize.*[1,1], ...
  'Cartesian','GPU', ...
  {'none'},0,1,0);

  if ( flgGroupProjections )
    % Create a temporary copy of the tiles scaled by the sin of the tilt angle
    % which results in stronger averaging at high tilts
    scaledTile = zeros(size(psTile),'single');

    % Since they will be averaged, first center and scale the total per 
    % prj intensities
    for iPrj = 1:d3
      psTile(:,:,iPrj) = psTile(:,:,iPrj) - mean(mean(psTile(:,:,iPrj)));
      psTile(:,:,iPrj) = psTile(:,:,iPrj) ./ rms(rms( psTile(:,:,iPrj)));
    end

    for iPrj = 1:d3
      scaledTile(:,:,TLT(iPrj,1)) = scaledTile(:,:,TLT(iPrj,1)) .* ...
                                                    abs(sind(TLT(iPrj,4)))+0.05;
    end

  end

  for iPrj = 1:d3

    if (flgGroupProjections)  

      % Take the full value at the projection of interest, sin(tiltangle) for the
      if iPrj == 1
        gTMP = psTile(:,:,1)   + scaledTile(:,:,2) + scaledTile(:,:,3);
      elseif iPrj == d3
        gTMP = psTile(:,:,d3)  + scaledTile(:,:,d3-1) + scaledTile(:,:,d3-2)
      else
        gTMP = psTile(:,:,iPrj)+ scaledTile(:,:,iPrj-1) + scaledTile(:,:,iPrj+1);
      end
    else
      
      gTMP = psTile(:,:,iPrj);
    end

%      % Take the full value at the projection of interest, sin(tiltangle) for the
%      if iPrj == 1
%      gTMP=(psTile(:,:,1) + ...
%           psTile(:,:,2).*0.55 +...
%           psTile(:,:,3).*0.25)./1.9;
%      elseif (iPrj > 1 && iPrj <  7) || (iPrj > d3-6 && iPrj < d3)
%     gTMP = (psTile(:,:,iPrj-1).*0.4+...
%             psTile(:,:,iPrj)+...
%            psTile(:,:,iPrj+1).*0.4)./1.8;
%      elseif (iPrj >= 7 && iPrj <  18) 
%      gTMP = (psTile(:,:,iPrj) + ...
%                   psTile(:,:,iPrj+1).*0.3)./1.3; 
%     elseif (iPrj > d3-17 && iPrj <= d3-6)
%      gTMP = (psTile(:,:,iPrj-1).*0.3+...
%                   psTile(:,:,iPrj))./1.3; 
%      elseif iPrj == d3
%      gTMP = (psTile(:,:,d3-2).*0.25+...
%                   psTile(:,:,d3-1).*0.55+...
%                   psTile(:,:,d3))./1.8;
%      else
%      gTMP = psTile(:,:,iPrj);
%      end
%    else
      
%      gTMP = psTile(:,:,iPrj);
%    end

    
   gTMP = gTMP - BH_movingAverage(gTMP,[nBgPix,nBgPix]);
   if (normFactor)
     bgSubPS(:,:,iPrj) = gTMP ./ BH_movingRMS(gTMP,floor([nRMSpix,nRMSpix]));
   else
     bgSubPS(:,:,iPrj) = gTMP;
   end
   clear gTMP
   
   if (rotationalAvg)
     
     rotTMP = bgSubPS(:,:,iPrj);
     for i = 0.5:0.5:360
        R = BH_defineMatrix([i,0,0],'Bah','forward');
        ROT1 = R(1).*rot1 + R(4).*rot2;
        ROT2 = R(2).*rot1 + R(5).*rot2;
        rotTMP = rotTMP + interpn(r1,r2,bgSubPS(:,:,iPrj),ROT1,ROT2,'linear',0);
     end   
      rotBgSubPS(:,:,iPrj) = rotTMP./720;
%       rotBgSubPS(:,:,iPrj) = rotBgSubPS(:,:,iPrj) ./ ...
%                                      BH_movingRMS(rotBgSubPS(:,:,iPrj), ...
%                                      [nBgPix,nBgPix]);
      clear rotTMP
   end

    
                        
  end
  


  cccStorage = cell(3,1);
  for iRefine = 1:3   
    
    iRefine
    if iRefine == 1
      % Initialize the best defocus from the global estimate for the first iter.
      maxDef = zeros(d3,1)+TLT(1,15);
      if (flgAstigmatism)
        maxAst = zeros(d3,1,'gpuArray')+TLT(1,12);
        maxAng = zeros(d3,1,'gpuArray')+TLT(1,13);
      else
        % Only using rotational average, so don't consider previosuly estimated
        % astigmatism.
        maxAst = zeros(d3,1,'gpuArray');
        maxAng = zeros(d3,1,'gpuArray');
      end
    end

    for iPrj = 1:d3
      defSearch{iRefine,iPrj} = maxDef(iPrj)-defRange{iRefine}:defInc{iRefine}:maxDef(iPrj)+defRange{iRefine};
    end
    
    cccStorage{iRefine} = zeros(length(defSearch{iRefine,1}),d3);
  

    for iDF = 1:length(defSearch{iRefine,1})
      if iRefine == 1
        % On first pass, search range is the same for all projections, so limit
        % calcs.
        df1 =  defSearch{iRefine,iPrj}(iDF) - 0;%maxAst(iPrj);
        df2 =  defSearch{iRefine,iPrj}(iDF) + 0;%maxAst(iPrj);
        [ Hqz ] = BH_ctfCalc(radialForCTF,Cs,WAVELENGTH, ...
                             [df1,df2,0],size(radialForCTF{1}), ...
                              AMPCONT,-1.0); 
      end

      for iPrj = 1:d3
        
        if iRefine > 1
    
        df1 =  defSearch{iRefine,iPrj}(iDF) - maxAst(iPrj);
        df2 =  defSearch{iRefine,iPrj}(iDF) + maxAst(iPrj);
        [ Hqz ] = BH_ctfCalc(radialForCTF,Cs,WAVELENGTH, ...
                             [df1,df2,maxAng(iPrj)],size(radialForCTF{1}), ...
                              AMPCONT,-1.0); 
        end
        
        % If not calc astigmatism use rotationally averaged always.
        if (iRefine == 1 || ~(flgAstigmatism)) && rotationalAvg
          [ iCCC ] = calc_CCC(radialForCTF{1},highCutoff, rotBgSubPS(:,:,iPrj), Hqz,0);
        else
          [ iCCC ] = calc_CCC(radialForCTF{1},highCutoff, bgSubPS(:,:,iPrj), Hqz,0);
        end
        cccStorage{iRefine}(iDF,iPrj) = gather(iCCC);
      end

      
    end
    % get the max scores for this iteration per projection
    maxDef = zeros(d3,1);
    for iPrj = 1:d3
      [~,c]=max(cccStorage{iRefine}(:,iPrj));
      maxDef(iPrj) = defSearch{iRefine,iPrj}(c);
    end
      

    if (flgAstigmatism)
      % If not leave maxAstig and maxAngle set to their initial values from the
      % TLT geometry.
      for iPrj = 1:d3
        n=1;
        iPrj
        for iAng = astigAngSearch{iRefine}
          for iDelDF = astigDefSearch{iRefine}
            iDelDfFull = iDelDF + maxAst(iPrj);
            iAngFull = iAng + maxAng(iPrj);

            df1 =  maxDef(iPrj) - iDelDfFull;
            df2 =  maxDef(iPrj) + iDelDfFull;

            % For values very close to zero, the search range may include
            % values |df1| < |df2| which is against convention.
            if abs(df1) >= abs(df2)           

              [ Hqz ] = BH_ctfCalc(radialForCTF,Cs,WAVELENGTH, ...
                                  [df1,df2,iAngFull],size(radialForCTF{1}), ...
                                  AMPCONT,-1.0);            

              [ iCCC ] = calc_CCC(radialForCTF{1},highCutoff, bgSubPS(:,:,iPrj), Hqz,1);                                           

              %fprintf('%d / %d coarse astigmatism search\n',n,size(initAstigCCC,1));
            else

              iCCC = -9999

            end

            astigCCC{iRefine}{iPrj}(n,:) = [iAngFull,iDelDfFull,iCCC]; 
            n = n + 1;      
          end
        end
      end

      for iPrj = 1:d3
        [~,c]=max(astigCCC{iRefine}{iPrj}(:,3));
        maxAst(iPrj) = astigCCC{iRefine}{iPrj}(c,2);
        maxAng(iPrj) = astigCCC{iRefine}{iPrj}(c,1);
      end

       maxAst
       maxAng
    end
  end
    
    

    
  


avgCCC=0;
for iPrj = 1:d3
  avgCCC = avgCCC+maxDef(iPrj);
end
avgCCC = avgCCC ./ d3;


 
  for iPrj = 1:d3

    defAstig = [maxDef(iPrj) - maxAst(iPrj), maxDef(iPrj) + maxAst(iPrj),maxAng(iPrj) ];
    [ Hqz ] = BH_ctfCalc(radialForCTF,Cs,WAVELENGTH,defAstig,paddedSize,AMPCONT,-1.0);
 
                                         
    [ ~, bandpass ] = calc_CCC(radialForCTF{1},highCutoff, bgSubPS(:,:,iPrj), Hqz, 1);
    
    bgSubPS(:,:,iPrj) = bgSubPS(:,:,iPrj) .* bandpass;
    if (rotationalAvg)
     rotBgSubPS(:,:,iPrj) = rotBgSubPS(:,:,iPrj) .* bandpass;
    end
  end
  
  

SAVE_IMG(MRCImage(gather(bgSubPS)),sprintf('fixedStacks/ctf/%s_bgOUT.mrc',fileName));
clear bgSubPS
if (rotationalAvg)
  SAVE_IMG(MRCImage(gather(rotBgSubPS)),sprintf('fixedStacks/ctf/%s_bgOUT_rotAvg.mrc',fileName));
end
clear rotBgSubPS

ang = TLT(:,4);
r1 = maxDef(TLT(:,1));
r2 = fit(ang, r1, 'smoothingSpline');
r3 = fit(ang, r2(ang), 'smoothingSpline');

ast1 = fit(ang,gather(maxAst(TLT(:,1))), 'smoothingSpline');
ang1 = fit(ang,gather(maxAng(TLT(:,1))), 'smoothingSpline');

      
figure('Visible','off'), ...
plot(ang,r3(ang),'bo',ang,zeros(1,d3)+TLT(1,15),'b--', ...
     ang,zeros(1,d3)+TLT(1,15)+defRange{1},'k--',...
     ang,zeros(1,d3)+TLT(1,15)-defRange{1},'k--');
  title(sprintf('CTF refine\nmeanDef %03.3f μm\nmeanDefPerTilt %03.3f μm ', TLT(1,15)*10^6,avgCCC*10^6));
  xlabel('Projection'); ylabel('defocus');
  
  saveas(gcf,sprintf('fixedStacks/ctf/%s_refine.pdf',fileName), 'pdf');

  if (flgAstigmatism)
    figure('Visible','off'), ...
    plot(ang,ast1(TLT(:,4)).*10^9,'bo',ang,(180/pi).*ang1(TLT(:,4)),'go');
      title('CTF refine astigmatism');
      xlabel('Projection'); ylabel('astig(nm) angle(deg)');

      saveas(gcf,sprintf('fixedStacks/ctf/%s_astig.pdf',fileName), 'pdf');
  end

    TLT(:,15) = r3(TLT(:,4));
    if (flgAstigmatism)
      TLT(:,12) = ast1(TLT(:,4));
      TLT(:,13) = ang1(TLT(:,4));
    end
    fileID = fopen(sprintf('fixedStacks/ctf/%s_ali%d_ctf_refine.tlt',STACK_PRFX,mapBackIter+1), 'w');
    fprintf(fileID,['%d\t%08.2f\t%08.2f\t%07.3f\t%07.3f\t%07.3f\t%07.7f\t%07.7f\t',...
             '%07.7f\t%07.7f\t%5e\t%5e\t%5e\t%7e\t%5e\t%5e\t%5e\t%5e\t%5e\t',...
             '%d\t%d\t%d\n'], TLT');
    fclose(fileID);   
 
end
clear

for i = 1:gpuDeviceCount
  gpuDevice(i);
end
end





function [ iCCC, bandpass ] = calc_CCC(radialForCTF,highCutoff, bgSubPS, Hqz, flgAst)


    ctfSQ = abs(Hqz);
    
    rV =  Hqz(1+size(Hqz, 1)/2,1+size(Hqz, 1)/2:end);
    freqVector  = radialForCTF(1+size(Hqz, 1)/2,1+size(Hqz, 1)/2:end);
    firstZero = find(rV > 0, 1,'first');
    secondZero = find(rV(firstZero+1:end) < 0, 1,'first') + firstZero;
    
    [~,firstMax]=min(abs(rV(firstZero:secondZero-1)- ...
                         rV(firstZero+1:secondZero))) ;
    firstMax = firstMax + firstZero;
    
%     if (flgAst)
%       bandpass = ( radialForCTF > freqVector(firstZero) & ...
%                    radialForCTF < highCutoff & ctfSQ < 3*rms(ctfSQ(:)));
%     else
% % %       bandpass = ( radialForCTF > freqVector(firstMax) & ...
% % %                    radialForCTF < highCutoff );
%     end
       
       lowRes =  radialForCTF > freqVector(firstZero) & radialForCTF < freqVector(firstMax);
       useRes =  radialForCTF > freqVector(firstZero) & radialForCTF < highCutoff;
       bandpass = single( useRes );
       bandpass(lowRes) = ctfSQ(lowRes).^4;
       
      
       ctfSQ = ctfSQ .* bandpass;
       
       bgSubPS = bgSubPS .* bandpass;
     

    
    iCCC = sum(sum((bgSubPS(useRes) .* ...
                    ctfSQ(useRes)))) ./ ...
                   ( numel(ctfSQ(useRes)).*...
                                        std2(ctfSQ(useRes)).*...
                                        std2(bgSubPS(useRes)) );

 
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
       doSplineInterp=1;

          mag =  (1+iDeltaZ./DFo).^0.5;
          
          estSize = 2048;
          ctf1 = BH_ctfCalc(pixelSize,TLT(iPrj,17),TLT(iPrj,18),DFo,estSize,TLT(iPrj,19),-1,1);
          ctf2 = BH_ctfCalc(pixelSize,TLT(iPrj,17),TLT(iPrj,18),iDeltaZ+DFo,estSize,TLT(iPrj,19),-1,1);
          ctf1 = ctf1(1:estSize/2);
          ctf2 = ctf2(1:estSize/2);
          firstZero = find(ctf1 > 0, 1, 'first');
          secondZero= find(ctf1(firstZero:end) < 0 , 1, 'first') + firstZero - 1;

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
         
         
%          testNoRefine = 0;
         if (testNoRefine)
           mag = 1
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
                                         % Using singleTaper here produces
                                         % a grid like artifact. Test
                                         % switch for EMC functions
          

   
          
    
        end
      end
    end
  
    
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
