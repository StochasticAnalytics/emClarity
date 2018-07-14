function [ ] = BH_ctf_Correct( PARAMETER_FILE, STACK_PRFX, PRECISION, THICKNESS, NWORKERS )
%CTF correction for tilt series using general geometry.
%   Correct for the CTF using a local approach, similar to strip based
%   periodogram, but with tiles that are smaller allowing for arbitrary
%   defocus gradients. 
%
%   The full stack is corrected, st if only a small region is to be used,
%   it would be faster to have trimmed the stack. This should be done
%   before ctf estimation though, st the correct origin is included in the
%   tilt information.
%

pBH = BH_parseParameterFile(PARAMETER_FILE);
try
  load(sprintf('%s.mat', pBH.('subTomoMeta')), 'subTomoMeta');
  mapBackIter = subTomoMeta.currentTomoCPR; clear subTomoMeta
catch
  mapBackIter = 0;
end

nWorkers = str2num(NWORKERS);
% As in other places a more definitive test should be run to balance any benefit 
% in oversampling (smoother rings from catesian --> polar) vs "dilution" of power
% from the smaller tile. Also should make tile size depend on the defocus explicitly
% to be sure to inlcude all high res information.
ctfSize = 2048;
tileSize = 420;
usableArea = str2num(THICKNESS)


try
  % make sure there isn't a refined version first.
  TLTNAME = sprintf('fixedStacks/ctf/%s_ali%d_ctf_refine.tlt',STACK_PRFX,mapBackIter+1);
  TLT = load(TLTNAME)
  fprintf('using refined TLT %s\n', TLTNAME);
 
catch
  TLTNAME = sprintf('fixedStacks/ctf/%s_ali%d_ctf.tlt',STACK_PRFX,mapBackIter+1);
  TLT = load(TLTNAME);
  fprintf('using TLT %s\n', TLTNAME);
end

try
  applyExposureFilter = pBH.('applyExposureFilter');
catch
  applyExposureFilter = 0;
end

% % % CTF_STACK = sprintf('%s/%s_ctf%s',pathName,fileName,extension);
!mkdir -p ctfStacks
inputStack = sprintf('aliStacks/%s_ali%d.fixed',STACK_PRFX,mapBackIter+1);
outputStack = sprintf('ctfStacks/%s_ali%d_ctf.fixed',STACK_PRFX,mapBackIter+1);
iHeader = MRCImage(inputStack,0);
STACK = single(getVolume(iHeader));

iHeader = getHeader(iHeader);
iPixelHeader = [iHeader.cellDimensionX/iHeader.nX, ...
                iHeader.cellDimensionY/iHeader.nY, ...
                iHeader.cellDimensionZ/iHeader.nZ];


[d1,d2,d3] = size(STACK)
nPrjs = d3;

[ evalMask, deltaZ ] = BH_multi_projectionMask( [d1,d2,d3;usableArea], TLT, 'GPU' );




% Local normalization doesn't address any large scale gradients in the
% images. Do a simple high pass over the lowest 7 frequencyBinns
bandNyquist = BH_bandpass3d([d1,d2,1],0,0,0,'GPU','nyquistHigh');

taperMask = fspecial('gaussian',[9,9],1.5);
for iPrj = 1:nPrjs
  
%   fprintf('running local normalization on %d/%d projections size %d pixels.\n',iPrj,nPrjs,tileSize);
  fprintf('running normalization on %d/%d projections.\n',iPrj,nPrjs);

  if strcmp(PRECISION, 'double')
    iProjection = double(gpuArray(STACK(:,:,TLT(iPrj,1))));
    iMask = double(convn(gpuArray(single(evalMask(:,:,TLT(iPrj,1)))),taperMask,'same'));
  else
    iProjection = gpuArray(STACK(:,:,TLT(iPrj,1)));
    iMask = convn(gpuArray(single(evalMask(:,:,TLT(iPrj,1)))),taperMask,'same');
  end
  
  iProjection = iProjection - mean(iProjection(evalMask(:,:,TLT(iPrj,1))));
  
  iProjection  = real(ifftn(fftn(iProjection.*iMask).*bandNyquist));
%   iProjection = iProjection - mean(iProjection(evalMask(:,:,TLT(iPrj,1))));
  % Look for any extreme outliers

  inFin = ~(isfinite(iProjection)); nInf = sum(inFin(:));
  if (nInf)
   fprintf('Removing %d (%2.4f) inf from prj %d\n',nInf,100*nInf/numel(iProjection),TLT(iPrj,1));
   iProjection(inFin) = 0;
  end
  
  iRms = rms(iProjection(evalMask(:,:,TLT(iPrj,1))));
  outliers = (iProjection > 6 * iRms); nOutliers = sum(outliers(:));
  
  if (nOutliers)
    fprintf('Truncating %d (%2.4f) outliers from prj %d\n',nOutliers,100*nOutliers/numel(iProjection),TLT(iPrj,1));
    % Set the outliers to a random value between 0.5 and 1 * 3*rms
    iProjection(outliers) = 6.*iRms.* (rand(size(iProjection(outliers)))-0.5);

%     iProjection(outliers) = sign(iProjection(outliers)).*3.*iRms.*((rand(size(iProjection(outliers)))./2)+0.5);
%     iProjection = iProjection - mean(iProjection(evalMask(:,:,TLT(iPrj,1))));
    iProjection = iProjection ./ ( rms(iProjection(evalMask(:,:,TLT(iPrj,1)))) .* cosd(TLT(iPrj,4)).^1.5);
  else
    iProjection = iProjection ./ ( iRms .* cosd(TLT(iPrj,4)).^1.5 );
  end
      
   
  
  STACK(:,:,TLT(iPrj,1)) = gather(single(iProjection.*iMask));
  clear iProjection iMask
end


clear bandNyquist

% Push to gpu when initializing workers
deltaZ = gather(deltaZ);
evalMask = gather(evalMask);
% SAVE_IMG(MRCImage(STACK),'outliers.mrc');
% SAVE_IMG(MRCImage(STACK), outputStack,iPixelHeader);
% SAVE_IMG(MRCImage(single(deltaZ)),'deltaZ.mrc');
% SAVE_IMG(MRCImage(single(evalMask)),'evalMask.mrc');
% TLT

try 
   ppool = parpool(nWorkers);
catch 
    delete(gcp);
    ppool = parpool(nWorkers);
end


[exposureFilter] = BH_exposureFilter(ctfSize, TLT,'cpu',1, 0);

for iPrj = 1:nPrjs
   % TLT(iPrj,1)
    % Assuming for the moment no duplicate tilts
    %%%dose = dosePerTilt(find(dosePerTilt(:,1) == TLT(iPrj,1)),3)
    fprintf('using dose %3.3f for projection at %2.2f\n',TLT(iPrj,11), TLT(iPrj,4));
    pFuture(iPrj) = parfeval(ppool,@runCTF,2, STACK(:,:,TLT(iPrj,1)), ...
                                              deltaZ(:,:,TLT(iPrj,1)), ...
                                              evalMask(:,:,TLT(iPrj,1)), ...
                                              TLT,iPrj,ctfSize,tileSize,...
                                              PRECISION,...
                                              exposureFilter(:,:,TLT(iPrj,1)));
end


for i = 1:nPrjs
  i
  [iPrj, ctfCorr, errorOut] = fetchNext(pFuture);

  STACK(:,:,TLT(iPrj,1)) = ctfCorr;
  
 
end

SAVE_IMG(MRCImage(STACK), outputStack,iPixelHeader);

end

function [ sliceOUT, eout ] = runCTF(sliceIN, dZ, evalMask, TLT ,iPrj, ...
                                     ctfSize,tileSize,precision,...
                                     expFilter)
%Apply the ctf correction.

evalMask = gpuArray(evalMask);
dZ = gpuArray(dZ);

if strcmp(precision, 'double')
  expFilter = double(gpuArray(expFilter));
else
  expFilter = gpuArray(expFilter);
end


eout = 'no error';
% Size to pad the tile to, in order to help deal with the polar nature of
% the CTF interpolating onto a rectangular grid.
fprintf('%f\n',TLT(iPrj,15))
CTFSIZE = ctfSize;
ddF = TLT(iPrj,12);
dPhi = TLT(iPrj,13);
D0 = TLT(iPrj,15);
PIXEL_SIZE = TLT(iPrj,16); 
Cs = TLT(iPrj,17);
WAVELENGTH = TLT(iPrj,18);
AMPCONT = TLT(iPrj,19);



[radialGrid,phi,~,~,~,~] = BH_multi_gridCoordinates(CTFSIZE.*[1,1],'Cylindrical','GPU',{'none'},1,0,0);
if strcmpi(precision, 'double')
  radialGrid = {double(radialGrid)./PIXEL_SIZE,0,double(phi)};
elseif strcmpi(precision,'single')
  radialGrid = {radialGrid./PIXEL_SIZE,0,phi};
end

% Tile to cut out each time. Balance of accuracy due to gradient (smaller)
% and due to power of signal (larger.)

% Tile extension that is tailed to zero on each dimension (pixels). Assumed
% to be much smaller than the area we trash after each tile.
apoSize = 6;

% Optimize fft alg to size specific.
if strcmpi(precision, 'double')
  test_fft = randn([CTFSIZE, CTFSIZE], 'gpuArray');
else
  test_fft = randn([CTFSIZE, CTFSIZE],'single', 'gpuArray');
end
fftw('planner', 'patient');
fft2(test_fft);
clear test_fft


% Must be smaller than tile, ideally < 1/2 tile size. Width must be reduced
% according to tilt angle to adjust for defocus gradient. This is not
% optimized, and may need to be adjusted in the future. 

STRIPWIDTH = 2*floor((64*abs(cosd(TLT(iPrj,4))).^1.5/2));
STRIPVERT  = floor(tileSize./2);


% Since we can do this on a tilted image, the model of a plane wave interacting
% piecewise at different heights would then suggest a narrowing tile width
% perpendicular to the tilt axis...maybe? condiser.

sliceIN = (gpuArray(sliceIN));
sliceOUT = zeros(size(sliceIN),'single','gpuArray');   

                       
incLow = ceil(tileSize./2);
incTop = tileSize - incLow;
border = ceil(incLow+apoSize)+1;

for i = 1+border: STRIPWIDTH/2 : size(sliceIN,1) - border
  for j = 1+border: STRIPVERT/2: size(sliceIN,2) - border
    if (evalMask(i,j))

     
      tile = (sliceIN(i - incLow: i+ incTop -1,...
                            j - incLow: j+ incTop -1));
      

      padVal = (CTFSIZE - tileSize)/2;
      if strcmp(precision,'double')
        tile = BH_padZeros3d(tile, [padVal,padVal], [padVal,padVal], ...
                                                          'GPU', 'doubleTaper');
      else
        tile = BH_padZeros3d(tile, [padVal,padVal], [padVal,padVal], ...
                                                          'GPU', 'singleTaper'); 
      end
       
      DF = D0 + (double(dZ(i,j)) .* PIXEL_SIZE);
      iDefocus = [DF - ddF, DF + ddF, dPhi];
      Hqz = BH_ctfCalc(radialGrid,Cs,WAVELENGTH,iDefocus,[CTFSIZE,CTFSIZE],AMPCONT,-0.5);
     
      ampFact = -1.0;
      % For now just let it be phase flipping, re-visit if needed.
      if (ampFact <= 0)
        % Original flipping phase with sharp transitions, trying multiplying by
        % the sqrt of the flipped function (past first zero) to create a
        % smoother filter that also downweights information near ctf zeros.
        
        % ctfTile = single(real(ifftn(fftn(tile).*(sign(Hqz)).*envelope)));
        ctfTile = single(real(ifftn(fftn(tile).*Hqz.*expFilter)));
      else
        Hqz = Hqz ./ (Hqz.^2 + ampFact^2);
        ctfTile = single(real(ifftn(fftn(tile).*Hqz.*envelope)));
      end


      sliceOUT(i - STRIPWIDTH/2+1:i+STRIPWIDTH/2, ...
         j - STRIPVERT /2+1:j+STRIPVERT /2) = ...
         ctfTile((CTFSIZE - STRIPWIDTH)/2+1:(CTFSIZE+STRIPWIDTH)/2  , ...
                 (CTFSIZE - STRIPVERT)/2+1 :(CTFSIZE+STRIPVERT)/2 );
    
    end
  end
end

                
sliceOUT =  gather(sliceOUT); 

clear sliceIN radialGrid bandpass Hqz ctfTile
       

end








