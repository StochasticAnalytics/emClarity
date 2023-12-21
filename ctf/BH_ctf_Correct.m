function [ ] = BH_ctf_Correct( PARAMETER_FILE, STACK_PRFX )
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
  mapBackIter = subTomoMeta.currentTomoCPR; 
catch
  mapBackIter = 0;
end

if isnan(str2double(STACK_PRFX))
  % It is a name, run here.
  nGPUs = 1;
  flgParallel = 0;
  STACK_LIST = {STACK_PRFX};
  ITER_LIST = {STACK_LIST};
else
  flgParallel = 1;
  nGPUs = pBH.('nGPUs');
  
  STACK_LIST_tmp = fieldnames(subTomoMeta.mapBackGeometry);
  STACK_LIST_tmp = STACK_LIST_tmp(~ismember(STACK_LIST_tmp,'tomoName'));
  ITER_LIST = cell(nGPUs,1);
  nST = 1; STACK_LIST = {};
  for iStack = 1:length(STACK_LIST_tmp)
    if subTomoMeta.mapBackGeometry.(STACK_LIST_tmp{iStack}).nTomos
      STACK_LIST{nST} = STACK_LIST_tmp{iStack};
      nST = nST + 1;
    end
  end
  clear STACK_LIST_tmp
  for iGPU = 1:nGPUs
    ITER_LIST{iGPU} = STACK_LIST(iGPU:nGPUs:length(STACK_LIST));
  end
end

pixelSize = pBH.('PIXEL_SIZE');
!mkdir -p ctfStacks

try 
  EMC_parpool(nGPUs);
catch 
  delete(gcp('nocreate'));
  EMC_parpool(nGPUs);
end


parfor iGPU = 1:nGPUs
  
  if ( flgParallel )
    useGPU = iGPU;
    gpuDevice(useGPU);
  else
    useGPU = BH_multi_checkGPU(-1);
    gpuDevice(useGPU);
  end
  
  for iTilt = 1:length(ITER_LIST{iGPU})
    

    STACK_PRFX = ITER_LIST{iGPU}{iTilt};

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
 
    inputStack = sprintf('aliStacks/%s_ali%d.fixed',STACK_PRFX,mapBackIter+1);
    outputStack = sprintf('ctfStacks/%s_ali%d_ctf.fixed',STACK_PRFX,mapBackIter+1);
    iMrcObj = MRCImage(inputStack,0);



    iHeader = getHeader(iMrcObj);
    iPixelHeader = [iHeader.cellDimensionX/iHeader.nX, ...
                    iHeader.cellDimensionY/iHeader.nY, ...
                    iHeader.cellDimensionZ/iHeader.nZ];

    d1 = iHeader.nX;
    d2 = iHeader.nY;
    nPrjs = iHeader.nZ;
    
    correctedStack = zeros(d1,d2,nPrjs,'single');
    
    for iPrj = 1:nPrjs
      
      CS = TLT(iPrj,17);
      WL = TLT(iPrj,18); 
      AMPCONT = TLT(iPrj,19);      
      ddF = TLT(iPrj,12);
      dPhi = TLT(iPrj,13);
      D0 = TLT(iPrj,15);


      fastFTSize = BH_multi_iterator([d1,d2],'fourier2d');
      padVal = BH_multi_padVal([d1,d2],fastFTSize);
      trimVal = BH_multi_padVal(fastFTSize,[d1,d2]);


      initImg = randn(fastFTSize,'single','gpuArray');
      f   = FFT(initImg);

      ctf = CTF(fastFTSize,pixelSize*10^10,'GPU');
      maxZ = 500;
      maxEval = cosd(TLT(iPrj,4)).*(d1/2) + maxZ./2*abs(sind(TLT(iPrj,4)));
      oX = ceil((d1+1)./2);
      oY = ceil((d2+1)./2);
      iEvalMask = floor(oX-maxEval):ceil(oX+maxEval);
  
      STRIPWIDTH = 512;
      STRIPWIDTH = STRIPWIDTH + mod(STRIPWIDTH,2);
      % take at least 1200 Ang & include the taper if equal to STRIPWIDTH
      tileSize   = floor(max(600./pixelSize, STRIPWIDTH + 28));
      tileSize = tileSize + mod(tileSize,2);
      %fprintf('stripwidth tilesize %d %d\n',STRIPWIDTH,tileSize);
      incLow = ceil(tileSize./2);
      incTop = tileSize - incLow;
  

      iProjection = BH_padZeros3d(getVolume(iMrcObj,[-1],[-1],TLT(iPrj,23),'keep'), ...
                                            padVal(1,:),padVal(2,:),'GPU','singleTaper');

      correctedPrj = zeros([d1,d2],'single','gpuArray');
      iProjectionFT = f.fwdFFT(iProjection);
  

      stripDefocusOffset = floor(STRIPWIDTH/2);
      for i = 1: STRIPWIDTH : d1

        if (i+tileSize-1) < d1
            endIDX = (i+tileSize-1);
            endCUT = i + STRIPWIDTH - 1 + 7;  
            trimmedSIZE = STRIPWIDTH;
        elseif any(ismember(i:d1,iEvalMask))
            endIDX = d1;
            endCUT = d1;
            trimmedSIZE = endCUT-i+1 -7;
        end    

        % The eval mask condition can be replaced once the per tomo condition
        % is trusted.
        if any(ismember(i:endIDX,iEvalMask))


          DF = D0 +(i + stripDefocusOffset - oX)*pixelSize*-1.*tand(TLT(iPrj,4));


          if ~( isempty(DF) )

            iDefocus = [DF - ddF, DF + ddF, dPhi];

             if pixelSize < 2.0e-10
               % use double precision - this is not enabled, but needs to be -
               % requires changes to radial grid as well.           
               ctf.new_img(iDefocus,CS,WL,AMPCONT,-1,-1);
             else
               ctf.new_img(iDefocus,CS,WL,AMPCONT,-1);
             end


            tile  = ctf.multiply(iProjectionFT);

            tile = BH_padZeros3d(real(f.invFFT(tile,2)), ...
                                 trimVal(1,:),trimVal(2,:),'GPU','single');

            % trim prior to pulling off gpu to minimize xfer
          else

            % No particles in this strip, so just replace with simple inversion
            % to keep the global image statistics ~ correct.

            tile = -1.*BH_padZeros3d(iProjection, trimVal(1,:),trimVal(2,:),'GPU','single');

          end

            %correctedStack(i + 7 : endCUT,:,TLT(iPrj,1)) = ...
             %             gather(tile(8:trimmedSIZE+7,:));

          correctedPrj(i:endIDX,:) = tile(i:endIDX,:); 

        else
        %fprintf('ignoring strip centered on %d for prj %d',i,TLT(iPrj,1));
        end
      end % end loop over strips
      correctedStack(:,:,TLT(iPrj,1)) = gather(correctedPrj);

    end % end loop over prjs
    
    SAVE_IMG(MRCImage(correctedStack), outputStack,iPixelHeader);

  end % end loop over tilt-series
end % end parfor

delete(gcp('nocreate'));


end

