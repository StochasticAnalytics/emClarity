function [ IMG_OUT, iPixelHeader, iOriginHeader, imgExt ] = BH_multi_loadOrBin( input_tilt_series_filename, samplingRate, DIMENSION, LOAD_VOL )
%Check to see if a cached binned image exists, either load or bin and load.
%   Switched to using imod's newstack and binvol to create binning and
%   removed inline binning from my workflow.

iPixelHeader = '';
iOriginHeader = '';
imgExt = '';
flgLoad = LOAD_VOL;
IMG_OUT = '';


try
  [imgPath, imgName, imgExt] = fileparts(input_tilt_series_filename);
catch
  input_tilt_series_filename
  error('Trouble getting fileparts for this tilt-series');
end

if isempty(imgPath)
  imgPath = '.';
end

rng('shuffle');
randIDX = randi([1,10^10],1);

if samplingRate > 1
  nameOUT = sprintf('cache/%s_bin%d%s', imgName, samplingRate, imgExt);
  doCalc = 0;
  if exist(nameOUT,'file')
    fprintf('Using cached file %s_bin%d%s\n', imgName, samplingRate,imgExt);
    [checkHeader,~] = system(sprintf('header %s > /dev/null',nameOUT));
    if (checkHeader)
      fprintf('File exists but appears to be corrupt %s_bin%d%s\n', imgName, samplingRate, imgExt);
      doCalc = 1;
    end
  else
    doCalc = 1;
  end
  
  if (doCalc)

    !mkdir -p cache
    switch DIMENSION
      case 3
        system(sprintf('binvol -BinningFactor %d -antialias 6 %s cache/%s_bin%d%s >  /dev/null', ...
                      samplingRate, input_tilt_series_filename, imgName, samplingRate, imgExt));
      case 2
        try
          tiltObj = MRCImage(input_tilt_series_filename,0);
        catch
          fprintf('If you have restarted somewhere it is possible the aligned stack is not found?\n');
          fprintf('Perhaps the value for CurrentTomoCpr is not correct?\n');
          error('Could not init an MRCImage for %s in loadOrBin',input_tilt_series_filename);
        end
        iHeader = getHeader(tiltObj);
        outputName = (sprintf('cache/%s_bin%d%s', imgName, samplingRate, imgExt));

        iPixelHeader = [iHeader.cellDimensionX/iHeader.nX .* samplingRate , ...
                        iHeader.cellDimensionY/iHeader.nY .* samplingRate, ...
                        iHeader.cellDimensionZ/iHeader.nZ .* samplingRate];
        
        iOriginHeader= [iHeader.xOrigin ./ samplingRate, ...
                        iHeader.yOrigin ./ samplingRate, ...
                        iHeader.zOrigin ./ samplingRate];
        
        pixelSize = iHeader.cellDimensionX/iHeader.nX; % Assuming X/Y the same and Z might be incorrect.
        
        % FIXME: forcing odd bin size so that transformations are the same as IMOD where the origin is between pixels for
        % even sized images.
        force_odd_dimension = true;
        [binSize, binShift] = BH_multi_calcBinShift([iHeader.nX, iHeader.nY], samplingRate, force_odd_dimension);

        % Gridding correction for the interpolation in the binning. Not
        % sure this is quite right, but it looks much better. TODO FIXME
        [ R ] = BH_multi_gridCoordinates([iHeader.nX,iHeader.nY],'Cartesian','GPU', {'none'},1,1,1);
        R = sinc(R).^-2;
        
        binSize = [binSize,iHeader.nZ];
        newStack = zeros(binSize,'single');
        for iPrj = 1:binSize(3)
          iProjection = gpuArray(OPEN_IMG('single',tiltObj,[],[],iPrj,'keep'));
          
          if (iPrj == 1)
            bhF = fourierTransformer(iProjection,'OddSizeOversampled');
          end
          
          iProjection = bhF.invFFT(bhF.fwdFFT(R.*iProjection,0,0,[1e-6,600,samplingRate*pixelSize,pixelSize]),2);
          
          iProjection = BH_resample2d(iProjection,[0,0,0],binShift,'Bah','GPU','forward',1/samplingRate,binSize(1:2),bhF);
          
          newStack(:,:,iPrj) = gather(iProjection);
        end
        
        SAVE_IMG(newStack,{outputName,'half'},iPixelHeader,iOriginHeader);
        clear newStack bpFilt iProjection
        %        system(sprintf('newstack -shrink %d -antialias 6 %s cache/%s_bin%d%s > /dev/null', ...
        %                                     samplingRate,input_tilt_series_filename, imgName, samplingRate,imgExt));
      otherwise
        error('DIMENSION should be 2 or 3\n')
    end
    
    
    
    
  end
  
  
  if (flgLoad)
    failedLoads = 0;
    while failedLoads < 3
      try
        fprintf('pwd is %s\n', pwd);
        fprintf(...
          'attempting to load cache/%s_bin%d%s\n', imgName, samplingRate,imgExt);
        m = MRCImage(sprintf(...
          'cache/%s_bin%d%s', imgName, samplingRate,imgExt));
        fprintf('Loaded the MRCImage\n');
        IMG_OUT = OPEN_IMG('single', m);
        fprintf('Loaded the volume\n');
        IMG_OUT = single(IMG_OUT);
        fprintf('Volume --> single\n');
        failedLoads = 3;
      catch
        failedLoads = failedLoads + 1
        pause(failedLoads.^3); % 1,8,27 second pauses
      end
    end
  end
else
  % So you could pass a -1 as sampling which would indicate to not load but that
  % syntax is only intended when resampling is required, so ignore the flag here
  % but throw a warning.
  fprintf('\n\nYou requested a sampling of -1 Nonsense!! loading anyway.\n\n');
  
  IMG_OUT = OPEN_IMG('single', input_tilt_series_filename);
end



end

