function [ OUTPUT ] = BH_multi_iterator( SIZES, OPERATION )
%Determine and add padding needed to handle large images iteratively
%   Detailed explanation goes here

% a rough balance between increasing size, and speed of fft, made this at home
% using cpu, update with GPU at work.
nB = [0,0,0];
% nextBest = [64,72,96,108,128,144,160,168,180,192,216,224,256,...
%             270,288,300,320,336,360,384,400,432,448,480,512];
nextBest = [64,72,80,84,90,96,108,112,120,126,128,144,160,162, ...
  168,180,192,216,224,240,256,270,288,320,324,336,...
  360,378,384,400,416,432,448,480,486,504,512];
% fall off between 1 and 0 in masking function
APODIZATION = 2.*6;

switch OPERATION
  
  case 'fourier'
    
    if any(SIZES(1,:) < 0)
      flgDescend = 1;
      SIZES(1,:) = abs(SIZES(1,:));
    else
      flgDescend = 0;
    end
    % the target size, usually set by gpu/cpu limits
    sizeTarget= SIZES(1,:);
    if ( flgDescend )
      for i = 1:3
        try
          nB(i) = nextBest(find(nextBest <= sizeTarget(i), 1, 'last'));
        catch
          nB(i) = nextBest(1);
        end
      end
    else
      for i = 1:3
        try
          nB(i) = nextBest(find(nextBest >= sizeTarget(i), 1, 'first'));
        catch
          nB(i) = nextBest(end);
        end
      end
    end
    
    
    if ~all(nB)
      error('next best not found in range 128-512 for [%d,%d,%d]', sizeTarget);
    end
    OUTPUT = nB;
    
  case 'fourier2d'
    % Found by
    % % %     for i = 64:2:3838*2
    % % %       if (sum(factor(i).*(factor(i) >= 5))< 10)
    % % %       f = [f,i];
    % % %       end
    % % %     end
    nextBest =  [64,72,80,84,90,96,108,112,120,126,128,144,160,162, ...
      168,180,192,216,224,240,252,256,270,288,320,324,336,...
      360,378,384,432,448,480,486,504,512,540,576,640,648,...
      720,756,768,810,864,896,960,972,1008,1024,1080,...
      1134,1152,1280,1296,1344,1440,1458,1512,1536,1620,1728,...
      1792,1920,1944,2016,2048,2160,2268,2304,2430,2560,...
      2592,2688,2880,2916,3024,3072,3240,3402,3456,3584,...
      3840,3888,4032,4096,4320,4374,4536,4608,4860,5120,...
      5184,5376,5760,5832,6048,6144,6480,6804,6912,7168,7290,...
      7488,7560,7840,7920,8192];
    sizeTarget= SIZES(1,:);
    nB = [0,0];
    for i = 1:2
      try
        nB(i) = nextBest(find(nextBest >= sizeTarget(i), 1, 'first'));
      catch
        nB(i) = 0;
      end
    end
    
    if ~all(nB)
      % For some reason "error" only takes scalar args for print formating.
      fprintf('next best not found in range 64-7290 for [%d,%d,%d]\n',target);
      error('asdf')
    end
    OUTPUT = nB;
    
  case 'convolution'
    % This could be optimized automatically to balance increased target size vs
    % number of iterations/post padding.
    
    % the target size, to get the most out of calcs, spend as much time on the
    % gpu as possible. With finer angular searches, the number of references
    % needs more memory, so a smaller size here means more transfers, but this
    % should be balanced by the finer angles (more comp)
    if all(SIZES(1,:) == 256)
      nextBest = [128,144,160,168,192,216,224,256];
    elseif all(SIZES(1,:) == 384)
      nextBest = [128,144,160,168,192,216,224,256,...
        288,300,320,336,360,384];
    elseif all(SIZES(1,:) == 432)
      nextBest = [128,144,160,168,192,216,224,256,...
        288,300,320,336,360,384,400,432];
    elseif all(SIZES(1,:) == 512)
      nextBest = [128,144,160,168,192,216,224,256,...
        288,300,320,336,360,384,400,432,480,512];
    elseif all(SIZES(1,:) > 512)
      nextBest = [128,144,160,168,192,216,224,256,...
        288,300,320,336,360,384,400,432,480,512,...
        540,576,640,648,720,756,768,810,864,896,960,972,1008,1024];
    end
    sizeImage      = SIZES(2,:);
    sizeTemplate   = SIZES(3,:); % the mask or kernel
    sizeParticle   = SIZES(4,:) ;% a subregion of sizeTemplate
    
    
    borderSizeCalc = floor((sizeTemplate + APODIZATION)./2);
    borderSizeKeep = borderSizeCalc + 2.*sizeParticle;
    
    abs(sum(sizeImage - sizeTemplate))
    sum(0.1.*sizeImage)
    if abs(sum(sizeImage - sizeTemplate)) < sum(0.1.*sizeImage)
      
      OUTPUT = [[0,0,0];[0,0,0] ;sizeImage; ...
        sizeImage; sizeImage ; [1,1,1]];
      return
    end
    score = zeros(length(nextBest),6);
    
    
    
    validCalc = repmat(nextBest',1,3) - 2.*repmat(borderSizeCalc,length(nextBest),1);
    validKeep = repmat(nextBest',1,3) - 2.*repmat(borderSizeKeep,length(nextBest),1);
    minIter= floor(repmat(sizeImage,length(nextBest),1)./validKeep);
    postPad= repmat(nextBest',1,3)- ...
      (repmat(sizeImage+borderSizeKeep,length(nextBest),1) -minIter.*(validKeep+1));
    
    % Added this so I can work with test cases where the volume to be
    % searched is the same size as the reference
    
    score(:,2:4) = repmat(nextBest',1,3) ./ postPad .* (minIter >= 0)
    
    [~, cX] = max(score(:,2)) ;
    [~, cY] = max(score(:,3)) ;
    [~, cZ] = max(score(:,4))  ;
    
    chunkSize = nextBest([cX,cY,cZ]);
    
    cY = cY + length(nextBest);
    cZ = cZ + 2.*length(nextBest);
    validAreaKeep = validKeep([cX,cY,cZ]);
    validAreaCalc = validCalc([cX,cY,cZ]);
    nIters = minIter([cX,cY,cZ])+1;
    postPAD= postPad([cX,cY,cZ]);
    
    
    OUTPUT = [borderSizeKeep;postPAD ;chunkSize; ...
      validAreaKeep; validAreaCalc ; nIters];
    
  case 'binning'
    
    binFactor = SIZES(1,:);
    sizeImage = SIZES(2,:);
    
    % diminishing return on increasing size beyond 512^3
    sizeTarget = min([512,512,512],sizeImage);
    
    % Choose the largest size that does fast fft and leaves no remainder at the
    % sampling.
    OUTPUT = [0,0,0];
    for i = 1:3
      
      rVal = 1;
      useTarget=true;
      n = 1;
      while rVal ~= 0
        if (useTarget)
          fftVal = find(nextBest <= sizeTarget(i), 1, 'last');
          rVal = rem(nextBest(fftVal), binFactor(i));
          useTarget=false;
        else
          try
            rVal = rem(nextBest(fftVal-n), binFactor(i));
            n = n + 1;
          catch
            error('Did not find factor of binvalue for bin %d targetsize %d', binFactor(i),sizeTarget(i))
          end
        end
      end
      
      OUTPUT(i) = nextBest(fftVal);
    end
    
  case 'extrapolate'
    
    targetSize = SIZES(1,:);
    kVal = SIZES(2,:);
    
    finalVal = [0,0,0];
    padVal = [0,0,0];
    intVal=true;
    while intVal
      intVal=false;
      nF = 2.^-kVal.*(targetSize-(1-2.^kVal));
      for i = 1:3
        if rem(ceil(nF(i)),nF(i))
          intVal = true;
          targetSize(i) = targetSize(i) + 1;
          padVal(i) = padVal(i) + 1;
        else
          finalVal(i) = nF(i);
        end
      end
    end
    
    OUTPUT = [finalVal; padVal];
  otherwise
    error('OPERATION is case-sensitive fourier, convolution or binning, not %s', OPERATION);
end


end


