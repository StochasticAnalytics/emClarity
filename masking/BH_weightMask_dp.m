function [ ctfWeights ] = BH_weightMask_dp(subTomoMeta, SIZE, SAMPLING,...
                                           GEOMETRY, PRECISION, METHOD,...
                                           varargin)
% Calculate the sampling of a number of subtomograms.
%
% If tiltGeometry, reconGeometry, and sTgeometry are structs, then calculate the
% weights for the full data set. Otherwise calculate a subset of weights
ctfScaleFactor = 1;
if nargin > 6
  ctfScaleFactor = varargin{1};
end
OverSIZE = 512;%512
if all(SIZE > 0)
  outputScaling = SIZE(1)/(512);
  SIZE = OverSIZE.*[1,1,1];%[512,512,512] ;
else
  error('Use regular weightmask for template matching for now');
  %optional override for template matching which may not always be cubic, wich
  %trades a little accuracy in the ctf mask for speed.
%   SIZE = abs(SIZE);
%   outputScaling = 1;
end

flgFirstPass = 1;

% Use this to control whether a group of ctfs for a single tomogram are made, or
% alternatively weights for the whole data set are calculated.
tomoList = GEOMETRY{1};
geometryFull = GEOMETRY{2};
nTomos = length(tomoList);

if (nTomos > 1)
  calcAllWeights = 1
  nCtfGroups = 1;
  ctfWeights = cell(2,1);
else
  calcAllWeights = 0;
  nCtfGroups = subTomoMeta.('ctfGroupSize').(tomoList{1})(1);
  ctfWeights = cell(1,nCtfGroups);
  
end

% Calculate at 512^3 and then reduce, so loop over nCtfGroups, pulling to main
% memory as you go.
[ rWeight ] = calc_rWeight( SIZE, PRECISION, METHOD);
TLT = subTomoMeta.('tiltGeometry').(tomoList{1});
pixelSize = TLT(1,16).*SAMPLING;
% Calc on the GPU if asked but pull to main memory so larger size wgts can
% be handled.
%[ exposureFilter ] = BH_exposureFilter( [SIZE(1:2)], TLT, METHOD, SAMPLING, 1  );   
%exposureFilter = gather(exposureFilter);

% % % if strcmp(PRECISION,'double')
% % %   exposureFilter = double(exposureFilter);
% % % end

[radialCTFCalc, phiCTFCalc, ~, ~, ~, ~ ] = ...
                  BH_multi_gridCoordinates(SIZE(1:2),'Cylindrical', ...
                                           METHOD, {'none'},1,1,0); 

% % % if strcmp(PRECISION, 'single')
% % %  radialCTFCalc = {single(radialCTFCalc ./ pixelSize),1,single(phiCTFCalc)};
% % % else
  % This should be the default.
  radialCTFCalc = {double(radialCTFCalc ./ pixelSize),1,double(phiCTFCalc)};
% % % end
% % % clear phiCTFCalc

for iCtfGroup = 1:nCtfGroups
  rec = cell(1+calcAllWeights,1);
  if strcmp(METHOD, 'GPU')
    rec{1} = zeros(SIZE,PRECISION,'gpuArray');
    if (calcAllWeights)
      rec{2} = zeros(SIZE,PRECISION,'gpuArray');
    end
    useGPU = 1;
    SIZE = gpuArray(SIZE);
  elseif strcmp(METHOD, 'cpu')
    rec{1} = zeros(SIZE,PRECISION);
    if (calcAllWeights)
      rec{2} = zeros(SIZE,PRECISION);
    end
    useGPU = 0;
  else
    error('METHOD must be  GPU or %s\n', 'cpu');
  end
  for iTomo = 1:1+(nTomos-1)*(calcAllWeights)


    reconGeometry = subTomoMeta.('reconGeometry').(tomoList{iTomo})./SAMPLING;
    reconShift = reconGeometry(2,:);
    
    % Also there are no offsets here, but this should be considered as in
    % syntheticMapback - also need to update things to save the tilt header.
    % Possibly just calculate the respective origins and lowerLeft vol as part of
    % the initialization.

    originPrj = subTomoMeta.('tiltGeometry').(tomoList{iTomo})(1,20:22)./SAMPLING;
    originPrj = ceil((originPrj+1)./2);
    originPrj(3) = 1;
    
    originVol = ceil((reconGeometry(1,1:3)+1)./2);
    lowerLeftVol = originPrj+reconShift-originVol;

    
    positionList = geometryFull.(tomoList{iTomo});
    positionList = positionList(positionList(:,26)~=-9999,:);
    TLT = subTomoMeta.('tiltGeometry').(tomoList{iTomo});
    nPrjs = size(TLT,1);
    nSubTomos = size(positionList,1) ;   

    


    % As of 20170109 CTF params are the same across a tilt series in the metadata.
    %defocus = TLT(1,15); defocus is one that could vary though.
    Cs = TLT(1,17);
    wavelength = TLT(1,18);
    ampContrast = TLT(1,19); % If experimenting with tilt dependend ampContrast this would not be correct
    
    [ exposureFilter ] = BH_exposureFilter( [SIZE(1:2)], TLT, METHOD, SAMPLING, 1  );
    if strcmp(PRECISION,'double')
      exposureFilter = double(exposureFilter);
    end
    for iSubTomo = 1:1+(nSubTomos-1)*(calcAllWeights)
      
      fprintf('%d/%d tomo %d/%d subtomo\n',iTomo,nTomos,iSubTomo,nSubTomos);
    %     xyzSubTomo = (positionList(iSubTomo,11:13)./1 + lowerLeftVol);
    %     prjVector = (xyzSubTomo- originPrj);
      
      if (calcAllWeights)
        xyzSubTomo = positionList(iSubTomo,11:13)./SAMPLING;       
        prjVector = xyzSubTomo - originVol + reconShift;
      else
        ctfGroupSize = subTomoMeta.('ctfGroupSize').(tomoList{iTomo})(2);
        
        xyzSubTomo = [(ctfGroupSize .* iCtfGroup) - ctfGroupSize/2, 0,0]./SAMPLING;
        prjVector = xyzSubTomo - originVol + reconShift;
        prjVector(2:3) = 0;
      end
     
      if ( calcAllWeights ) 
        iGold = positionList(iSubTomo,7);
      else
        iGold = 1;
      end
      for iPrj = 1:nPrjs


        data = rWeight;

        % The central section is located at -1.*tiltAngle in the 3dFT since we
        % rotate the specimen back by this amount in reconstruction.
        rTilt = BH_defineMatrix([90,TLT(iPrj,4),-90],'Bah','invVector');
        if (calcAllWeights)
          rSubTomo = reshape(positionList(iSubTomo,17:25),3,3);
        else
          rSubTomo = eye(3);
        end
        % To calculate the defocus, we need to rotate to where the projection is
        % in real space TLT(iPrj,4)

        rProjection = BH_defineMatrix([90,TLT(iPrj,4),-90],'Bah','forwardVector');
        prjCoords = rProjection*prjVector';
        
        iDf = prjCoords(3).*pixelSize +TLT(iPrj,15);
        iDefocus = [iDf - TLT(iPrj,12), iDf + TLT(iPrj,12), TLT(iPrj,13)];
        % Note the transpose (=inverse since the rotation matrix is orthogonal) is
        % taken because the stored matrix is for interpolation
        r = rSubTomo'*rTilt;

        % By default these are single PRECISION, but we truncate them anyhow, so
        % leave as single.
        [X,Y,Z,~,~,~] = BH_multi_gridCoordinates([SIZE(1:2),1],'Cartesian',METHOD,...
                                  {'single',r,[0,0,0]','invVector',1,1},0,1,0);


         % assuming ampContrast = 0.1, using -0.15 results in a weight with 
        % (0.1^0.15)^2~ 0.5 at zero freqency. Allows some recovery of low freq without
        % creating too severe a blur
        [Hqz, HqzUnMod] = BH_ctfCalc(radialCTFCalc,Cs,wavelength,iDefocus,SIZE(1:2),ampContrast,-0.5,-1); 

         %Default is double
         if strcmpi(PRECISION,'single')
           Hqz = single(abs(Hqz.*HqzUnMod)).^ctfScaleFactor;
         elseif strcmpi(PRECISION, 'double')
           Hqz = double(abs(Hqz.*HqzUnMod)).^ctfScaleFactor;
         end
       
         if SAMPLING > 1
           Hqz = Hqz .^ (1*(SAMPLING - 1)^-3);           
         end

        data = data.*Hqz.*exposureFilter(:,:,TLT(iPrj,1));

        data = data(:);

        % Shift from image to array coordinates and set any out of bounds values to
        % the origin.
        originXYZ = ceil(((SIZE + 1)./2));
        X = X(:)+originXYZ(1); 
        Y = Y(:)+originXYZ(2); 
        Z = Z(:)+originXYZ(3); 


        outOfBounds = logical(( X < 1 | X > SIZE(1) ) + (Y < 1 | Y > SIZE(2)) + (Z < 1 | Z > SIZE(3)));
        X(outOfBounds) = 1;%originXYZ(1);
        Y(outOfBounds) = 1;%originXYZ(2);
        Z(outOfBounds) = 1;%originXYZ(3);



        l = sub2ind(SIZE,round(X),round(Y),round(Z));
        clear X Y Z
        % much faster and effecitive for this simple sort.
        [B,I] = sort(l);
        idu = I(logical(B(1:end-1) - B(2:end)));

        tiltScale = 1- ( abs(sind(TLT(iPrj,4))).*0.2 );
        rec{iGold}(l(idu)) = rec{iGold}(l(idu)) + tiltScale.*data(idu);
        clear l idu
  
      end % end loop over projections

      for iGold = 1:1 + calcAllWeights
        % Zero out the large value from out of bounds conditions
        rec{iGold}(1) = 0;
      end

    end % end loop over subTomos
  end %end loop over Tomos



  g = BH_multi_gaussian3d(16.*[1,1,1],1.5);
  gf=fftshift(fftn(ifftshift(g)));
  g=real(fftshift(ifftn(ifftshift(gf.^3)))); clear gf
  g = g ./ sum(g(:));
  
  if (useGPU)
    g = gpuArray(g);
  end
 %SAVE_IMG(MRCImage(gather(rec{1})),'tmpR1.mrc');
 %SAVE_IMG(MRCImage(gather(rec{2})),'tmpR2.mrc');
  for iGold = 1:1+calcAllWeights
      
    ctfWeights{iGold,iCtfGroup} = convn(single(rec{iGold}),g,'same'); rec{iGold} = [];
%     ctfWeights{iGold,iCtfGroup} = convn(ctfWeights{iGold,iCtfGroup},g,'same'); 
%     ctfWeights{iGold,iCtfGroup} = convn(ctfWeights{iGold,iCtfGroup},g,'same');
    
  
    if (outputScaling ~= 1)
      ctfWeights{iGold,iCtfGroup} = BH_reScale3d( ctfWeights{iGold,iCtfGroup}...
                                      ,'', sprintf('%f',outputScaling), METHOD);
    end
    
    if (flgFirstPass)
      % Only calc this once.
      rad = fftshift(BH_bandpass3d(size(ctfWeights{iGold,iCtfGroup}),...
                                                       0,0,0,'cpu','nyquist'));
      flgFirstPass = 0;                                               
    end
      
    ctfWeights{iGold,iCtfGroup} = gather(ctfWeights{iGold,iCtfGroup}) .* rad;
    ctfWeights{iGold,iCtfGroup} = (ctfWeights{iGold,iCtfGroup} - ...
                          min(ctfWeights{iGold,iCtfGroup}(rad > 0.0)))./((nPrjs)./2-0.5);
  

  end


end % end loop over ctfGroups
clear rad g rec Hqz HqzUnMod B I data outOfBounds
end

function  [ rWeight ] = calc_rWeight( SIZE, PRECISION, METHOD)

    rWeight = ((abs([-1*floor((SIZE(1))/2):0,1:floor((SIZE(1)-1)/2)])'));
    if strcmp(METHOD,'GPU')
        rWeight = gpuArray(rWeight);
    end
    rOrig = ceil((SIZE(1)+1)./2);
    % % % % imod tilt zero freq = 0.2 * first non zero component
    rWeight(rOrig ) = 0.2;
    [rCut] = find(rWeight == floor(0.45*SIZE(1)));
    pixelFallOff = rCut(1) ;
    taperLow = 0.5+0.5.*cos((((1:pixelFallOff)).*pi)./(length((1:pixelFallOff+1))));

    pixelFallOff = SIZE(1)-rCut(2)+1 ;
    taperTop = 0.5+0.5.*cos((((1:pixelFallOff)).*pi)./(length((1:pixelFallOff+1))));
    rWeight(1:rCut(1)) = rWeight(1:rCut(1)).*flip(taperLow)';
    rWeight(rCut(2):end) = rWeight(rCut(2):end).*taperTop';


    %rWeight = rWeight + 1./rWeight.^2;
    % resample2d only handles scaling right now, so pad to z=3
    rWeight = repmat((rWeight), 1, SIZE(2),1);
    if strcmpi(PRECISION,'single')
      rWeight = single(rWeight);
    else
      % This should be the default.
      rWeight = double(rWeight);
    end
end

