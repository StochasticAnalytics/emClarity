function [ ctfWeights ] = BH_weightMask_dp(subTomoMeta, SIZE, SAMPLING,...
                                           GEOMETRY, PRECISION, METHOD,...
                                           varargin)
% Calculate the sampling of a number of subtomograms.
%
% If tiltGeometry, reconGeometry, and sTgeometry are structs, then calculate the
% weights for the full data set. Otherwise calculate a subset of weights
ctfScaleFactor = 1
if nargin > 6
  ctfScaleFactor = varargin{1};
end
OverSIZE = 512;

%512
% % % if all(SIZE > 0)
% % %   outputScaling = SIZE(1)/(512);
% % %   SIZE = OverSIZE.*[1,1,1];%[512,512,512] ;
% % % else
% % %   error('Use regular weightmask for template matching for now');
% % %   %optional override for template matching which may not always be cubic, wich
% % %   %trades a little accuracy in the ctf mask for speed.
% % % %   SIZE = abs(SIZE);
% % % %   outputScaling = 1;
% % % end
outputScaling = false;
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
% % % [ rWeight ] = calc_rWeight( SIZE, PRECISION, METHOD);

TLT = subTomoMeta.('tiltGeometry').(tomoList{1});
pixelSize = TLT(1,16).*SAMPLING;


% % % [radialCTFCalc, phiCTFCalc, ~, ~, ~, ~ ] = ...
% % %                   BH_multi_gridCoordinates(SIZE(1:2),'Cylindrical', ...
% % %                                            METHOD, {'none'},1,1,0); 
% % % 
% % % 
% % % % This should be the default.
% % % radialCTFCalc = {radialCTFCalc ./ pixelSize,1,phiCTFCalc};


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


% % %     reconGeometry = subTomoMeta.('reconGeometry').(tomoList{iTomo})./SAMPLING;
    reconGeometry = subTomoMeta.('reconGeometry').(tomoList{iTomo});

    reconShift = reconGeometry(2,:);
    
    % Also there are no offsets here, but this should be considered as in
    % syntheticMapback - also need to update things to save the tilt header.
    % Possibly just calculate the respective origins and lowerLeft vol as part of
    % the initialization.

% %     originPrj = subTomoMeta.('tiltGeometry').(tomoList{iTomo})(1,20:22)./SAMPLING;
% % %     originPrj = ceil((originPrj+1)./2);
    originPrj = subTomoMeta.('tiltGeometry').(tomoList{iTomo})(1,20:22);

% % %     originPrj(3) = 1;
    
    originVol = ceil((reconGeometry(1,1:3)+1)./2);
    lowerLeftVol = originPrj+reconShift-originVol;

    
    positionList = geometryFull.(tomoList{iTomo});
    positionList = positionList(positionList(:,26)~=-9999,:);
    TLT = subTomoMeta.('tiltGeometry').(tomoList{iTomo});
    nPrjs = size(TLT,1);
    nSubTomos = size(positionList,1) ;   

    zCoords = sort(positionList(:,13));
    

    % The thickness in the tomogram is NOT the thickness in the projection
    % if the specimen is tilted at zero-tilt angle. Need to think of a
    % better way to calculate this.
% %     if (length(zCoords) < 5)
% %       iThickness = 100
% %     else
% %       iThickness =  iqr(zCoords).*pixelSize*10^9
% % %       meanMax = mean(zCoords(1:5)).*pixelSize*10^9;
% % %       meanMin = mean(zCoords(end-5:end)).*pixelSize*10^9;
% %     end
% %     
% %     if iThickness > 400
% %       fprintf('capping thickness to 400\n');
% %       iThickness = 400;
% %     end
    iThickness = 75;

    

% % %     [ exposureFilter ] = BH_exposureFilter( SIZE(1:2), TLT, METHOD, SAMPLING, 1 );

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
% % %       for iPrj = 1:nPrjs
% % % 
% % %          
% % %         iCs = TLT(iPrj,17);
% % %         iWavelength = TLT(iPrj,18);
% % %         iPhaseShift = TLT(iPrj,19);   
% % %         iDefocus = TLT(iPrj,15);   
% % %         iddF = TLT(iPrj,12);
% % %         idPHI = TLT(iPrj,13);
% % % 
% % %         data = rWeight;
% % % 
% % %         % The central section is located at -1.*tiltAngle in the 3dFT since we
% % %         % rotate the specimen back by this amount in reconstruction.
% % %         rTilt = BH_defineMatrix([90,TLT(iPrj,4),-90],'Bah','invVector');
        if (calcAllWeights)
          rSubTomo = reshape(positionList(iSubTomo,17:25),3,3);
        else
          rSubTomo = eye(3);
        end
% % %         % To calculate the defocus, we need to rotate to where the projection is
% % %         % in real space TLT(iPrj,4)
% % % 
% % %         rProjection = BH_defineMatrix([90,TLT(iPrj,4),-90],'Bah','forwardVector');
% % %         prjCoords = rProjection*prjVector';
% % %         
% % %         iDf = prjCoords(3).*pixelSize + iDefocus;
% % %         defVect = [iDf - iddF, iDf + iddF, idPHI];
% % %         % Note the transpose (=inverse since the rotation matrix is orthogonal) is
% % %         % taken because the stored matrix is for interpolation
% % %         r = rSubTomo'*rTilt;
% % % 
% % %         % By default these are single PRECISION, but we truncate them anyhow, so
% % %         % leave as single.
% % %         [X,Y,Z,~,~,~] = BH_multi_gridCoordinates([SIZE(1:2),1],'Cartesian',METHOD,...
% % %                                   {'single',r,[0,0,0]','invVector',1,1},0,1,0);
% % % 
% % % 
% % %          % assuming ampContrast = 0.1, using -0.15 results in a weight with 
% % %         % (0.1^0.15)^2~ 0.5 at zero freqency. Allows some recovery of low freq without
% % %         % creating too severe a blur
% % %         [Hqz, ~] = BH_ctfCalc(radialCTFCalc,iCs,iWavelength,defVect,SIZE(1:2),iPhaseShift,-1.0); 
% % %         
% % % 
% % %         
% % % %          %Default is double
% % % %          if strcmpi(PRECISION,'single')
% % % %            Hqz = single(abs(Hqz.*HqzUnMod)).^ctfScaleFactor;
% % % %          elseif strcmpi(PRECISION, 'double')
% % % %            Hqz = double(abs(Hqz.*HqzUnMod)).^ctfScaleFactor;
% % % %          end
% % %        
% % %          %Default is double
% % %          if strcmpi(PRECISION,'single')
% % %            Hqz = single(abs(Hqz).^2);
% % %          elseif strcmpi(PRECISION, 'double')
% % %            Hqz = double(abs(Hqz).^2);
% % %          end
% % %          if SAMPLING > 1
% % %            Hqz = Hqz .^ (1*(SAMPLING - 1)^-3);           
% % %          end

         
% % %         fractionOfDose = TLT(iPrj,14)/mean(TLT(:,14));
% % %         fractionOfElastics = exp(-1.*iThickness/( cosd(TLT(iPrj,4))*400 ));
%     
%         fprintf('fractionOfDose %2.2f fractionOfElastic %2.2f at angle %2.2f\n',...
%                  fractionOfDose, fractionOfElastics,TLT(iPrj,4));

        % TODO does this make sense to run as power of one or two?
        exposureFilterPower = 1;
% % %         data = data.*Hqz.*exposureFilter(:,:,TLT(iPrj,1)).^exposureFilterPower.*(fractionOfDose*fractionOfElastics);


% % %         data = data(:);
% % % 
% % %         % Shift from image to array coordinates and set any out of bounds values to
% % %         % the origin.
% % %         originXYZ = ceil(((SIZE + 1)./2));
% % %         X = X(:)+originXYZ(1); 
% % %         Y = Y(:)+originXYZ(2); 
% % %         Z = Z(:)+originXYZ(3); 
% % % 
% % % 
% % %         outOfBounds = logical(( X < 1 | X > SIZE(1) ) + (Y < 1 | Y > SIZE(2)) + (Z < 1 | Z > SIZE(3)));
% % %         X(outOfBounds) = 1;%originXYZ(1);
% % %         Y(outOfBounds) = 1;%originXYZ(2);
% % %         Z(outOfBounds) = 1;%originXYZ(3);
% % % 
% % % 
% % % 
% % %         l = sub2ind(SIZE,round(X),round(Y),round(Z));
% % %         clear X Y Z
% % %         % much faster and effecitive for this simple sort.
% % %         [B,I] = sort(l);
% % %         idu = I(logical(B(1:end-1) - B(2:end)));
% % % 
% % %         tiltScale = 1;%- ( abs(sind(TLT(iPrj,4))).*0.2 );
% % %         rec{iGold}(l(idu)) = rec{iGold}(l(idu)) + tiltScale.*data(idu);
% % %         clear l idu
% % %   
% % %       end % end loop over projections
      [maxSize,maxCoord] =  max(SIZE);
      if (2.*maxSize > 512)
        padScaling = maxSize ./ 512;
        paddedSize = 512.*[1,1,1];
      else
        paddedSize = ceil(2.*maxSize).*[1,1,1];
        padScaling = 1/2;
      end

      [ SF3D ] = fourierCtfRecTex(paddedSize, [positionList(iSubTomo,1:10),...
                                                xyzSubTomo.*SAMPLING,...
                                                positionList(iSubTomo,14:16),...
                                                rSubTomo(1:9),positionList(iSubTomo,26)],...
                                                TLT, reconGeometry, originPrj, ...
                                                SAMPLING);
                                              


% % % 
% % %       % For now just recreate the Hermitian pair
% % %       nX = paddedSize(1);
% % %       nZ = paddedSize(3);
% % %       oX = floor(nX/2)+1;
% % %       isOdd = mod(nX,2);
      [ rec{iGold}  ] = BH_multi_makeHermitian(SF3D, paddedSize, padScaling);
      

      
% % %       try    
% % %         tmpArray = zeros(paddedSize,'single','gpuArray');
% % %         tmpArray(oX-1+isOdd:end,:,:) = SF3D;
% % %         tmpArray(1:oX-2+isOdd,:,:) = SF3D(oX-1+isOdd:-1:2,:,nZ:-1:1);
% % %         clear SF3D
% % % 
% % %         rec{iGold} =  BH_reScale3d(tmpArray,'',padScaling,'GPU');
% % %         clear tmpArray
% % %       catch
% % %         fprintf('SIZE %d paddedSize %d scaleFactor %d\n',SIZE(1),paddedSize(1),padScale);
% % %         error('I broke in generating the Hermitian mates!');
% % %       end
% % % %       rec{iGold}(oX-1+isOdd:end,:,:)  = SF3D;
% % % %       rec{iGold}(1:oX-2+isOdd,:,:) = SF3D(oX-1+isOdd:-1:2,:,nZ:-1:1);
% % % %       clear SF3D
      
      

%       % For now just recreate the Hermitian pair
%       nX = size(rec{iGold} ,1);
%       nZ = size(rec{iGold} ,3);
%       oX = floor(nX/2)+1;
%       isOdd = mod(nX,2);

%       rec{iGold}(oX-1+isOdd:end,:,:)  = SF3D;
%       rec{iGold}(1:oX-2+isOdd,:,:) = SF3D(oX-1+isOdd:-1:2,:,nZ:-1:1);
%       clear SF3D
      
          
% % %       for iGold = 1:1 + calcAllWeights
% % %         % Zero out the large value from out of bounds conditions
% % %         rec{iGold}(1) = 0;
% % %       end

    end % end loop over subTomos
  end %end loop over Tomos



% % %   g = BH_multi_gaussian3d(16.*[1,1,1],1.25);
% % %   gf=fftshift(fftn(ifftshift(g)));
% % %   g=real(fftshift(ifftn(ifftshift(gf.^3)))); clear gf
% % %   g = g ./ sum(g(:));
% % %   
% % %   if (useGPU)
% % %     g = gpuArray(g);
% % %   end
  

 %SAVE_IMG(MRCImage(gather(rec{1})),'tmpR1.mrc');
 %SAVE_IMG(MRCImage(gather(rec{2})),'tmpR2.mrc');

  for iGold = 1:1+calcAllWeights
      
    if (outputScaling)
      ctfWeights{iGold,iCtfGroup} = gather(BH_reScale3d( rec{iGold}...
                                      ,'', sprintf('%f',outputScaling), METHOD)); rec{iGold} = [];
    else
      ctfWeights{iGold,iCtfGroup} = gather(rec{iGold}); rec{iGold} = [];
    end
    


% % % % %     ctfWeights{iGold,iCtfGroup} = convn(single(rec{iGold}),g,'same'); rec{iGold} = [];

    
  
% % % %     if (outputScaling ~= 1)
% % % %       ctfWeights{iGold,iCtfGroup} = BH_reScale3d( ctfWeights{iGold,iCtfGroup}...
% % % %                                       ,'', sprintf('%f',outputScaling), METHOD);
% % % %     end
% % % %     if (flgFirstPass)
% % % %       % Only calc this once.
% % % %       rad = fftshift(BH_bandpass3d(size(ctfWeights{iGold,iCtfGroup}),...
% % % %                                                        0,0,0,'cpu','nyquist'));
% % % %       flgFirstPass = 0;                                               
% % % %     end
% % % %       
% % % %     ctfWeights{iGold,iCtfGroup} = gather(ctfWeights{iGold,iCtfGroup}) .* rad;
% % % %     
% % % % 
% % % %     ctfWeights{iGold,iCtfGroup} = (ctfWeights{iGold,iCtfGroup} - ...
% % % %                           min(ctfWeights{iGold,iCtfGroup}(rad > 0.0)))./((nPrjs)./2-0.5);
% % % %   
% % % %   % The mean value should probably be < 0.5 due to exposure filtering and
% % % %   % inelastic losses. Figure out how to put this on an absolute scale.
% % % %   %FIXME
% % % %   
% % % %   m = BH_movingAverage(ctfWeights{iGold,iCtfGroup},[7,7,7]);
% % % % 
% % % %   ctfWeights{iGold,iCtfGroup}(m < (mean(m(m>0)) + 0.5*std(m(m>0)))) = 0;
% % % %   
% % % %   meanPositiveValues = mean(mean(mean(ctfWeights{iGold,iCtfGroup}(ctfWeights{iGold,iCtfGroup} > 1e-2))));
% % % %   
% % % %   ctfWeights{iGold,iCtfGroup} = gather(ctfWeights{iGold,iCtfGroup} .* (0.5/meanPositiveValues));


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

