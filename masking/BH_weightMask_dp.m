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

    iThickness = 75;

   

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

        if (calcAllWeights)
          rSubTomo = reshape(positionList(iSubTomo,17:25),3,3);
        else
          rSubTomo = eye(3);
        end


        % TODO does this make sense to run as power of one or two?
        exposureFilterPower = 1;

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

      [ rec{iGold}  ] = BH_multi_makeHermitian(SF3D, paddedSize, padScaling);

      
      



    end % end loop over subTomos
  end %end loop over Tomos





  for iGold = 1:1+calcAllWeights
      
    if (outputScaling)
      ctfWeights{iGold,iCtfGroup} = gather(BH_reScale3d( rec{iGold}...
                                      ,'', sprintf('%f',outputScaling), METHOD)); rec{iGold} = [];
    else
      ctfWeights{iGold,iCtfGroup} = gather(rec{iGold}); rec{iGold} = [];
    end
    




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

