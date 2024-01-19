function [ CCC_STORAGE] = BH_multi_angularSearch( ANGLE_STEP, ...
  PEAK_LIST, ...
  IN_PLANE_SEARCH, ...
  iClassImg, iClassWdg, ...
  ref_FT, refWDG, ...
  refRotAvg_FT, ...
  volMask, bandpassFilt, ...
  padCalc, padREF,...
  peakMask, peakCOM, IDX, ...
  refSym)
%UNTITLED Summary of this function goes her
%   Detailed explanation goes here

% Make a normalization factor for wedge weighting. This could be done in the
% begining, but to test I'll put it here so I don't have to change functino I/o

peakBinary = (peakMask >= 0.01);
volBinary  = (volMask >= 0.01);



% If searching the full in plane range, limit based on symmetry. Important for
% wedge bias and also helps with speed.
if IN_PLANE_SEARCH(1) == -180
  limitSymmetry = 1;
  fprintf('limiting to symmetry constrained in-plane search.\n')
else
  limitSymmetry = 0;
end



iClassTrim = iClassImg(padREF(1,1)+1 : end - padREF(2,1), ...
  padREF(1,2)+1 : end - padREF(2,2), ...
  padREF(1,3)+1 : end - padREF(2,3) );


%[ iClassImg ] = BH_bandLimitCenterNormalize(unMaskedClassImage, bandpassFilt, volMask);



if length(PEAK_LIST) > 1
  % Search angles based on previously found peaks.
  peakList = PEAK_LIST;
  
  % Either a refinment around top three or top one peaks, or checking all
  % inplane angles for top ten axially averaged.
  if (size(peakList,2) == 10)
    % Refinementd
    % Get unique references for peakList wedge weight normalization
    referenceList = unique(peakList(:,1));
    nRefs = length(referenceList);
    iClassImg2 = cell(nRefs,1);
    for iRef = 1:length(bandpassFilt)
      iRef
      size( iClassTrim)
      size(bandpassFilt{iRef})
      [ iClassImg2{iRef} ] = BH_bandLimitCenterNormalize(iClassTrim.*peakMask, bandpassFilt{iRef}, peakBinary,padCalc,'single');
    end
    angCount=1;
    
    nAngles = 46.*size(peakList,1)+1
    pause(4)
    cccStorage = zeros(nAngles, 13, 'double', 'gpuArray');
    
    for iAngle = 1 : size(peakList,1)
      
      
      iRef = peakList(iAngle, 1);
      phi  = peakList(iAngle, 2); phiInc   = peakList(iAngle,5);
      theta= peakList(iAngle, 3); thetaInc = peakList(iAngle,6);
      psi  = peakList(iAngle, 4); psiInc   = peakList(iAngle,7);
      
      if size(peakList,1) == 1 % This is the final refinement
        %         superSample = 1;
        inPlaneSearch = psi-psiInc : psiInc : psi+psiInc;
        polarSearch   = theta-thetaInc :thetaInc : theta+thetaInc;
        azimuthalSearch= phi-phiInc    : phiInc    : phi  + phiInc;
      else
        inPlaneSearch  = psi-psiInc : psiInc    : psi  + psiInc;
        polarSearch    = theta-thetaInc: thetaInc  : theta+ thetaInc;
        azimuthalSearch= phi-2*phiInc    : phiInc    : phi  + 2*phiInc;
      end
      
      
      
      
      
      searchList = zeros(46,3);
      nSearch = 1;
      for iPhi = azimuthalSearch
        for iTheta = polarSearch
          for iPsi = inPlaneSearch
            searchList(nSearch, :) = [iPhi, iTheta, iPsi];
            nSearch = nSearch + 1;
          end
        end
      end
      
      
      
      for iRefine = 1:nSearch-1
        
        
        RotMat = BH_defineMatrix(searchList(iRefine,:),'Bah', 'forward');
        
        
        [ rotRef ] = BH_resample3d(ref_FT(:,:,:,iRef),RotMat, ...
          peakList(iAngle,8:10), ...
          'Bah', 'GPU', 'forward');
        if isa(refWDG,'cell')
          rotWDG = ifftn(BH_resample3d(refWDG{iRef},RotMat, ...
            peakList(iAngle,8:10), ...
            'Bah', 'GPU', 'forward'));
        else
          rotWDG = refWDG;
        end
        
        rotRef = rotRef(padREF(1,1)+1 : end - padREF(2,1), ...
          padREF(1,2)+1 : end - padREF(2,2), ...
          padREF(1,3)+1 : end - padREF(2,3) );
        
        
        
        
        rotRef_FT =  ...
          BH_bandLimitCenterNormalize(rotRef.*volMask, bandpassFilt{iRef}, volBinary,padCalc,'single');
        rotRef_FT2 =  ...
          BH_bandLimitCenterNormalize(rotRef.*peakMask, bandpassFilt{iRef}, peakBinary,padCalc,'single');
        
        rotRef_FT = conj(rotRef_FT);
        rotRef_FT2= conj(rotRef_FT2);
        
        
        % find translational shift using rotationally averaged tightly masked
        % reference, to reduce chance of drift to alternate lattice sites.
        try
          [ estPeakCoord ] =  BH_multi_xcf_Translational(  ...
            iClassImg2{iRef}, rotRef_FT2, ...
            peakMask, peakCOM);
        catch
          iRef
          
          
          error('sdfsd')
        end
        
        % apply only a tranlational shift to the particle
        [ rotClassImg ]  = BH_resample3d(iClassImg, [0,0,0], ...
          estPeakCoord,'Bah','GPU','inv');
        
        rotClassImg = rotClassImg(padREF(1,1)+1 : end - padREF(2,1), ...
          padREF(1,2)+1 : end - padREF(2,2), ...
          padREF(1,3)+1 : end - padREF(2,3) );
        
        
        [ rotClassImg1 ]  = BH_bandLimitCenterNormalize(...
          rotClassImg.*volMask, ...
          bandpassFilt{iRef} , volBinary,...
          padCalc,'single');
        %         [ rotClassImg2 ]  = BH_bandLimitCenterNormalize(...
        %                                                    rotClassImg.*peakMask, ...
        %                                                    bandpassFilt, peakBinary, ...
        %                                                    padCalc,'single');
        clear rotClassImg
        
        % now calc CCC, setting sampling shift to zero
        [ iCCC, iWeight ] = ...
          BH_multi_xcf_Rotational( rotClassImg1, rotRef_FT, ...
          iClassWdg, rotWDG, ...
          peakMask);
        
        %         [ finalPeakCoord ] =  BH_multi_xcf_Translational(  ...
        %                                                    rotClassImg2, rotRef_FT2, ...
        %                                                    peakMask, peakCOM);
        cccStorage(angCount,:) = [iRef, IDX, searchList(iRefine,:), iCCC, ...
          iWeight, peakList(iAngle,8:10)+estPeakCoord, ...
          phiInc,thetaInc,psiInc];
        
        angCount = angCount + 1;
      end % end inPlane
    end % end search over best peaks
    
  else
    % top ten
    % Get unique references for peakList wedge weight normalization
    
    %     referenceList = unique(peakList(:,1));
    %     nRefs = length(referenceList);
    angCount=1;
    
    nRefs = size(ref_FT,4);
    iClassImg2 = cell(nRefs,1);
    for iRef = 1:nRefs
      [ iClassImg2{iRef} ] = BH_bandLimitCenterNormalize(iClassTrim.*peakMask, bandpassFilt{iRef}, peakBinary,padCalc,'single');
    end
    nAngles = nRefs.*size(peakList,1).*length(IN_PLANE_SEARCH);
    cccStorage = zeros(nAngles, 10, 'double', 'gpuArray');
    
    for iAngle = 1:size(peakList,1)
      % The assumption is the best reference for the axially averaged is
      % also the best ref otherwise. Maybe not true.
      %iRef = peakList(iAngle, 1);
      phi  = peakList(iAngle, 2);
      theta= peakList(iAngle, 3);
      
      for iInPlane = IN_PLANE_SEARCH
        psi = iInPlane;
        
        evaluateRef = ones(1,nRefs);
        if (limitSymmetry)
          evaluateRef = evaluateRef.*((abs(psi).*evaluateRef) < 180 ./ refSym);
        end
        for iRef = 1:nRefs
          if (evaluateRef(iRef))
            RotMat = BH_defineMatrix([phi, theta, psi - phi],'Bah', 'forward');
            
            [ rotRef ] = BH_resample3d(ref_FT(:,:,:,iRef), ...
              RotMat,peakList(iAngle,4:6), ...
              'Bah', 'GPU', 'forward');
            
            
            if isa(refWDG,'cell')
              rotWDG = ifftn(BH_resample3d(refWDG{iRef},RotMat, ...
                peakList(iAngle,4:6), ...
                'Bah', 'GPU', 'forward'));
            else
              rotWDG = refWDG;
            end
            
            rotRef = rotRef(padREF(1,1)+1 : end - padREF(2,1), ...
              padREF(1,2)+1 : end - padREF(2,2), ...
              padREF(1,3)+1 : end - padREF(2,3) );
            
            
            
            
            rotRef_FT =  ...
              BH_bandLimitCenterNormalize(rotRef.*volMask, bandpassFilt{iRef} , volBinary,padCalc,'single');
            rotRef_FT2 =  ...
              BH_bandLimitCenterNormalize(rotRef.*peakMask, bandpassFilt{iRef} , peakBinary,padCalc,'single');
            
            rotRef_FT = conj(rotRef_FT);
            rotRef_FT2= conj(rotRef_FT2);
            
            
            
            
            
            % find translational shift
            [ estPeakCoord ] =  BH_multi_xcf_Translational(  ...
              iClassImg2{iRef}, rotRef_FT2, ...
              peakMask, peakCOM);
            
            % apply only a tranlational shift to the particle
            [ rotClassImg ]  = BH_resample3d(iClassImg, [0,0,0], ...
              estPeakCoord,'Bah','GPU','inv');
            
            rotClassImg = rotClassImg(padREF(1,1)+1 : end - padREF(2,1), ...
              padREF(1,2)+1 : end - padREF(2,2), ...
              padREF(1,3)+1 : end - padREF(2,3) );
            
            
            [ rotClassImg1 ]  = BH_bandLimitCenterNormalize(rotClassImg.*volMask, ...
              bandpassFilt{iRef} , volBinary,padCalc,'single');
            %             [ rotClassImg2 ]  = BH_bandLimitCenterNormalize(rotClassImg.*peakMask, ...
            %                                                              bandpassFilt{iRef} , peakBinary,padCalc,'single');
            clear rotClassImg
            
            
            % now calc CCC, setting sampling shift to zero
            [ iCCC, iWeight ] = ...
              BH_multi_xcf_Rotational( rotClassImg1, rotRef_FT, ...
              iClassWdg, rotWDG,...
              peakMask);
            
            %             [ finalPeakCoord ] =  BH_multi_xcf_Translational(  ...
            %                                                        rotClassImg2, rotRef_FT2, ...
            %                                                        peakMask, peakCOM);
            
            
            cccStorage(angCount,:) = [iRef, IDX, phi, theta, psi - phi, iCCC, ...
              iWeight,  estPeakCoord + ...
              peakList(iAngle,4:6) ];
            angCount = angCount + 1;
          end
        end
      end % end inPlane
    end % end search over best peaks
  end
  
else % search angles based on grideSearchAngles
  
  referenceList = 1:size(ref_FT,4);
  nRefs = length(referenceList);
  iClassImg2 = cell(nRefs,1);
  
  
  for iRef = 1:nRefs
    
    [ iClassImg2{iRef} ] = BH_bandLimitCenterNormalize(iClassTrim.*peakMask, bandpassFilt{iRef}, peakBinary,padCalc,'single');
  end
  angCount=1;
  nAngles = nRefs.*sum(ANGLE_STEP(:,2)+1).*length(IN_PLANE_SEARCH);
  cccStorage = zeros(nAngles, 10, 'double', 'gpuArray');
  for iAngle = 1:size(ANGLE_STEP,1)
    
    theta = ANGLE_STEP(iAngle,1);
    
    % Calculate the increment in phi so that the azimuthal sampling is
    % consistent and equal to the out of plane increment.
    
    phiStep = ANGLE_STEP(iAngle,3);
    
    % To prevent only searching the same increments each time in a limited
    % grid search, radomly offset the azimuthal angle by a random number
    % between 0 and 1/2 the azimuthal increment.
    azimuthalRandomizer = rand(1)*phiStep/2;
    
    for iAzimuth = 0:ANGLE_STEP(iAngle,2)
      phi = rem(phiStep * (iAzimuth + azimuthalRandomizer),360) ;
      
      % For axially averaged this is always zero.
      for iInPlane = IN_PLANE_SEARCH
        psi = iInPlane;
        [phi,theta,psi];
        
        for iRef = 1:nRefs
          
          
          
          RotMat = BH_defineMatrix([phi, theta, psi ],'Bah', 'forward');
          
          
          [ rotRef ] = BH_resample3d(refRotAvg_FT(:,:,:,iRef),RotMat, ...
            [0,0,0], 'Bah', 'GPU', 'forward');
          
          
          
          if isa(refWDG,'cell')
            rotWDG = ifftn(BH_resample3d(refWDG{iRef},RotMat, ...
              [0,0,0], ...
              'Bah', 'GPU', 'forward'));
          else
            rotWDG = refWDG;
          end
          
          rotRef = rotRef(padREF(1,1)+1 : end - padREF(2,1), ...
            padREF(1,2)+1 : end - padREF(2,2), ...
            padREF(1,3)+1 : end - padREF(2,3) );
          
          
          try
            rotRef_FT =  ...
              BH_bandLimitCenterNormalize(rotRef.*volMask, bandpassFilt{iRef} , volBinary,padCalc,'single');
            rotRef_FT2 =  ...
              BH_bandLimitCenterNormalize(rotRef.*peakMask, bandpassFilt{iRef} , peakBinary,padCalc,'single');
          catch
            size(volMask)
            size(peakMask)
            size(bandpassFilt{iRef} )
            size(rotRef)
            error('size mismatch in multi_angular search')
          end
          rotRef_FT = conj(rotRef_FT);
          rotRef_FT2= conj(rotRef_FT2);
          
          
          % find translational shift
          [ estPeakCoord ] =  BH_multi_xcf_Translational(  ...
            iClassImg2{iRef}, rotRef_FT2, ...
            peakMask, peakCOM);
          
          % apply only a tranlational shift to the particle
          [ rotClassImg ]  = BH_resample3d(iClassImg, [0,0,0], ...
            estPeakCoord,'Bah','GPU','inv');
          
          rotClassImg = rotClassImg(padREF(1,1)+1 : end - padREF(2,1), ...
            padREF(1,2)+1 : end - padREF(2,2), ...
            padREF(1,3)+1 : end - padREF(2,3) );
          
          
          [ rotClassImg1 ]  = BH_bandLimitCenterNormalize(rotClassImg.*volMask, ...
            bandpassFilt{iRef} , volBinary,padCalc,'single');
          %         [ rotClassImg2 ]  = BH_bandLimitCenterNormalize(rotClassImg.*peakMask, ...
          %                                                          bandpassFilt{iRef} , peakBinary,padCalc,'single');
          clear rotClassImg
          
          
          % now calc CCC, setting sampling shift to zero
          [ iCCC, iWeight ] = ...
            BH_multi_xcf_Rotational( rotClassImg1, rotRef_FT, ...
            iClassWdg,rotWDG, ...
            peakMask);
          
          %         [ finalPeakCoord ] =  BH_multi_xcf_Translational(  ...
          %                                                    rotClassImg2, rotRef_FT2, ...
          %                                                    peakMask, peakCOM);
          %
          
          cccStorage(angCount,:) = [iRef, IDX, phi, theta, psi, iCCC, ...
            iWeight,  estPeakCoord];
          
          angCount = angCount + 1;
        end % peak search over refs
        
        
      end % end inPlane
    end % azimuthal
    
  end % polar
  
  % end of else clause, which is a search over axially averaged reference.
end

%[ cccStorage ] = BH_multi_peakSearch(referenceList, cccStorage);
cccStorage = cccStorage(( ~(sum(isnan(cccStorage),2)) ),:);
cccStorage = cccStorage(( cccStorage(:,6) ~= 0 ),:);

CCC_STORAGE = sortrows(gather(cccStorage),-6);

clear  rotRef rotWDG refWDG rotRef_FT rotRef_FT2 clear rotClassImg1 iClassImg2
end % end angularSearch function

