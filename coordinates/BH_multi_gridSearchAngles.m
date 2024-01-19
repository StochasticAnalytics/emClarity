function [  nIN_PLANE, IN_PLANE_SEARCH, angleStep, nAngles ] = ...
  BH_multi_gridSearchAngles( ANGLE_SEARCH)
%Consolodating function, calculate angular sampling.
%
%
%   Called by:
%
%   BH_templateSearch3d -
%     doesn't use: nAngles, ANGLE_LIST
%
%   BH_alignClass3d -
%     doesn't use: nIN_PLANE, ANGLE_INCREMENT
%
%   BH_alignRaw3d -
%     doesn't use: nIN_PLANE, ANGLE_INCREMENT
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Goals & limitations:
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   TODO:
%     - handle helical grid search
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

OUT_OF_PLANE = ANGLE_SEARCH(1:2);
IN_PLANE = ANGLE_SEARCH(3:4);

% For monolayer which can generally be at 0 or 180 +/- create search option
% around these two, indicated by a negative value for the out of plane range

if (OUT_OF_PLANE(1) < 0)
  OUT_OF_PLANE(1)= abs(OUT_OF_PLANE(1));
  biPolarSearch = 1;
else
  biPolarSearch = 0;
end


if (length(ANGLE_SEARCH) == 5)
  symmetryConstrainedSearch = ANGLE_SEARCH(5);
else
  symmetryConstrainedSearch = 0;
end


if all(IN_PLANE)
  
  
  psiStep = IN_PLANE(2);
  IN_PLANE_SEARCH = -IN_PLANE(1):psiStep:IN_PLANE(1);
  %psiStep = double(psiStep); IN_PLANE_SEARCH = double(IN_PLANE_SEARCH);
  if ~ismember(0,IN_PLANE_SEARCH)
    IN_PLANE_SEARCH(end+1) = 0;
  end
  
  if (symmetryConstrainedSearch)
    for iSym = 1:symmetryConstrainedSearch-1
      IN_PLANE_SEARCH = [IN_PLANE_SEARCH,IN_PLANE_SEARCH + iSym.*(360/symmetryConstrainedSearch)];
    end
  end
  nIN_PLANE = length(IN_PLANE_SEARCH);
  angleStep(1,:) = [0,0,0,0,psiStep];
else
  IN_PLANE_SEARCH = 0;
  if (symmetryConstrainedSearch)
    for iSym = 1:symmetryConstrainedSearch-1
      IN_PLANE_SEARCH = [IN_PLANE_SEARCH,IN_PLANE_SEARCH + iSym.*(360/symmetryConstrainedSearch)];
    end
    
  end
  
  psiStep = 0;
  nIN_PLANE = 1; % Always search the unrotated sample
  angleStep(1,:) = [0,0,0,0,0];
end


if all(OUT_OF_PLANE)
  
  if rem(OUT_OF_PLANE(1), OUT_OF_PLANE(2))
    error('Out of plane range not divisible by out of plane step.')
  else
    topPolar = OUT_OF_PLANE(1)/OUT_OF_PLANE(2);
    thetaStep = OUT_OF_PLANE(2);
  end
  
  if (biPolarSearch)
    polarAngles = 0:OUT_OF_PLANE(2):OUT_OF_PLANE(1);%(0:topPolar).*thetaStep;
    polarAngles = [polarAngles, flip(180-polarAngles)];
  else
    polarAngles = 0:OUT_OF_PLANE(2):OUT_OF_PLANE(1);%(0:topPolar).*thetaStep
  end
  
  for iPolarAngle = 1:length(polarAngles)
    theta = polarAngles(iPolarAngle);
    
    if ( theta == 0 || theta == 180)
      nAzimuthal = 0;
      phiStep = 0.5.* sind(thetaStep)^-1.0 .* thetaStep;
    else
      phiStep = sind(theta/1)^-1.0 .* thetaStep;
      if (theta - thetaStep == 0) || (theta + thetaStep == 180)
        % strict even spacing leaves the first out of plane undersampled
        phiStep = 0.5 * phiStep;
        
      end
      
      
      nAzimuthal = floor(360/phiStep);
      
      
    end
    % first position is psiStep independent of this.
    if (iPolarAngle)
      angleStep(iPolarAngle,:) = [theta, nAzimuthal, ...
        phiStep, thetaStep, psiStep];
    else
      angleStep(iPolarAngle,:) = [theta, nAzimuthal, ...
        phiStep, 0, psiStep];
    end
  end
else
  % Setting top polar limits the angular search to at most the in plane angles,
  % as the azimuth is also zero for iPolarAngle = 1.
  topPolar = 0;
end

% number of angles in the
% full search, rotationally averaged, first in plane, refinement1, refinement 2.
nAngles = zeros(5,1);
if sum(any(angleStep))
  
  for i = 1:size(angleStep,1)
    if angleStep(i,1) == 0
      
      % theta = 0 in plane only
      nAngles(1) = nAngles(1) + nIN_PLANE;
    else
      
      nAngles(1) = nAngles(1) + (angleStep(i,2) .*  nIN_PLANE);
    end
    nAngles(2) = nAngles(2) +  angleStep(i,2);
  end
else
  % Translational only search
  nAngles = nAngles + 1;
end

nAngles(3) = 10 .* nIN_PLANE;
nAngles(4) = 175;
nAngles(5) = 630;




end % end of gridSearchAngles function

