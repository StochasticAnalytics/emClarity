function [ subTomoMeta ] = BH_recordAngularSampling( subTomoMeta,cycleNumber, angleStep, inPlaneSearch)
% Record the angles sampled
%   This is just the distribution, but the angles for every subtomo are going to
%   be different due to the azimuthalRandomizer.

nAngles = numel(inPlaneSearch) * sum(angleStep(:,2) + 1)

anglesSampled = zeros(nAngles,3);
nANG = 1;

for iAngle = 1:size(angleStep,1)
  
  theta    = angleStep(iAngle,1);
  % % %   thetaInc = angleStep(iAngle,4);
  % Calculate the increment in phi so that the azimuthal sampling is
  % consistent and equal to the out of plane increment.
  
  phiInc = angleStep(iAngle,3);
  
  
  % To prevent only searching the same increments each time in a limited
  % grid search, radomly offset the azimuthal angle by a random number
  % between 0 and 1/2 the azimuthal increment.
  
  azimuthalRandomizer = (rand(1)-0.5)*phiInc;
  
  for iAzimuth = 0:angleStep(iAngle,2)
    phi = rem((phiInc * iAzimuth)+azimuthalRandomizer,360);
    
    for iInPlane = inPlaneSearch
      psi    = iInPlane;
      % % %     psiInc = angleStep(iAngle,5);
      
      try
        anglesSampled(nANG,:) = [phi,theta,psi -phi];
        nANG = nANG + 1;
      catch
      end
    end
  end
end

subTomoMeta.(cycleNumber).('anglesSampled') = anglesSampled;

end

