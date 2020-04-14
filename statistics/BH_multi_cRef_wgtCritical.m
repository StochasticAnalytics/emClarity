function [sf3d,wienerThreshold] = BH_multi_cRef_wgtCritical(sf3d)
%Downweight critically undersampled regions when pref. orientation is a
%problem.
%   Detailed explanation goes here

  startingMax = max(sf3d(:));
  valAtZero = max(10,0.1*startingMax); % ~ value at zero sampling (a bit less after the subtraction to keep the
                  % value at minNumSampled unchanged with a smooth transition.
                  
  wienerThreshold = (1.5.*(median(sf3d(sf3d(:)>valAtZero))-valAtZero));



  minNumSampled = 0.2.*median(sf3d(sf3d(:)>10)); % value below where a penalty is add (very little until low numbers)
  minFactor = 75/minNumSampled;
  minWeight = gpuArray(10); % decreasing this increase the downweighting as you move from 0 to minNumSampled

    % There shouldb't be any less than zero but due to the quality weighting
  % there could be
  sf3d(sf3d < 1) = 1;
  m = (sf3d < minNumSampled);

  sf3d(m) = (sf3d(m)+1) + valAtZero.^(minWeight.^((minFactor.*sf3d(m)+1).^-1));
  sf3d(m) = sf3d(m)-valAtZero.^(minWeight.^(minFactor*minNumSampled+1).^-1)+1;
  m = sf3d > startingMax;
  sf3d(m) = startingMax + log(sf3d(m));
  

      
end

