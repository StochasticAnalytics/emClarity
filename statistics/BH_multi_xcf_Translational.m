function [ PEAK_COORD ] =  BH_multi_xcf_Translational(  rotPART_FT, REF_FT, ...
  peakMask, PEAK_COM)

%Consolodating function, calculate wedge weight and cross correlation
%
%
%   Called by:
%
%   BH_alignRaw3d
%
%   BH_alignClass3d
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
%   TODO: Add something to report when no valid peak is found.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% regular cross-correlation
iCCCmap = fftshift(real(ifftn(rotPART_FT.*REF_FT)));

minVal = min(iCCCmap,[],'all');
iCCCmap(peakMask < 0.95) = minVal;


peakCOM = PEAK_COM;


[maxVal, maxInd] = max(iCCCmap(:));
[maxX, maxY, maxZ] = ind2sub(size(iCCCmap),maxInd);
PEAK_COORD = [maxX, maxY, maxZ];

% Sometimes there is an "empty class", ignore this in xfc
if (~isnan(gather(maxVal)) && (maxVal ~= 0))
  
  for iCOM = 1:2
    
    if (peakCOM)
      if iCOM == 2
        peakCoord = round(PEAK_COORD);
      else
        % from max within peak window
        peakCoord = PEAK_COORD;
      end
      
      peakLOW = peakCoord - peakCOM;
      peakTOP = peakCoord + peakCOM;
      
      [cmX, cmY, cmZ] = ndgrid(gpuArray(-1*peakCOM(1):peakCOM(1)), ...
        gpuArray(-1*peakCOM(2):peakCOM(2)), ...
        gpuArray(-1*peakCOM(3):peakCOM(3)) );
      
      try
        boX = iCCCmap(peakLOW(1):peakTOP(1), ...
          peakLOW(2):peakTOP(2), ...
          peakLOW(3):peakTOP(3));
        
        % Center of mass calc only makes sense for positive values
        boX = boX - minVal;

              
      cMass = [ sum(sum(sum(boX.*cmX))) ; ...
        sum(sum(sum(boX.*cmY))) ; ...
        sum(sum(sum(boX.*cmZ))) ] ./ sum(boX(:));
      
      PEAK_COORD = peakCoord + cMass';
      catch
        PEAK_COORD = peakCoord;
        % peakCOM = gather(peakCOM);
        % peakCoord = gather(peakCoord);
        % peakLOW = gather(peakLOW);
        % peakTOP = gather(peakTOP);
        % iCCCmap = gather(iCCCmap);
        % peakMask = gather(peakMask);
        % save('xfcCalcErr.mat','peakCOM','peakCoord','peakLOW','peakTOP',...
        %   'iCCCmap','peakMask');
        % error('failed to box out COM calc, saving troubleshooting variables in xfcClacErr.mat')
      end
      

      
    end
  end
else
  fprintf('maxVal in iCCCmap is %f\n', maxVal);
  
  PEAK_COORD = [0,0,0];
end
PEAK_COORD = PEAK_COORD - ceil((size(iCCCmap)+1)./2);

if any(isnan(PEAK_COORD))
  fprintf('%f %f %f\n', PEAK_COORD);
  fprintf('Setting nans to 0\n');
  PEAK_COORD =  PEAK_COORD.* ~isnan(PEAK_COORD);
end



clearvars -except PEAK_COORD

end % end of multi_xcf

