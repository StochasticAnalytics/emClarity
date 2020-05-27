function [ PEAK_COORD, mip ] =  BH_multi_xcf_Translational_2(  rotPART_FT, REF_FT, ...
                                                        wdgMask, refWdg,...
                                                        peakMask, PEAK_COM, ...                                                        
                                                        mip)
                                                       
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


% Option for different correlation functions, for now switching from Xcf to Mcf
% fixed, but only here (not rotational

try
iCCCmap = fftshift(real(ifftn( ...
            (refWdg.*rotPART_FT .* REF_FT .* wdgMask) ./ ...
            (sqrt( ( sum( abs(REF_FT   .*wdgMask).^2,'all'   ) .* ...
                     sum( abs(rotPART_FT.*refWdg).^2,'all')))./numel(wdgMask)))));
catch
  save('xcf_probs.mat');
  error('asdf')
end
               

if isnumeric(mip.x)
  updateVals = mip.x < iCCCmap(mip.mask);
  mip.x(updateVals) = iCCCmap(updateVals);
else
  mip.x = iCCCmap(mip.mask);
end
mip.N = mip.N + 1;

% The center of mass calc doesn't make sense when there are negative mass
% values, so originally, I set all less than 0 = 0. This ignores plently of
% useful information, so instead, now shift all intensitys to be >= 0
% % % iCCCmap(iCCCmap < 0) = 0;
iCCCmap = iCCCmap - min(iCCCmap(:));

if ~isempty(eraseMask)
  % Zero out tighter zone for cases of repeating lattice
  iCCCmap(eraseMask) = 0;
end

peakCOM = PEAK_COM;


[maxVal, maxInd] = max(iCCCmap(:));
[maxX, maxY, maxZ] = ind2sub(size(iCCCmap),maxInd);
PEAK_COORD = [maxX, maxY, maxZ];

% Sometimes there is an "empty class", ignore this in xfc
if (~isnan(gather(maxVal)) && (maxVal ~= 0)) 

  for iCOM = 1:2
    
    if (peakCOM)
      if iCOM == 2

  %         figure, imshow3D(gather(iCCCmap)); pause(20)


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
      catch
        peakCOM = gather(peakCOM);
        peakCoord = gather(peakCoord);
        peakLOW = gather(peakLOW);
        peakTOP = gather(peakTOP);
        iCCCmap = gather(iCCCmap);
        peakMask = gather(peakMask);
        save('xfcCalcErr.mat','peakCOM','peakCoord','peakLOW','peakTOP',...
                              'iCCCmap','peakMask');
        error('failed to box out COM calc, saving troubleshooting variables in xfcClacErr.mat')
      end


    cMass = [ sum(sum(sum(boX.*cmX))) ; ... 
              sum(sum(sum(boX.*cmY))) ; ...
              sum(sum(sum(boX.*cmZ))) ] ./ sum(boX(:));
            
    PEAK_COORD = peakCoord + cMass';

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



clearvars -except PEAK_COORD mip

end % end of multi_xcf

