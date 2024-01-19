function [ iCCC, iWeight ] =  BH_multi_xcf_Rotational( rotPART_FT, REF_FT, ...
  wdgMask, refWdg, ...
  wCCC)

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

iWeight = 1;

% Testing, might be better later to just make a separate logical peak mask
% once at the begining and pass it in. For now this means far fewer
% changes.
if isnan(gather(refWdg(1)))
  refWdg = 1;
end

if isnan(gather(wdgMask(1)))
  wdgMask = 1;
end

if isa(wCCC, 'cell')
  flgWeightCCC = 1;
else
  flgWeightCCC = 0;
end



if (flgWeightCCC)
  iCCC = 0;
  
  for iShell = 3:length(wCCC);
    
    % reference is the conjugate of filtered fft
    numerator = real(refWdg(wCCC{iShell}).*rotPART_FT(wCCC{iShell}) .* ...
      REF_FT(wCCC{iShell}) .* wdgMask(wCCC{iShell}));
    denominator = ( sum(sum(sum( abs(REF_FT(wCCC{iShell}).* ...
      wdgMask(wCCC{iShell})).^2   ))) .* ...
      sum(sum(sum( abs(rotPART_FT(wCCC{iShell}).*...
      refWdg(wCCC{iShell})).^2))) );
    
    
    iCCCshell    = real(sum(sum(sum(numerator))) ./ sqrt(denominator));
    if wCCC{1}
      iCCC = iCCC + abs((iCCCshell^2 / (iCCCshell + wCCC{2}^-1.5))^3);
    else
      iCCC = iCCC + (iCCCshell^2 / (iCCCshell + wCCC{2}^-1.5))^3;
    end
  end
else
  % reference is the conjugate of filtered fft
  numerator = real(refWdg.*rotPART_FT .* REF_FT .* wdgMask);
  denominator = ( sum(sum(sum( abs(REF_FT.*wdgMask).^2   ))) .* ...
    sum(sum(sum( abs(rotPART_FT.*refWdg).^2))) );
  
  
  iCCC    = sum(sum(sum(numerator))) ./ sqrt(denominator);
  
end
clear rotPART_FT REF_FT  REF_WDG iWedgeMask peakMask numerator denominator

end % end of multi_xcf

