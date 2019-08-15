  function [ ctfDepth ] = BH_ctfCalcError( pixelSize, Cs, Lambda, Defocus, ...
                                            CTFSIZE, AMPCONT, ...
                                            resCutOff,thicknessAng,dampeningMax,cycleNumber)
%For a given resolution wanted, calculate the point where destructive
%interference drops the CTF amplitude to some given threshold.
%   Detailed explanation goes here

% dimensions are the wanted resolution in Ang, the tomogram thickness in
% pixels, and the cutoff for allowable dampening.
NYQ = ceil((CTFSIZE+1)/2);
pixelSize = pixelSize*10^10;

if resCutOff <= 2*pixelSize
  resCutOff = 2.55*pixelSize
end

[ rad ] = BH_multi_gridCoordinates([CTFSIZE,1],'Cylindrical','GPU',{'none'},1,0,1);
rad = {rad./pixelSize,0,rad.*0};

% ctf1 = BH_ctfCalc(rad,Cs,Lambda,Defocus,CTFSIZE,AMPCONT,-1);

% The max difference from the mean defocus is half the thickness
maxDiff = thicknessAng./2;

% If you just add the beginning and the end, you can have a second lobe (of
% a full phase reversal) at high res for low (1ish) defocus. Add together a
% range to avoid this.
nCtf = 0;
ctf1 = [];
length([-maxDiff:0.1*maxDiff:maxDiff].*10^-10)
for jDelDef = [-maxDiff:0.1*maxDiff:maxDiff].*10^-10
  if isempty(ctf1)
    ctf1 = BH_ctfCalc(rad,Cs,Lambda,Defocus+jDelDef,CTFSIZE,AMPCONT,-1,1);    
  else
    ctf1 = ctf1 + BH_ctfCalc(rad,Cs,Lambda,Defocus+jDelDef,CTFSIZE,AMPCONT,-1,1);
  end
  nCtf = nCtf + 1;
end


% Working along first dimension, so logical indexing ignoring 2d is okay.

noErrorIDX = find( abs(ctf1(1:NYQ)./nCtf) > dampeningMax,1,'last');
noErrorRes = gather(1 / rad{1}(noErrorIDX));

% c2 = BH_ctfCalc(rad,Cs,Lambda,Defocus-(maxDiff*10^-10),CTFSIZE,AMPCONT,-1);
% c3 = abs((ctf1(1:NYQ)+c2(1:NYQ))./2);
% lastIDX = find(c3 > dampeningMax,1,'last')
% highestRes = 1 / rad{1}(lastIDX)

ctfDepth = -1;

if noErrorRes < resCutOff
  ctfDepth = 2*maxDiff*10^-10;
  fprintf('\nat this pixelSize %3.3f the mean defocus is sufficient\n',pixelSize);
else
  % Hard set the minimum ctfDepth to be 10nm - based on other published
  % results, this may be too conservative
  while maxDiff > 5
    maxDiff = maxDiff .* .9;
    nCtf = 1;
    ctf1 = [];
    for jDelDef = [-maxDiff:0.1*maxDiff:maxDiff].*10^-10
      if isempty(ctf1)
        ctf1 = BH_ctfCalc(rad,Cs,Lambda,Defocus+jDelDef,CTFSIZE,AMPCONT,-1);    
      else
        ctf1 = ctf1 + BH_ctfCalc(rad,Cs,Lambda,Defocus+jDelDef,CTFSIZE,AMPCONT,-1);
      end
      nCtf = nCtf + 1;
    end
    % Working along first dimension, so logical indexing ignoring 2d is okay.
    noErrorIDX = find(abs(ctf1(1:NYQ)./nCtf) > dampeningMax,1,'last');
    noErrorRes = 1 / rad{1}(noErrorIDX);
    
    if any(noErrorRes < resCutOff)
      ctfDepth = 2*maxDiff*10^-10;
      break
    end
  end
end


if (cycleNumber)
  if ctfDepth < 0
    fprintf('The optimal ctfDepth was not found\n');
    fprintf('Inputs %3.3e pix %3.3e cs %3.3e wl %3.3e def %3.3e resTarget %3.3e tomoDepth\n',...
             pixelSize*10^-10, Cs, Lambda, Defocus,resCutOff,thicknessAng);
    ctfDepth = min(thicknessAng/30 * 10^-9,resCutOff(1) * 10 ^-8)
  elseif ctfDepth < 0.5*10e-9
    fprintf('\n\nCapping ctfDepth to 5 nm from a calc %3.3f nm\n\n',ctfDepth*10^9);
    ctfDepth = 10e-9;
  end
else
  % Cap cycle zero to a max of 3 tilts
  if ctfDepth < thicknessAng*10^-10/3 || ctfDepth == -1
    fprintf('Cycle 0, cap ctfDepth from %3.3e nm  to %3.3e nm \n',ctfDepth*10^9,thicknessAng/30);
    ctfDepth = (thicknessAng/30)*10^-9;
  end
end
