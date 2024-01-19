function [ ctfMask, Hqz] = BH_ctfCalc(PixelSize, Cs, Lamda, Defocus, ...
  CTFSIZE, AMPCONT, Phase_Only, varargin)
% Calculate a ctf and with the given modification to information at resolution
% lower than the first peak. Also return the unmodified ctf.

% FIXME: I think the need for double was for overflow in the phase
% we use single in cisTEM, so it must just be  scaling issue. Or better yet, switch over to mexCTF

precision = 'single';
flgComplex = 0;
calcOneD = 0;
thisZero = 0;
doHalfGrid = 0;
if nargin == 8
  if isnumeric(varargin{1})
    if varargin{1} < 0
      precision = 'double';
      calcOneD = 0;
    else
      maxZ = varargin{1};
      precision = 'double';
      calcOneD = 1;
    end
  else
    precision = 'double';
    calcOneD = 0;
    flgComplex = 1;
  end
elseif nargin ==9
  thisZero = varargin{2};
else
  calcOneD = 0;
end


if AMPCONT < 0
  AMPCONT = abs(AMPCONT);
  % This is used to damp higher frequecies in the CTF refinement during
  % tomoCPR. This is probably not a good way to do this and should be
  % removed. FIXME
  flgDampen = 1;
else
  flgDampen = 0;
end

if isa(PixelSize,'cell')
  radialGrid = PixelSize{1};
  preShiftedOrigin = PixelSize{2};
  phi = PixelSize{3};
  calcRad = 0;
  
  
  % While switching to the default behavior of using ang, check any input radial grid
  % and make sure it is in Angstrom. Check the middle of the frequency
  % range that way the input could be shifted or not.
  if (radialGrid(ceil((size(radialGrid,1)+1)/4),1) > 1e6)
    radialGrid = radialGrid .* 10^-10;
  end
  
  if numel(preShiftedOrigin) == 2
    doHalfGrid = preShiftedOrigin(2);
    preShiftedOrigin = preShiftedOrigin(1);
  end
else
  PIXEL_SIZE = PixelSize*10^10;
  preShiftedOrigin = 0;
  calcRad = 1;
end

CS = Cs;
WL = Lamda;

if numel(Defocus) == 1
  
  df1 = Defocus;
  df2 = Defocus;
  phi0 = 0;
elseif numel(Defocus) == 3
  df1 = Defocus(1);
  df2 = Defocus(2);
  phi0= Defocus(3);
end



CS = CS * 10^10;
WL = WL * 10^10;
df1 = df1 * 10^10;
df2 = df2 * 10^10;




if numel(CTFSIZE) == 1
  CTFSIZE(2) = CTFSIZE(1);
end

if (calcOneD)
  CTFSIZE(2) = 1;
end

if  ( calcRad )
  
  if strcmpi(precision, 'single')
    [radialGrid,phi,~,~,~,~] = ...
      BH_multi_gridCoordinates(CTFSIZE(1:2),'Cylindrical','GPU',{'none'},1,0,0);
  else
    
    [radialGrid,phi,~,~,~,~] = ...
      BH_multi_gridCoordinates(CTFSIZE(1:2),'Cylindrical','GPU',{'none'},1,0,0);
    radialGrid = double(radialGrid);
    phi = double(phi);
  end
  
  radialGrid = radialGrid ./ PIXEL_SIZE;
  
end

% Any additional phase shift due to the phase plate is stored with the
% amplitude contrast. This will produce very small errors (< 1% for .07,0.1)
% for older versions that expect just the amplitude contrast ratio, rather
% than the phase shift.
% % % if (abs(AMPCONT - 1.0) < 1e-3)
% % %    precomputed_amplitude_contrast_term = pi / 2.0;
% % % else
% % %   precomputed_amplitude_contrast_term = atan2(AMPCONT,sqrt(1.0 - AMPCONT^2));
% % % end

% df1 should be defocus of greater mag and phi0 -90/90

% phasePerturbation = pi.*(0.5.*CS.*WL^3.*(radialGrid).^4 + DF.*WL.*(radialGrid).^2);
dfTerm = 0.5.*( (df1+df2) + (df1-df2)*cos(2.*(phi-phi0)) );
phasePerturbation = pi.*(0.5.*CS.*WL^3.*(radialGrid).^4 + ...
  WL.*(radialGrid).^2 .* dfTerm);
% dPdQ = 2*pi*CS*WL^3.*radialGrid.^3 + 2*WL.*radialGrid.*dfTerm;
if ( flgComplex )
  ctfMask =  exp(-1i.*(phasePerturbation-atan2(AMPCONT,sqrt(1+AMPCONT))));
  Hqz = exp(+1i.*(phasePerturbation-atan2(AMPCONT,sqrt(1+AMPCONT))));
  return
else
  %   Hqz = (sqrt(1-AMPCONT^2).*sin(phasePerturbation) - AMPCONT.*cos(phasePerturbation));
end

Hqz = sin(phasePerturbation - AMPCONT);

nanCheck = isnan(Hqz);
if gather(sum(nanCheck(:)))
  Hqz(nanCheck) = 0;
end


if Phase_Only < 0
  
  
  % FIXME I need to know if I am a half grid or els this fails!
  
  
  oX = floor(CTFSIZE(1)/2)+1;
  oY = floor(CTFSIZE(2)/2)+1;
  
  
  if (calcOneD)
    if (preShiftedOrigin)
      rV = Hqz(oX:end);
    else
      rV = Hqz(1:oX); % should this be oX-1?
    end
  elseif (doHalfGrid)
    if ( preShiftedOrigin)
      rV = Hqz(1:end,oY);
    else
      rV = Hqz(1:end,1);
    end
  else
    if ( preShiftedOrigin)
      rV = Hqz(oX:end,oY);
    else
      rV = Hqz(1:oX,1);
    end
  end
  
  firstZero = find(rV > 0, 1,'first');
  if isempty(firstZero) || firstZero < floor(0.1.*CTFSIZE(1))
    firstMin = floor(CTFSIZE(1)/2)-6;
  else
    [~,firstMin]=min(abs(rV(7:firstZero-1)-rV(8:firstZero)));
    firstMin = firstMin + 6;
  end
  
  
  if ( preShiftedOrigin && ~calcOneD)
    if doHalfGrid
      freqMin  = radialGrid(firstMin,ceil((CTFSIZE(1)+1)./2));
      maxRes = 0.5./radialGrid(ceil((CTFSIZE(1)+1)/2),1);
    else
      try
        freqMin  = radialGrid(ceil((CTFSIZE(1)+1)./2)+firstMin,ceil((CTFSIZE(2)+1)./2));
      catch
        ceil((CTFSIZE(1)+1)./2)
      end
      maxRes = 0.5./radialGrid(1,ceil((CTFSIZE(2)+1)/2));
    end
    
  else
    freqMin  = radialGrid(firstMin,1);
    freqZero = radialGrid(firstZero,1);
    maxRes = 0.5./radialGrid(ceil((CTFSIZE(1)+1)/2),1);
  end
  
  
  if (thisZero > 0)
    lowCut = 1./(0.1*freqMin+0.9*freqZero);
    if isempty(lowCut)
      lowCut = 2*maxRes;
    end
    if thisZero > 1
      bFactor = thisZero;
    else
      bFactor = 100;
    end
    ctfMask = BH_bandpass3d(size(Hqz),0,800,lowCut,'GPU',maxRes);
    % This term is straight from dTegunov's deconv
    snr = 10.^3.*exp((-2.2.*bFactor).*radialGrid);
    ctfMask = ctfMask .* Hqz ./ (Hqz.^2 + 1./snr);
    snr = [];
  else
    if (flgDampen)
      envelope = (exp(-20.*(radialGrid.*(0.5/max(radialGrid(:)))).^1.25) +0.1)./1.1;
    else
      envelope = 1;
    end
    
    try
      ctfMask = (envelope).*(sign(Hqz).*(radialGrid <= freqMin ).*abs(Hqz).^abs(Phase_Only) + (radialGrid > freqMin).*Hqz);
    catch
      size(rV)
      size(radialGrid)
      firstMin
      freqMin
      firstZero
      thisZero
      maxRes
      error('sfd')
    end
  end
  
  
elseif Phase_Only == 1
  
  ctfMask = sign(Hqz);
  
else
  ctfMask = Hqz;
end



end
