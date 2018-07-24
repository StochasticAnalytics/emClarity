function [ ctfMask, ctfUnMod] = BH_ctfCalc(PixelSize, Cs, Lamda, Defocus, ...
                                CTFSIZE, AMPCONT, Phase_Only, varargin)
% Calculate a ctf and with the given modification to information at resolution
% lower than the first peak. Also return the unmodified ctf.
toAngstrom = 0;
precision = 'single';
flgComplex = 0;
if nargin > 7
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
else
  calcOneD = 0;
end


if AMPCONT < 0
  AMPCONT = abs(AMPCONT);
  flgDampen = 1;
else
  flgDampen = 0;
end

if isa(PixelSize,'cell')
  radialGrid = PixelSize{1};
  preShiftedOrigin = PixelSize{2};
  phi = PixelSize{3};
  calcRad = 0;
  % While double precision is still required at high sampling rates and
  % higher defocus, switching from meters --> angstroms allows single
  % precision with less error, especially useful for non-corrective (i.e.
  % band limited comparisons)
  if numel(preShiftedOrigin) == 2
    toAngstrom = 1;
    preShiftedOrigin = preShiftedOrigin(1);
  end
else
  PIXEL_SIZE = PixelSize;
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

if ( toAngstrom )
  CS = CS * 10^10;
  WL = WL * 10^10;
  df1 = df1 * 10^10;
  df2 = df2 * 10^10;
%   PIXEL_SIZE = PIXEL_SIZE * 10^10;
end

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
%     [radialGrid,phi,~,~,~,~] = ...
%      BH_multi_gridCoordinates(CTFSIZE(1:2),'Cylindrical','cpu',{'none'},1,0,0);
    % Consider adding a flag for wimpier GPUs to do the double precision on
    % cpu if it seems warranted. For now, this is a small calc so don't
    % worry about it.
    [radialGrid,phi,~,~,~,~] = ...
     BH_multi_gridCoordinates(CTFSIZE(1:2),'Cylindrical','GPU',{'none'},1,0,0);
    radialGrid = double(radialGrid);
    phi = double(phi);
  end
  
    radialGrid = radialGrid ./ PIXEL_SIZE;

end





% df1 should be defocus of greater mag and phi0 -90/90

% phasePerturbation = pi.*(0.5.*CS.*WL^3.*(radialGrid).^4 + DF.*WL.*(radialGrid).^2);
dfTerm = 0.5.*( (df1+df2) + (df1-df2)*cos(2.*(phi-phi0)) );
phasePerturbation = pi.*(0.5.*CS.*WL^3.*(radialGrid).^4 + ...
                         WL.*(radialGrid).^2 .* dfTerm);
% dPdQ = 2*pi*CS*WL^3.*radialGrid.^3 + 2*WL.*radialGrid.*dfTerm;                                              
if ( flgComplex )
  ctfMask = exp(-1i.*(phasePerturbation+(pi/2+asin(AMPCONT))));
  ctfUnMod = exp(1i.*(phasePerturbation+(pi/2+asin(AMPCONT))));
  return
else
  Hqz = (sqrt(1-AMPCONT^2).*sin(phasePerturbation) - AMPCONT.*cos(phasePerturbation));
end
% dHqz = dPdQ.*(sqrt(1-AMPCONT^2).*cos(phasePerturbation)+AMPCONT.*sin(phasePerturbation));
% At small pixel size, it is necessary to use double precision to get the
% corners. Could also just bandpass to Nyquist or simply check for NaNs
nanCheck = isnan(Hqz);
if gather(sum(nanCheck(:)))
    Hqz(nanCheck) = 0;
    % Should probably issue some kind of warning.
    
%   randFill = rand(size(Hqz),'single','gpuArray');
%   Hqz(nanCheck) = randFill(nanCheck);
end


if Phase_Only < 0

      if ( preShiftedOrigin )
        rV = Hqz(ceil((CTFSIZE(1)+1)./2),ceil((CTFSIZE(2)+1)./2):end);
      else
        rV = Hqz(1,1:end);
      end


      
      firstZero = find(rV > 0, 1,'first');
      
      if isempty(firstZero)
        %fprintf('No zero crossings in CTF!\n');
        ctfMask = Hqz;
        ctfUnMod = Hqz;
      elseif (firstZero < 9)
        %fprintf('First zero is less than 9 pix\n');
        ctfMask = Hqz;
        ctfUnMod = Hqz; 
 
      else
        
        [~,firstMin]=min(abs(rV(7:firstZero-1)-rV(8:firstZero)));
        firstMin = firstMin + 6;
        
     
        if ( preShiftedOrigin )
          freqMin  = radialGrid(ceil((CTFSIZE(1)+1)./2),firstMin);
        else
          freqMin  = radialGrid(1, firstMin);
        end

        if (flgDampen)
%           envelope = (exp(-100.*(radialGrid.*10^-10).^3) ); 
          envelope = (exp(-20.*(radialGrid.*(0.5/max(radialGrid(:)))).^1.25) +0.1)./1.1;
%           figure, plot(radialGrid(1,1:ceil((CTFSIZE(1)+1)./2)).*10^-10,...
%                        envelope(1,1:ceil((CTFSIZE(1)+1)./2)))
        else
          envelope = 1;
        end
        
% % %         if (calcOneD)
% % %           DTF = sign(sin(pi.*WL.*radialGrid.^2.*maxZ)./(4.*pi^2.*radialGrid.^2));
% % %           DTF(1) = 0;
% % %         else
% % %           DTF = 1;
% % %         end
        DTF = 1;
        try
          ctfMask = (DTF.*envelope).*(sign(Hqz).*(radialGrid <= freqMin ).*abs(Hqz).^abs(Phase_Only) + (radialGrid > freqMin).*Hqz);
          
          ctfUnMod= Hqz;
          
        catch
          SAVE_IMG(MRCImage(gather(single(Hqz))),'ctfCalcErrLine101.mrc');
          fID = fopen('ctfCalcErrLine101.txt','w');
          fprintf(fID,'df1 %f\ndf2 %f\nph0 %f\nfirstzero %f\nfreqmin %f\nfirstmin %f\n',...
           df1,df2,phi0,firstZero,freqMin,firstMin);
          fclose(fID)
        end
      end
      
elseif Phase_Only == 1
    % Binary phase flip filter
    %ctfMask = (Hqz <= 0) - (Hqz > 0);
    %ctfMask = ctfMask .* (radialGrid <= 0.25) + (radialGrid>0.25);
    

                                                   
                                                
    ctfMask = sign(Hqz);
  
else
  ctfMask = Hqz;
end

%     mask = ...
%        ((radialGrid < 0.45) + ...
%          (radialGrid >= 0.45) .* exp(-1.*((radialGrid - 0.45).^2)./((2*2/f1)^2)));
%     mask = ...     
%        ((radialGrid > HIGHFREQ*1.15) + (radialGrid < 0.35) -1 +...
%          (radialGrid >= 0.35) .* exp(-1.*((radialGrid - 0.35).^2)./((2*2/f1)^2)));
%     mask = mask .* (mask > .00001); 
% 
%     ctfMask = ctfMask .* mask;

end
