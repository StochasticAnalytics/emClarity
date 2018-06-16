function [ cRefFilter, cRefBfactorFilter, targetFuncShells ] = ...
                                         BH_multi_cRef( fscINFO, radialGrid, ...
                                                        bFactor, flgAlignRef, ...
                                                        varargin)
%Calculate cRef bandpass filters.
%   Detailed explanation goes here

if (nargin == 5)
  calcTargetFunction = 1;
else
  calcTargetFunction = 0;
  targetFuncShells = '';
end
% First check to see if any cones were calculated in addition to spherical
% shells
coneOverlap = zeros(size(radialGrid),'uint8','gpuArray');
% % % if any(fscINFO{2}(:,2:end))
% % %   nCones = fscINFO{7};
% % %   % Don't include the spherical shells
% % %   coneOffset = 1;
% % % else
% % %   nCones = 0;
% % %   coneOffset = 0;
% % % end

% Since this is applied to the raw volumes just use the spherical
nCones = 0;
coneOffset = 0;

coneList = fscINFO{8};
halfAngle= fscINFO{9};

% No longer using for any sharpening ( to delete)
cRefFilter = zeros(size(radialGrid),'single');
%cRefBfactorFilter = zeros(size(radialGrid),'single');
cRefBfactorFilter = '';
  % oversampled frequency axis
  osX = fscINFO{4};
  
% Include an MTF correction for falcon II, based on Henderson's 2014 findings,
% move this to a param file option including other detectors
mtfX = 0:0.5/(length(osX)-1):0.5;
flgPrintUsage =1;
for iFilter = 1+coneOffset:1+nCones
 iFilter
  % Set any negative FSC values to zero prior to fitting to prevent knots in the
  % cRef curve.
  nonNegFSC = fscINFO{2}(:,iFilter);
  nonNegFSC = nonNegFSC .* (nonNegFSC > 0);
  fitFSC = fit(fscINFO{1}(:,iFilter), nonNegFSC,'cubicSpline');
  
  if iFilter > 1
    
    coneOrientation = BH_defineMatrix(coneList{iFilter-1},'Bah','invVector');
    [ radius,~,height,~,~,~ ] = ...
                                BH_multi_gridCoordinates( size(radialGrid), ...
                                                          'Cylindrical', ...
                                                          'GPU', ...
                                                          {'single',...
                                                          coneOrientation,...
                                                          [0,0,0]', ...
                                                          'invVector',...
                                                          1,1},...
                                                          0, 0, 0 );

    iConeMask = (rad2deg(atan2(radius,abs(height))) < halfAngle{iFilter-1});
    clear radius height
    coneOverlap = coneOverlap + uint8(iConeMask);
  else
    iConeMask = 1;
    
  end
  
  % This should be true everywhere except alignemnt for FSC calculation.
  if (flgAlignRef)
    % gaussian fall off to make sure cRef goes to zero if the curve gets wonky.
    % force form fsc = 0.5 
    forceMask = fscINFO{5}{iFilter};
  else
    % force from fsc = 0.143
    forceMask = fscINFO{6}{iFilter};
  end
  
  if (calcTargetFunction && iFilter == 1)
    [targetFuncShells] = calc_shells(radialGrid, forceMask, osX)
  end
    switch fscINFO{3}{3}
      case 0
        adHocMTF = 1;
        if flgPrintUsage 
          fprintf('\n\nUsing MTF 0\n\n');
          flgPrintUsage = 0;
        end
      case 1
        detector = -100;
        capVal = .05;
        adHocMTF =((exp(detector.*osX.^1.25)+capVal)./capVal).^-1;        

        if flgPrintUsage 
          fprintf('\n\nUsing MTF new\n\n');
          flgPrintUsage = 0;
        end     
      case 2
        detector = -25;
        capVal = .06;
        adHocMTF = ((exp(detector.*osX.^1.25)+capVal)./capVal).^-1;
        
        if flgPrintUsage 
          fprintf('\n\nUsing MTF orig\n\n');
          flgPrintUsage = 0;
        end
      otherwise
        detector = -1.*round(fscINFO{3}{3});
        capVal = fscINFO{3}{3}-round(fscINFO{3}{3});
        adHocMTF = ((exp(detector.*osX'.^1.25)+capVal)./capVal).^-1;
        
        if flgPrintUsage 
          fprintf('\n\nUsing MTF exp %d %3.3f\n\n',detector,capVal);
          flgPrintUsage = 0;
        end       
        
    end


  cRef  =  fit(osX, sqrt(abs( 2.*fitFSC(osX)./(1+fitFSC(osX)))) .* ...
                 adHocMTF .* forceMask, 'cubicSpline'); 


% % %   cRef= fit(osX, ((exp(-25.*mtfX'.^1.25)+0.06)./1.06).^-1 .* ...
% % %                  sqrt(abs( 2.*fitFSC(osX)./(1+fitFSC(osX)))) .*  ...
% % %                    forceMask, 'cubicSpline');
  cRefFilter = cRefFilter + ...
               iConeMask .* reshape(cRef(radialGrid),size(radialGrid));  
%   figure, plot(osX,cRef(osX)); title(sprintf('fit %d',iFilter));

% % %   if (bFactor)
% % %     
% % %     cRef= fit(osX, exp(bFactor.*osX.^2) .* ...
% % %                    ((exp(-25.*mtfX'.^1.25)+0.06)./1.06).^-1 .* ...
% % %                    sqrt(abs( 2.*fitFSC(osX)./(1+fitFSC(osX)))) .*  ...
% % %                    forceMask, 'cubicSpline');
% % %     cRefBfactorFilter = cRefBfactorFilter+ ...
% % %                        iConeMask .* reshape(cRef(radialGrid),size(radialGrid));      
% % %   end
  clear iConeMask
end

if (nCones)
  coneOverlap = single(coneOverlap);
%  figure, imshow3D(coneOverlap);
  divZeroMask = (coneOverlap~=0);
  cRefFilter(divZeroMask)  = cRefFilter(divZeroMask) ./coneOverlap(divZeroMask);
% % %   cRefBfactorFilter(divZeroMask)  = cRefBfactorFilter(divZeroMask) ./coneOverlap(divZeroMask);
  
  [ gaussKernel ] = gpuArray(BH_multi_gaussian3d(7, 1.5 ));

  cRefFilter = ifftshift(gather(convn(fftshift(cRefFilter), gaussKernel, 'same')));
% % %   cRefBfactorFilter = ifftshift(gather(convn(fftshift(cRefBfactorFilter), gaussKernel, 'same')));
end

cRefFilter = gather(cRefFilter);
% % % cRefBfactorFilter = gather(cRefBfactorFilter);
                                  
clear fitFSC osX forceMask cRef radialGrid fscINFO                                            

end

function [targetFuncShells] = calc_shells(radialGrid, forceMask, osX)
                       

                     
padDIM = size(radialGrid)
rad = radialGrid(1,1:ceil(padDIM(2)/2),1);
% Lump all lower resolution info into one shell and don't go out further than
% necessary

highResCutoff = osX(find(forceMask < 1e-3, 1, 'first'));
lowResCutoff = (10+2/highResCutoff)^-1;
fprintf('%f high %f low\n',1/highResCutoff,1/lowResCutoff);
highIDX = gather(find(rad > highResCutoff, 1, 'first'))
lowIDX = gather(find(rad > lowResCutoff,1,'first'))
idxVector = gpuArray(uint32([1:numel(radialGrid)]'));

% check a few bin sizes
if (highIDX - lowIDX)/5 < 30
  binInc = 5;
else 
  binInc = floor((highIDX-lowIDX)/30);
end  
% The first two positions are a "header" for 0/1 to take abs value, and the
% "wiener filter constant" = 1/(N included)^3/2
targetFuncShells = cell(binInc+2,1);



nIDX = lowIDX;
nShell = 4;

iMask = (0 <= radialGrid & radialGrid < rad(lowIDX));
targetFuncShells{1} = uint32(1); % for now default to taking abs value
targetFuncShells{2} = (0);
targetFuncShells{3} = gather(idxVector(iMask(:)));

while nIDX < highIDX


    iMask = ( rad(nIDX) <= radialGrid & radialGrid < rad(nIDX+binInc));
    
    targetFuncShells{nShell} = gather(idxVector(iMask(:)));
    targetFuncShells{2} = targetFuncShells{2} + length(targetFuncShells{nShell});
        nShell = nShell + 1;
        nIDX = nIDX + binInc;
end


end % end of the calc_shells function/

