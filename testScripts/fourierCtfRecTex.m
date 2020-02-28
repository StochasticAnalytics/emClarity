function  [ reconstruction ] = fourierCtfRecTex(wantedSize, positionList, TLT, reconGeometry, originPrj, varargin)


   % % % % 
  % For testing override inputs
  % % % %

  if length(varargin) > 1
    flgTesting = varargin{1,2};
  else
    flgTesting = false;
  end
  
 
% Insert a check forces a cubic dimension during reconstruction, followed
  % by rescaling
  d1 = wantedSize(1);
  d2 = wantedSize(2);
  d3 = wantedSize(3);
  
  samplingRate = 1;

  exposureWeight = 0;
  if nargin == 6 && ~flgTesting
    iThickness = 75;
    samplingRate = varargin{1};
    fractionOfDose = TLT(:,14)./mean(TLT(:,14));
    fractionOfElastics = exp(-1.*iThickness./( cosd(TLT(:,4))*400 ));
    fractionOfElastics = fractionOfElastics ./ max(fractionOfElastics);
    doHalfGrid = 1;
    centerGrid = 1;
    [ exposureWeight ] = BH_exposureFilter( wantedSize(1:2), TLT, 'GPU', samplingRate, centerGrid, doHalfGrid);  
    for iPrj = 1:size(exposureWeight,3)
      exposureWeight(:,:,TLT(iPrj,1)) = exposureWeight(:,:,TLT(iPrj,1)).^2 .* ...
                                        (fractionOfDose(iPrj).*(fractionOfElastics(iPrj)));
    end
  end


  
  if (d1 ~= d2 || d1 ~= d3)
    error('fix the cubic check ya twit');
  end


  if ( flgTesting )
    % These will just go in the object properties
    rawtlt = -60:3:60;
    nTilts = length(rawtlt);
    TLT = zeros(nTilts,23,'single','gpuArray');
    TLT(:,4) = rawtlt;
    TLT(:,1) = 1:length(rawtlt);
    % 12,13 = astigmatism = 0
    % def,pix,cs,wl,ampcont
    TLT(:,[15,16,17,18,19]) = repmat([-1.2e-6,2.0e-10,2.7e-3,1.969e-12,1e-1],nTilts,1);
    TLT(:,11) = [[41:-1:20,1:19].*3]';
    iThickness = 75;
    samplingRate = 1;
    fractionOfDose =1;
    fractionOfElastics = exp(-1.*iThickness/( cosd(TLT(:,4))*400 ));
    doHalfGrid = 1;
    centerGrid = 1;
    reconGeometry = wantedSize;
    reconShift = [0,0,0];
    originPrj = wantedSize;
    originPrj = ceil((originPrj+1)./2);
    originPrj(3) = 1;
    % Position list is just for the wanted subTomo
    originVol = ceil((reconGeometry(1,1:3)+1)./2);
    lowerLeftVol = originPrj+reconShift-originVol;
    [ exposureWeight ] = BH_exposureFilter( wantedSize(1:2), TLT, 'GPU', samplingRate, centerGrid, doHalfGrid); 
    for iPrj = 1:size(exposureWeight,3)
      exposureWeight(:,:,TLT(iPrj,1)) = exposureWeight(:,:,TLT(iPrj,1)).^2 .* ...
                                        (fractionOfDose.*(fractionOfElastics(iPrj)));
    end    % This matrix would rotate a vector into the particles reference frame,
    % we want to rotate back, so take the transpose to invert.
    particleRotMat = BH_defineMatrix([0,0,0],'Bah','forward'); 
    prjVector = [0,0,0] - originVol + reconShift + [0.0,0.0,1.0];
  else
    % These will just go in the object properties
    reconShift = reconGeometry(2,:);
    originPrj = ceil((originPrj+1)./2);
    originPrj(3) = 1;
    % Position list is just for the wanted subTomo
    originVol = ceil((reconGeometry(1,1:3)+1)./2);
    lowerLeftVol = originPrj+reconShift-originVol;


    % This matrix would rotate a vector into the particles reference frame,
    % we want to rotate back, so take the transpose to invert.
    particleRotMat = reshape(positionList(17:25),3,3)';
    
    % Taken from synthetic mapback, not sure why I'm adding the 1 here.
    prjVector = positionList(11:13) - originVol + reconShift + [0.0,0.0,1.0];
    % We want the origin of the subTomogram in the tilted image to calculate
    % a defocus offset.
  end

  nZpad = 9;
  if mod(nZpad+1,2)
    error('nZpad is assmued odd');
  end
  
  [ taper ] = BH_multi_calcTaper(floor(nZpad./2));
  taper = [ flip(taper), 1, taper ];
  taper = taper./sum(taper(:));


     
  oX = floor(d1/2)+1;
  oY = floor(d2/2)+1;
  oZ = floor(d3/2)+1;


  [X,Y,Z] = BH_multi_gridCoordinates([d1,d2,nZpad],'Cartesian',...
                                           'GPU',{'none'},0,1,0,{'halfgrid'});
                                         
  [rad,phi] = BH_multi_gridCoordinates([d1,d2],'Cylindrical',...
                                           'GPU',{'none'},1,1,0,{'halfgrid'});    

  % Assuming the pixel size is constant.
  rad = {(rad./ (samplingRate.*TLT(1,16).*10^10)),[1,1],1.*phi};
  clear phi

  spatialToPixel = 0.5;
  
  xIn = X + 1;
  yIn = Y + oY;
  zIn = Z + floor(nZpad/2);
  
  wIn = ones(size(X));
  for iTaper = 1:length(taper)
    wIn(:,:,iTaper) = taper(iTaper);
  end

  
  X = X + spatialToPixel;
  Y = Y + spatialToPixel;
  Z = Z + spatialToPixel;

  reconstruction = zeros([oX,d2,d3],'single','gpuArray');

  weights = reconstruction;


for iPrj = 1:size(TLT,1)
 
 
  % Need a defocus offset based on XYZ position in the tomogram
  rTilt =  BH_defineMatrix(TLT(iPrj,4),'TILT','forwardVector') ;
  prjCoords = rTilt * prjVector';
  iDefocusOffset = (prjCoords(3).*(TLT(iPrj,16))+TLT(iPrj,15));
  
  iDefocus = [iDefocusOffset - TLT(iPrj,12), ...
              iDefocusOffset + TLT(iPrj,12), ...
              TLT(iPrj,13)];
            

  % TODO add an option to check and oversample the CTF, then crop in real space, to fix aliasing probs.

  iCTF = (BH_ctfCalc(rad,TLT(iPrj,17),TLT(iPrj,18),iDefocus,[oX,d2],TLT(iPrj,19),-1)).^2;
  
  % In some cases we may be before the first zero. TODO fixme
  if (mean(iCTF(:)) < 0.25 )
    iCTF = iCTF.*0 + 1;
  end
%    iCTF = mexCTF(true,true,int16(d1),int16(d2),single(TLT(iPrj,16)*10^10),single(TLT(iPrj,18)*10^10),single(TLT(iPrj,17)*10^3),single(iDefocus(1)*-1*10^10),single(iDefocus(2)*-1*10^10),single(30),single(0.1));
  if (length(exposureWeight) > 1)
    iWeight = wIn.*repmat(exposureWeight(:,:,TLT(iPrj,1)),1,1,nZpad);
    iCTF    = iWeight .* repmat(iCTF,1,1,nZpad);
  else
    iCTF = repmat(iCTF,1,1,nZpad).*wIn;
  end
  


  
  % The sample is rotated by tiltA, which is like rotating the beam by
  % -tiltA. We want to insert in the plane perpendicular to the beam.
  xform = particleRotMat *rTilt; 

    Xnew = X.*xform(1) + (Z.*xform(7)) ;
    Ynew = Y;
    Znew = X.*xform(3) + (Z.*xform(9)) ;
  
%   iZ = iOffset - nOffCenter - 1 + spatialToPixel;
%   Xnew = X.*xform(1) + Y.*xform(4) + (iZ*xform(7)) ;
%   Ynew = X.*xform(2) + Y.*xform(5) + (iZ*xform(8)) ;
%   Znew = X.*xform(3) + Y.*xform(6) + (iZ*xform(9)) ;

  % For any points rotated into the unsampled half, generate Hermitian
  % mate. The CTF is real, so no need to also negate that, however, for
  % reconstruc
  hermitianSym = find(Xnew < 0);
  Xnew(hermitianSym) = -1.*Xnew(hermitianSym);
  Ynew(hermitianSym) = -1.*Ynew(hermitianSym);
  Znew(hermitianSym) = -1.*Znew(hermitianSym);

  clear hermitianSym

  % Need to use the conjugate!!!

  Xnew = floor(Xnew + 1) ;
  Ynew = floor(Ynew + oY);
  Znew = floor(Znew + oZ) ;
  
  
  xOutOfBounds = Xnew(:,1,:) <= 0 & Xnew(:,1,:) > size(reconstruction,1);
  yOutOfBounds = Ynew(:,1,:) <= 0 & Ynew(:,1,:) > size(reconstruction,2);
  zOutOfBounds = Znew(:,1,:) <= 0 & Znew(:,1,:) > size(reconstruction,3);


  % Not sure what is the most efficient way to do this. All at once, better mem harder sort, or incrementally.

  % For rotation about a single axis, we only need to worry about one value
  % of that axis to deterimine unique positions  
  [~, ia] = unique([reshape(Xnew(:,1,:),oX*nZpad,1),reshape(Znew(:,1,:),oX*nZpad,1)],'rows');
  k = false([oX*nZpad,1],'gpuArray');
  k(ia) = true;
  k(xOutOfBounds | yOutOfBounds | zOutOfBounds) = false;
  
  k = repmat(reshape(k,oX,1,nZpad),1,d2,1);
  k = find(k);
%   A = [Xnew(k),Ynew(k),Znew(k)];
 
%   A = unique([Xnew,Ynew,Znew],'rows');

  Xnew = Xnew(k);
  Ynew = Ynew(k);
  Znew = Znew(k);
  


%   keepVal = Xnew(k) > 0 & Ynew(k) > 0 & Znew(k) > 0 &...
%             Xnew(k) <= size(reconstruction,1) & ...
%             Ynew(k) <= size(reconstruction,2) & ...
%             Znew(k) <= size(reconstruction,3);
%   A = A(keepVal,:);
        linearIDX = sub2ind([oX,d2,d3],Xnew,Ynew,Znew);
%   A = A(keepVal,:);
%         linearIDX = sub2ind([oX,d2,d3],A(:,1),A(:,2),A(:,3));



  % A now has a list of uniqe 3D coords in the reconstruction 
  % Those coords now need to be xformed back into the projection (transpose, should it also be negative?)
  % Optimized for single-axis
  Xnew = Xnew - 1;
  Znew = Znew - oZ;
  xInterp = Xnew.*xform(1) + Znew .* xform(3) + 1;
  zInterp = Xnew.*xform(7) + Znew .* xform(9) + floor(nZpad/2);
%   xInterp = (A(:,1)-1).*xform(1) + (A(:,2)-oY) .*xform(2) + (A(:,3)-oZ) .* xform(3) + 1;
%   yInterp = (A(:,1)-1).*xform(4) + (A(:,2)-oY) .*xform(5) + (A(:,3)-oZ) .* xform(6) + oY;
%   zInterp = (A(:,1)-1).*xform(7) + (A(:,2)-oY) .*xform(8) + (A(:,3)-oZ) .* xform(9) + floor(nZpad/2);




  reconstruction(linearIDX) =  reconstruction(linearIDX)  + ...
                                                    interpn(xIn,yIn,zIn,iCTF,xInterp,Ynew,zInterp,'linear',0);
                                                  
  weights(linearIDX) =  weights(linearIDX)  + ...
                                                    interpn(xIn,yIn,zIn,iWeight,xInterp,Ynew,zInterp,'linear',0);


%   weights(linearIDX) =  weights(linearIDX)  + 1;
      


end

clear iCTF xIn yIn 
  clear xInterp yInterp iWeight

doConv = false;

if doConv
  
% TODO these should be set up top, and possible outside.
gaussDiam = 5;
g = gpuArray(BH_multi_gaussian3d(gaussDiam.*[1,1,1],0.75));

reconstruction = BH_multi_makeHermitian(reconstruction,gaussDiam,1);
weights = BH_multi_makeHermitian(weights + 0.01,gaussDiam,1);
reconstruction = convn(reconstruction,g,'same')./convn(weights,g,'same');
recosntruction = reconstruction(1+gaussDiam:end,:,:);

else                                        
reconstruction = reconstruction ./ (weights + 0.01);
end
% % Trim the fat
% reconstruction = reconstruction(gaussDiam+1:end-gaussDiam, ...
%                                 gaussDiam+1:end-gaussDiam, ...
%                                 gaussDiam+1:end-gaussDiam);

clear weights g bp

end % end of function

