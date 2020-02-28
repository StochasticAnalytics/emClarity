function  [ reconstruction ] = fourierCtfRec(wantedSize, positionList, TLT, reconGeometry, originPrj, varargin)


   % % % % 
  % For testing override inputs
  % % % %
  flgTesting = false;
  
 
% Insert a check forces a cubic dimension during reconstruction, followed
  % by rescaling
  d1 = wantedSize(1);
  d2 = wantedSize(2);
  d3 = wantedSize(3);
  
  samplingRate = 1;
  
  exposureWeight = 0;
  if nargin == 6
    iThickness = 75;
    samplingRate = varargin{1};
    fractionOfDose = TLT(:,14)/mean(TLT(:,14));
    fractionOfElastics = exp(-1.*iThickness/( cosd(TLT(:,4))*400 ));
    doHalfGrid = 1;
    centerGrid = 1;
    [ exposureWeight ] = BH_exposureFilter( wantedSize(1:2), TLT, 'GPU', samplingRate, centerGrid, doHalfGrid);  
    for iPrj = 1:size(exposureWeight,3)
      exposureWeight(:,:,TLT(iPrj,1)) = exposureWeight(:,:,TLT(iPrj,1)) .* ...
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
    reconGeometry = wantedSize;
    reconShift = [0,0,0];
    originPrj = wantedSize;
    originPrj = ceil((originPrj+1)./2);
    originPrj(3) = 1;
    % Position list is just for the wanted subTomo
    originVol = ceil((reconGeometry(1,1:3)+1)./2);
    lowerLeftVol = originPrj+reconShift-originVol;

    % This matrix would rotate a vector into the particles reference frame,
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
  

     
  oX = floor(d1/2)+1;
  oY = floor(d2/2)+1;
  oZ = floor(d3/2)+1;

  spreadVal = [1];
  %spreadVal =  [0.4904    0.3733    0.1363]
  nOffCenter = length(spreadVal) - 1;

  % Replicate the taper 
  spreadVal = [flip(spreadVal(2:end)),spreadVal];

  [X,Y] = BH_multi_gridCoordinates([d1,d2],'Cartesian',...
                                           'GPU',{'none'},0,1,0,{'halfgrid'});
                                         
  [rad,phi] = BH_multi_gridCoordinates([d1,d2],'Cylindrical',...
                                           'GPU',{'none'},1,1,0,{'halfgrid'});    

  % Assuming the pixel size is constant.
  rad = {rad(:) ./ (samplingRate.*TLT(1,16)),1,phi(:)};
  clear phi

  spatialToPixel = 0.5;
  
  X = X(:) + spatialToPixel;
  Y = Y(:) + spatialToPixel;

  reconstruction = zeros([oX,d2,d3],'single','gpuArray');

  weights = reconstruction;

for iPrj = 1:size(TLT,1)
  

  % Need a defocus offset based on XYZ position in the tomogram
  rTilt = BH_defineMatrix([0,TLT(iPrj,4),0],'SPIDER','forwardVector');
  prjCoords = rTilt * prjVector';

  iDefocusOffset = prjCoords(3).*(TLT(iPrj,16));
  
  iDefocus = [iDefocusOffset + TLT(iPrj,15) - TLT(iPrj,12), ...
              iDefocusOffset + TLT(iPrj,15) + TLT(iPrj,12), ...
              TLT(iPrj,13)];
            
          
            
            
  iCTF = (BH_ctfCalc(rad,TLT(iPrj,17),TLT(iPrj,18),iDefocus,[oX,d2],TLT(iPrj,19),-1)).^2;
  if (length(exposureWeight) > 1)
    iCTF = iCTF(:).*reshape(exposureWeight(:,:,TLT(iPrj,1)),oX.*d2,1);
  else
    iCTF = iCTF(:);
  end
  
  % The sample is rotated by tiltA, which is like rotating the beam by
  % -tiltA. We want to insert in the plane perpendicular to the beam.
  xform = particleRotMat *rTilt; 

  for iOffset = 1:nOffCenter*2+1

    iZ = iOffset - nOffCenter - 1 + spatialToPixel;
    Xnew = X.*xform(1) + Y.*xform(4) + (iZ*xform(7)) ;
    Ynew = X.*xform(2) + Y.*xform(5) + (iZ*xform(8)) ;
    Znew = X.*xform(3) + Y.*xform(6) + (iZ*xform(9)) ;
    
    % For any points rotated into the unsampled half, generate Hermitian
    % mate. The CTF is real, so no need to also negate that, however, for
    % reconstruc
    hermitianSym = Xnew < 0;
    Xnew(hermitianSym) = -1.*Xnew(hermitianSym);
    Ynew(hermitianSym) = -1.*Ynew(hermitianSym);
    Znew(hermitianSym) = -1.*Znew(hermitianSym);
    
    % Need to use the conjugate!!!
    
    Xnew = Xnew + 1 ;
    Ynew = Ynew + oY ;
    Znew = Znew + oZ ;

    for dZ = [0:1]
      zCoord = (floor(Znew)+dZ);
      z_inBounds = (zCoord <= d3 & zCoord > 0);   
      weight_z = abs(1.0 - abs(Znew-zCoord));

      for dY = [0:1]
        yCoord = (floor(Ynew)+dY);
        weight_yz = abs(1.0 - abs(Ynew-yCoord)) .* weight_z;

        y_inBounds = (z_inBounds & yCoord <= d2 & yCoord > 0);

        for dX = [0:1]
  
          xCoord = (floor(Xnew)+dX);
          weight_xyz = weight_yz .* (spreadVal(iOffset) *  abs(1.0 - abs(Xnew-xCoord))) ;
          
          x_inBounds = ( y_inBounds & xCoord <= oX & xCoord > 0 );          

          linearIDX = sub2ind([oX,d2,d3],xCoord(x_inBounds), ...
                                         yCoord(x_inBounds), ...
                                         zCoord(x_inBounds));

          weight_xyz = weight_xyz(x_inBounds);

          xCTF = weight_xyz.*iCTF(x_inBounds);
          
                                                     
          toKeep = true(length(linearIDX),1,'gpuArray');

   
     % Currently ~ 2.7s for a 257^3 - half that time is in this block.
     % Unique is the biggest killer, then assignments out.
     % If all parts of loop could be assigned to a single output grid which
     % is already unique, this could save most of this time.
          while ~isempty(linearIDX)
            
            [~,ai] = unique(linearIDX);
            ci = linearIDX(ai);

            reconstruction(ci) = reconstruction(ci) + xCTF(ai);
            weights(ci) = weights(ci) + weight_xyz(ai);

            toKeep(ai) = false;
            
            linearIDX = linearIDX(toKeep);
            weight_xyz = weight_xyz(toKeep); 
            xCTF = xCTF(toKeep);
            toKeep = toKeep(toKeep);
            
          end



        end
      end
    end 
    
  end % end of offset loop
  


end

% TODO these should be set up top, and possible outside.
gaussDiam = 5;
g = BH_multi_gaussian3d(gaussDiam.*[1,1,1],1.5);
% % TODO this should be worked in to avoid extra allogation.
% reconstruction = BH_padZeros3d(reconstruction,gaussDiam.*[1,1,1], ...
%                                               gaussDiam.*[1,1,1], ...
%                                               'GPU','single');
% weights = BH_padZeros3d(weights,gaussDiam.*[1,1,1], ...
%                                               gaussDiam.*[1,1,1], ...
%                                               'GPU','single');                                            
% % Mirror along the origin                                            
% reconstruction(1:gaussDiam,:,:) =  reconstruction(2*gaussDiam:-1:gaussDiam+1,:,size(reconstruction,3):-1:1); 
% weights(1:gaussDiam,:,:) =  weights(2*gaussDiam:-1:gaussDiam+1,:,size(weights,3):-1:1);                                           

reconstruction = convn(reconstruction,g,'same')./convn((weights+.01),g,'same');
% % Trim the fat
% reconstruction = reconstruction(gaussDiam+1:end-gaussDiam, ...
%                                 gaussDiam+1:end-gaussDiam, ...
%                                 gaussDiam+1:end-gaussDiam);

clear weights

end % end of function

