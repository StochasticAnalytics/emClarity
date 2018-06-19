function [ Gc1,Gc2,Gc3,g1,g2,g3 ] = BH_multi_gridCoordinates( SIZE, SYSTEM, METHOD, ...
                                                  TRANSFORMATION, ...
                                                  flgFreqSpace, ...
                                                  flgShiftOrigin, flgRad, ...
                                                  varargin)
%Return grid vectors in R3 for various coordinate systems.
%   Create grid vectors of dimension SIZE, that are either Cartesian,
%   Cylindrical, or Spherical. Optionally only return a matrix with radial
%   values. Grids are centered with the origin at ceil(N+1/2). 
if strcmpi(METHOD,'GPU')
  SIZE = gpuArray(single(SIZE));
else
  SIZE = single(SIZE);
end

% Option to use pre-created vectors which is surprisingly expensive to
% create.
makeVectors = 1;
if nargin == 8
  makeVectors = 0;
  x1 = varargin{1}{1};
  y1 = varargin{1}{2};
  z1 = varargin{1}{3};
end
  
if numel(SIZE) == 3
  sX = SIZE(1) ; sY = SIZE(2) ; sZ = SIZE(3);
  
  flg3D = 1;
elseif numel(SIZE) == 2
  sX = SIZE(1) ; sY = SIZE(2) ; 
  if strcmpi(METHOD,'GPU'); sZ = gpuArray(single(1)); else sZ = single(1);end
  flg3D = 0;
else
  error('Only 2D or 3D grid vectors are supported.');
end

% flgShiftOrigin handles whether to place the origin at the center of the image
% (TRUE) or not (FALSE; use for FFTs) and optionally to handle offsets due to
% different conventions in other software. E.g. IMOD defines the center as NX/2
% which creates a different value for even/odd images.
if length(flgShiftOrigin) == 4
  % the origin in IMODs case for an even image is -0.5 relative to mine.
  % Switching to force odd size - 20171201
  conventionShift = [0,0,0];
%   conventionShift = flgShiftOrigin(2:4) .* (1-mod([sX,sY,sZ],2));
  flgShiftOrigin = flgShiftOrigin(1);
else
  conventionShift = [0,0,0];
end
flgTrans = 1;
flgGridVectors = 0;
flgSymmetry = 0;
flgSequential = 0;

symInc = 0;
symIDX = 0;

flgMask = 0;
if iscell(TRANSFORMATION) 
  
  switch TRANSFORMATION{1}
    
    case 'none'
      flgTrans = 0;
      R = [1,0,0;0,1,0;0,0,1];
      dXYZ = [0,0,0]';
      DIR = 'forwardVector';
      MAG = {1};
      % The majority of function calls that are not in a resample/rescale
      % program call this case, and don't expect a cell output.
      
      
    case 'gridVectors'
      flgGridVectors = 1;
      R = [1,0,0;0,1,0;0,0,1];
      dXYZ = TRANSFORMATION{3};
      DIR = 'forwardVector';
      MAG = {TRANSFORMATION{6}};
      
    case 'single'      
      if numel(TRANSFORMATION{2}) == 9
        R = reshape(TRANSFORMATION{2},3,3);
      else
        R = reshape(TRANSFORMATION{2},2,2);
      end
        
      dXYZ = TRANSFORMATION{3};
      if length(dXYZ) == 2
        dXYZ = [dXYZ;0];
      end
      DIR = TRANSFORMATION{4};
      if (TRANSFORMATION{5} > 1)
        flgSymmetry = TRANSFORMATION{5};
        symInc = 360 / flgSymmetry;
        symIDX = 0:flgSymmetry-1;
        Gc1 = cell(flgSymmetry,1);
        Gc2 = cell(flgSymmetry,1); 
        Gc3 = cell(flgSymmetry,1);
      else
        symIDX = 1;
        symInc = 0;
      end
      
      if strcmpi(DIR, 'inv') || strcmpi(DIR,'forwardVector')
        MAG = {TRANSFORMATION{6}};
      else
        MAG = {1./TRANSFORMATION{6}}; % faster to just do A(I) but left as {{}} for clarity
      end
      
      if length(TRANSFORMATION) == 7
        if islogical(TRANSFORMATION{7})
          binaryVol = TRANSFORMATION{7};
          flgMask = 1;
        end
      end
      


      
    case 'sequential'
      flgSequential = 1;
      
      if ( size(TRANSFORMATION,1) > 2 )
        error('Initial only implement up to 2 sequential transformations.\n')
      elseif (TRANSFORMATION{1,5} ~= 1)
        error('no Symmetry operations on sequential transformations.\n')
      elseif ~strcmpi(TRANSFORMATION{1,4}, 'forward')
        error('only forward transformations supported for sequential.\n')
      else
        nTrans = 2;
      end
      
      for iTrans = 1:nTrans
        R_seq{iTrans} = reshape(TRANSFORMATION{iTrans,2},3,3);
      
        dXYZ_seq{iTrans} = TRANSFORMATION{iTrans,3};
    
       
        % Convention is only forward so np need to consider flipping
        MAG_seq{iTrans} = TRANSFORMATION{iTrans,6};
             
      end
      % For now assuming no symmetry operation on sequential transformations
      flgSymmetry = TRANSFORMATION{1,5};
      DIR = TRANSFORMATION{1,4};

      symInc = 360 / flgSymmetry;
      symIDX = 0:flgSymmetry-1;
      Gc1 = {}; Gc2 = {}; Gc3 ={};
      
     
      
    otherwise
      error(['TRANSFORMATION must be a cell,',...
            '(none,gridVectors,single,sequential),',...
            'Rotmat, dXYZ, forward|inv, symmetry\n']);
  end
end


if ( makeVectors )
  if strcmpi(METHOD, 'GPU')
   % sX = gpuArray(sX) ; sY = gpuArray(sY) ; sZ = gpuArray(sZ);
     if flgShiftOrigin == 1
  %      x1 = [-1*floor((sX)/2):0,1:floor((sX-1)/2)];
  %      y1 = [-1*floor((sY)/2):0,1:floor((sY-1)/2)];
       x1 = [-1*floor((sX)/2):floor((sX-1)/2)];
       y1 = [-1*floor((sY)/2):floor((sY-1)/2)];
       if (flg3D); z1 = [-1*floor((sZ)/2):floor((sZ-1)/2)]; end

  %      if (flg3D); z1 = [-1*floor((sZ)/2):0,1:floor((sZ-1)/2)]; end
     elseif flgShiftOrigin == -1
       x1 = [1:sX];
       y1 = [1:sY];
       if (flg3D); z1 = [1:sZ]; end
     elseif flgShiftOrigin == -2
       x1 = fftshift([1:sX]);
       y1 = fftshift([1:sY]);
       if (flg3D); z1 = fftshift([1:sZ]); end     
     else
       x1 = [0:floor(sX/2),-1*floor((sX-1)/2):-1];      
       y1 = [0:floor(sY/2),-1*floor((sY-1)/2):-1];  
       if (flg3D); z1 = [0:floor(sZ/2),-1*floor((sZ-1)/2):-1]; end     
     end
  elseif strcmpi(METHOD, 'cpu')
     if flgShiftOrigin == 1
       x1 = [-1*floor((sX)/2):0,1:floor((sX-1)/2)];
       y1 = [-1*floor((sY)/2):0,1:floor((sY-1)/2)];
       if (flg3D); z1 = [-1*floor((sZ)/2):0,1:floor((sZ-1)/2)]; end
     elseif flgShiftOrigin == -1
       x1 = [1:sX];
       y1 = [1:sY];
       if (flg3D); z1 = [1:sZ]; end  
     elseif flgShiftOrigin == -2
       x1 = fftshift([1:sX]);
       y1 = fftshift([1:sY]);
       if (flg3D); z1 = fftshift([1:sZ]); end       
     else
       x1 = [0:floor(sX/2),-1*floor((sX-1)/2):-1];      
       y1 = [0:floor(sY/2),-1*floor((sY-1)/2):-1];  
       if (flg3D); z1 = [0:floor(sZ/2),-1*floor((sZ-1)/2):-1]; end   
     end
  else
    error('METHOD should be cpu or GPU')
  end
end
  
% Make any needed shifts for convention
x1 = x1 - conventionShift(1);
y1 = y1 - conventionShift(2);
if (flg3D); z1 = z1 - conventionShift(3); end

if strcmpi(DIR, 'inv') || strcmpi(DIR, 'forwardVector')
  x1 = x1 - dXYZ(1);
  y1 = y1 - dXYZ(2);
  if (flg3D); z1 = z1 - dXYZ(3); end
  
  shiftDir = -1;
else
  shiftDir = 1;
end

if (flgFreqSpace)
  x1 = x1 ./ sX;
  y1 = y1 ./ sY;
  if (flg3D); z1 = z1 ./ sZ; end
end


if ~(flg3D)
  z1 = 0;
end
% save grid vectors prior to any transformations
g1 = x1; g2 = y1; g3 = z1;

if (flgGridVectors)
  G1 = ''; G2 = ''; G3 = '';
  return
end

% Rescale the vectors prior to making gridVectors 
if (flgSequential)
  x1 = x1.*MAG_seq{1};
  y1 = y1.*MAG_seq{1};
  z1 = z1.*MAG_seq{1};
else
  x1 = x1.*MAG{1};
  y1 = y1.*MAG{1};
  z1 = z1.*MAG{1};  
end

% No matter the case, the cartesian grids are needed
[X,Y,Z] = ndgrid(x1,y1,z1);


% Optionally evaluate only a smaller masked region
if (flgMask) 
  X = X(binaryVol);
  Y = Y(binaryVol);
  Z = Z(binaryVol);
end

% for non-sequential transformations this only loops once.
for iTrans = 1:1+(flgSequential)
  
  if (flgSequential)
    rAsym = R_seq{iTrans};  
    if (iTrans == 1)
      dXyzAsym = 0;
      % Instead of shifting then shifting back, just note the original shift
    else
       % adding the R2' because the first term doesn't need to be multiplied by
       % R2 (in the first action under symmetry loop) but making a change here
       % which involves extra multiplications is okay, since this function is
       % used much less than 'single' style resampling.
       dXyzAsym = ( R_seq{iTrans}'*R_seq{iTrans-1} * ...
                    dXYZ_seq{iTrans-1}.*MAG_seq{iTrans-1} + ...
                    R_seq{iTrans-1} * dXYZ_seq{iTrans}.*MAG_seq{iTrans} );
                  
      X = Gc1{1}.*MAG_seq{iTrans};
      Y = Gc2{1}.*MAG_seq{iTrans};
      Z = Gc3{1}.*MAG_seq{iTrans};
    end
  else
    rAsym = R;
    % Note that if symmetric, dXYZ changes each loop after this point
    dXyzAsym = dXYZ.*MAG{1};
  end

  for iSym = symIDX

    % Only in plane symmetries considered anywhere so inv|forward shouldn't
    % matter.

    R = rAsym * BH_defineMatrix([iSym.*symInc,0,0],'Bah','inv');

    % Any forward transformations of the grids
    if (flgTrans)
      dXYZ = shiftDir .* R*dXyzAsym;

      Xnew = X.*R(1) + Y.*R(4) + Z.*R(7) - dXYZ(1);
      Ynew = X.*R(2) + Y.*R(5) + Z.*R(8) - dXYZ(2);
      if (flg3D)
        Znew = X.*R(3) + Y.*R(6) + Z.*R(9) - dXYZ(3);
      else
        Znew = 0;
      end
    else
      Xnew = X; Ynew = Y ; Znew = Z;
    end

    % Only return the radial grid if requested
    if (flgRad)
      G1 = sqrt(Xnew.^2 + Ynew.^2 + Znew.^2);
      G2 = '';
      G3 = '';

    else
      switch SYSTEM
        case 'Cartesian'
          G1 = Xnew; G2 = Ynew; G3 = Znew;
        case 'Spherical'
          G1 = sqrt(Xnew.^2 + Ynew.^2 + Znew.^2);
          G2 = atan2(Ynew,Xnew);
          % set from [-pi,pi] --> [0,2pi]
          G2(G2 < 0) = G2(G2 < 0) + 2.*pi;
          % [0,pi]
          G3 = acos(Z./G1);

        case 'Cylindrical'

          G1 = sqrt(Xnew.^2 + Ynew.^2);
          G2 = atan2(Ynew,Xnew);
          % set from [-pi,pi] --> [0,2pi]
          G2(G2 < 0) = G2(G2 < 0) + 2.*pi;
          G3 = Znew;

        otherwise
          error('SYSTEM must be Cartesian, Spherical, Cylindrical')
       end
    end
    % Only use as cell if symmetry is requested
    if (flgSymmetry)
      % 
      
      Gc1{iSym+1} = G1;
      Gc2{iSym+1} = G2;
      Gc3{iSym+1} = G3;
    else
      Gc1 = G1 ; clear G1
      Gc2 = G2 ; clear G2
      Gc3 = G3 ; clear G3
    end
  end % loop over symmetric transformations
end % loop over sequential transformations

  

clear X Y Z  Xnew Ynew Znew x1 y1 z1

end


