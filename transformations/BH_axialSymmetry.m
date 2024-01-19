function [ symIMG ] = BH_axialSymmetry( IMG, SYMMETRY, ANG_SHIFT,  ...
  METHOD, XYZ_SHIFT )
%Apply axial (z) symmetry to an image.
%   Goal is to generate either a symmetric reference that is more planar, or a
%   rotationally smeared out reference to use to speed up a grid search by first
%   searching over the polar and azimuthal angles with the rotationally averaged
%   ref, then searching the in-plane angles at the best position only.

overrideRadial = 0;
if SYMMETRY == -720
  overrideRadial = true;
  SYMMETRY = 720;
end
if SYMMETRY > 0
  if SYMMETRY > 60
    if ~(overrideRadial)
      flgRadial = true;
    else
      flgRadial = true;
    end
    % Generating a rotationally smeared out reference, not just applying axial symmetry
    % Create a radial mask similar to ramp weight for WBP
    [ radialCylinder ,~,~,~,~,~ ] = BH_multi_gridCoordinates( size(IMG), ...
      'Cylindrical',...
      METHOD, ...
      {'none'}, ...
      0, 0, 0 );
  else
    flgRadial = false;
  end
  
  symInc = 360/SYMMETRY;
  % An additional axial rotation to orient the reference, use positive sense
  % for indexing
  if ANG_SHIFT < 0
    ANG_SHIFT = 360 + ANG_SHIFT;
  end
  
  if SYMMETRY > 6
    % Old slower way
    symIMG = zeros(size(IMG), 'single','gpuArray');
    for iSym = 0:SYMMETRY-1
      symIMG = symIMG + BH_resample3d(IMG,[ANG_SHIFT+(iSym.*symInc),0,0],XYZ_SHIFT, 'Bah', 'GPU', 'forward');
    end
  else
    % More memory but faster
    symIMG = BH_resample3d(IMG,[ANG_SHIFT,0,0],XYZ_SHIFT, {'Bah',SYMMETRY,'linear'}, 'GPU', 'forward');
  end
  
else
  
  if SYMMETRY == -1 || SYMMETRY == -2 || SYMMETRY == -3
    symIMG = flip(img,abs(SYMMETRY));
  else
    error('mirror symmetry must be -1,-2,-3')
  end
  
end

if (flgRadial)
  symIMG = real(ifftn(fftn(symIMG).*radialCylinder));
end










