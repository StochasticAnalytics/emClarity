function COM = EMC_centerOfMass(IMAGE, ORIGIN)
%
% COM = EMC_centerOfMass(IMAGE, MASK, ORIGIN)
% Compute the center of mass of real space 2d/3d IMAGE.
%
% Input:
%   IMAGE (numeric):    2d/3d image.
%
%   ORIGIN (int):       Origin convention; 0, 1 or 2;
%                       The center of mass (COM) is relative to this origin.
%                       See EMC_coordVectors for more details.
%
% Output:
%   COM (row vector):   Center of mass of the IMAGE; 2d:[x,y] or 3d:[x,y,z]
%                       NOTE: it has the same precision and method as the IMAGE.
%                             Use cast and EMC_setMethod to cast/push|gather.
%
% Other EMC-files required:
%   EMC_getClass, EMC_coordVectors
%

% Created:  18Jan2020, R2019a
% Version:  v.1.0.  unittest (TF, 4Feb2020).
%           v.1.1.  switch from computing the ndgrid to vector broadcasting.
%                   It is slightly slower with large IMAGEs, but it is more
%                   memory efficient as the ndgrids don't need to be computed (TF, 18Jan2020).
%           v.1.2.  negative values are correctly handled. IMAGE with unique values
%                   (even 0) are allowed (TF, 28Jan2020).
%

%% MAIN
if ~isnumeric(IMAGE)
    error('EMC:IMAGE', 'IMAGE should be numeric, got %s', class(IMAGE))
elseif isvector(IMAGE)
    error('EMC:IMAGE', 'IMAGE should be a 2d or 3d array, got vector')
end

if ~isscalar(ORIGIN) || ~isnumeric(ORIGIN)
    error('EMC:origin', 'ORIGIN should be an integer')
elseif ORIGIN ~= 1 && ORIGIN ~= 2 && ORIGIN ~= 0
    error('EMC:origin', 'ORIGIN should be 0, 1 or 2, got %02f', ORIGIN)
end

imgSize = size(IMAGE);
[precision, isOnGpu] = EMC_getClass(IMAGE);
if isOnGpu
    [vX, vY, vZ] = EMC_coordVectors(imgSize, 'gpu', {'origin', ORIGIN; 'precision', precision});
else
    [vX, vY, vZ] = EMC_coordVectors(imgSize, 'cpu', {'origin', ORIGIN; 'precision', precision});
end

IMAGE = IMAGE - min(IMAGE, [], 'all');  % min to 0
total = sum(IMAGE, 'all');
if total == 0
    if EMC_is3d(imgSize)
        COM = [sum(vX)/imgSize(1), sum(vY)/imgSize(2), sum(vZ)/imgSize(3)];
    else
        COM = [sum(vX)/imgSize(1), sum(vY)/imgSize(2)];
    end
    return
end

if EMC_is3d(imgSize)
    COM = [sum(IMAGE.*vX', 'all'), sum(IMAGE.*vY, 'all'), sum(IMAGE.*reshape(vZ,1,1,[]), 'all')] ./ total;
else
    COM = [sum(IMAGE.*vX', 'all'), sum(IMAGE.*vY, 'all')] ./ total;
end

if any(isnan(COM))
    error('EMC:IMAGE', 'IMAGE has a NaN and/or Inf value, resulting into a COM of %s', mat2str(COM))
end

end  % EMC_centerOfMass
