function COM = EMC_centerOfMass(IMAGE, ORIGIN)
%
% COM = EMC_centerOfMass(IMAGE, ORIGIN)
% Compute the center of mass of real space 2d/3d IMAGE.
%
% Input:
%   IMAGE (numeric):	2d/3d IMAGE.
%
%   ORIGIN (int):     	Origin convention; 0, 1 or 2;
%                      	See EMC_coordVectors for more details.
%
% Output:
%   COM (row vector):   Center of mass of the IMAGE; 2d:[x,y] or 3d:[x,y,z]
%                       NOTE: it has the same precision and method as the IMAGE.
%                             Use EMC_setPrecision and EMC_setMethod to cast/push|gather.
%
% Note:
%   - The algorithm is quite fast (if IMAGE is a gpuArray), but not memory efficient
%     (use of coordinates grids). One could compute the grids one by one with repmat
%     to be more memory efficient (but slower); for gX: grid = repmat(vX', [1,IMAGEsize(2:3)]).
%
% Other EMC-files required:
%   EMC_getClass, EMC_coordVectors
%

% Created:  18Jan2020, R2019a
% Version:  v.1.0.  unittest (TF, 4Feb2020).
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

IMAGEsize = size(IMAGE);
[precision, isOnGpu] = EMC_getClass(IMAGE);
if isOnGpu
    [vX, vY, vZ] = EMC_coordVectors(IMAGEsize, 'gpu', {'origin', ORIGIN; 'precision', precision});
else
    [vX, vY, vZ] = EMC_coordVectors(IMAGEsize, 'cpu', {'origin', ORIGIN; 'precision', precision});
end

% Normalize to have only positive values.
% NOTE: subtract by min requires more code (if min < 0; if same value) 
%       and speed difference is not significant compared to the rest).
IMAGE = IMAGE ./ max(IMAGE, [], 'all');
IMAGEsum = sum(IMAGE, 'all');
if IMAGEsum == 0
    error('EMC:IMAGE', 'IMAGE is empty (sum of pixels is zero)')
end

if EMC_is3d(IMAGEsize)
    [gX, gY, gZ] = ndgrid(vX, vY, vZ);
    COM = [sum(IMAGE.*gX, 'all'), sum(IMAGE.*gY, 'all'), sum(IMAGE.*gZ, 'all')] ./ IMAGEsum;
else
    [gX, gY] = ndgrid(vX, vY);
    COM = [sum(IMAGE.*gX, 'all'), sum(IMAGE.*gY, 'all')] ./ IMAGEsum;
end

if any(isnan(COM))
    error('EMC:IMAGE', 'IMAGE has a NaN and/or Inf value, resulting into a COM of %s', mat2str(COM))
end

end  % EMC_centerOfMass
