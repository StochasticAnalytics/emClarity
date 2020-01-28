function [is3d, SIZE, ndim] = EMC_is3d(SIZE)
%
% [is3d, SIZE, ndim] = EMC_is3d(SIZE)
% Check whether or not SIZE is describing the size of a 3d array.
%
% Input:
%   SIZE (int vector):  Vector describing the size of an array; 3d:[x,y,z], 2d:[x,y] or 1d:[x,1]|[1,y].
%                       NOTE: Column vectors are also accepted.
%                       NOTE: Although 2d unit vectors ([1, 1]) are describing scalars, they are accepted.
%
% Output:
%   is3d (bool):        Wheter or not SIZE is describing the size of a 3d array.
%
%   SIZE (int vector):  z squeezed SIZE ([x,y,1] -> [x,y])
%
%   ndim (int):         Number of element in SIZE. Can be 2 or 3.
%                       NOTE: As with ndims(SIZE), if z=1 -> ndim=2
%
% Example:
%   [is3d, SIZE, ndim] = EMC_is3d([10,10,10]);  output: is3d=true,	SIZE=[10,10,10], ndim=3
%   [is3d, SIZE, ndim] = EMC_is3d([1,10,10]);   output: is3d=true,	SIZE=[10,1,10],  ndim=3
%   [is3d, SIZE, ndim] = EMC_is3d([10,10,1]);   output: is3d=false,	SIZE=[10,10],    ndim=2
%   [is3d, SIZE, ndim] = EMC_is3d([1,10,1]);    output: is3d=false,	SIZE=[1,10],     ndim=2
%   [is3d, SIZE, ndim] = EMC_is3d([1,10]);      output: is3d=false,	SIZE=[1,10],     ndim=2
%   [is3d, SIZE, ndim] = EMC_is3d([10,1]);      output: is3d=false,	SIZE=[10,1],     ndim=2
%

% Created:  18Jan2020
% Version:  v.1.0   better docstring (TF, 23Jan20).
%           v.1.1   unittest (TF, 23Jan20).
%

if ~isnumeric(SIZE) || ~isvector(SIZE) || any(SIZE < 1) || any(rem(SIZE, 1))
    error('EMC_is3d:SIZE', 'SIZE should be a numeric vector of positive integers')
end

ndim = numel(SIZE);
if ndim == 3
	if SIZE(3) == 1
     	SIZE = SIZE(1:2);
       	is3d = false;
        ndim = 2;
	else
      	is3d = true;
	end
elseif ndim == 2
  	is3d = false;
elseif ndim == 1
    error('EMC_is3d:SIZE', ['SIZE should have at least 2 elements. To describe vectors, the SIZE should ', ...
          'be [1, N] for row vectors or [N, 1] for column vector.']);
else
  	error('EMC_is3d:SIZE', 'SIZE has more element than maximum supported (3), got %d', ndim);
end

end  % EMC_is3d
