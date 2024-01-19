function [ KERNEL ] = BH_multi_gaussian3d(SIZE, StdDev )
%Create a 3d gaussian kernel for convolution based smoothing in real space.
%   Detailed explanation goes here

if numel(SIZE) == 1
  hood = [SIZE,SIZE,SIZE];
elseif numel(SIZE) ~= 3
  error('SIZE must be 1 or 3 vector');
else
  hood = SIZE;
end

if StdDev < 0
  StdDev = abs(StdDev);
  fourthPow = 2;
else
  fourthPow = 1;
end
[ X,Y,Z,~,~,~ ] = BH_multi_gridCoordinates( hood, 'Cartesian', 'cpu', ...
  {'none'}, 0, 1, 0 );
% Zeroth order gaussian derivative
gfunc = @(x,y,z,s)((2.*pi.*s.^2*fourthPow).^(-3/2).*...
  exp(-0.5*(x./s).^(2*fourthPow)).*...
  exp(-0.5*(y./s).^(2*fourthPow)).* ...
  exp(-0.5*(z./s).^(2*fourthPow)));


KERNEL = gfunc(X,Y,Z,StdDev);
% normalize
KERNEL = KERNEL .* sum(KERNEL(:)).^-1;
clear X Y Z x y z hood c
end

