function [Kx, Ky, Kxy ] = BH_multi_gaussian2d(SIZE, StdDev, varargin )
%Create a 2d gaussian kernel for convolution based smoothing in real space.
%   Note that the real space versions are normalized to have an area of 1
%   under the full kernel (or they should be anyway) while the fourier
%   space are not normalized so that the max transfer for a 0 order is 1 at
%   the origin.

Kx = '';
Ky = '';
Kxy = '';

flgFourier = 0;
if any(SIZE < 0)
  SIZE = abs(SIZE);
  flgFourier = 1;
end

if numel(SIZE) == 1
  hood = [SIZE,SIZE,1];
elseif numel(SIZE) ~= 2
  error('SIZE must be 1 or 2 vector');
else
  hood = SIZE;
end


if nargin > 2
  ORDER = varargin{1};
else
  ORDER = 0;
end

% If real space, shift origin to center, if fourier leave origin at corner.
[ X,Y,~,~,~,~ ] = BH_multi_gridCoordinates( hood, 'Cartesian', 'cpu', ...
  {'none'}, flgFourier, ...
  1-flgFourier, 0 );




if (flgFourier)
  % Zeroth order gaussian derivative
  gfunc = @(x,y,s)(exp(-0.5*(x.^2+y.^2).*(s.*2.*pi).^(2)));
  normFactor = @(x,y,r,Ord)((2i.*pi).^Ord);
else
  gfunc = @(x,y,s)(exp(-0.5*(x.^2+y.^2)./s.^2));
  normFactor = @(x,y,s,Ord)((2.*pi).^-1.*(s^2).^(-(Ord+1)));
end


switch ORDER
  
  case 0
    
    if ( flgFourier )
      Kx = normFactor(X,Y,1,ORDER).*gfunc(X,Y,StdDev);
    else
      Kx = normFactor(X,Y,StdDev,ORDER).*gfunc(X,Y,StdDev);
    end
    
  case 1
    
    if ( flgFourier )
      Kx = -X.*normFactor(X,Y,X,ORDER).*gfunc(X,Y,StdDev);
      Ky = -Y.*normFactor(X,Y,Y,ORDER).*gfunc(X,Y,StdDev);
    else
      Kx = normFactor(X,Y,StdDev,ORDER).*gfunc(X,Y,StdDev);
      Ky = normFactor(X,Y,StdDev,ORDER).*gfunc(X,Y,StdDev);
    end
    
  case 2
    
    if ( flgFourier )
      Kx  = X.^2.*normFactor(X,Y,X,ORDER).*gfunc(X,Y,StdDev);
      Ky  = Y.^2.*normFactor(X,Y,Y,ORDER).*gfunc(X,Y,StdDev);
      Kxy = X.*Y.*normFactor(X,Y,(X.*Y),ORDER).*gfunc(X,Y,StdDev);
    else
      Kx  = (-1+(X./StdDev).^2).*normFactor(X,Y,StdDev,ORDER).*gfunc(X,Y,StdDev);
      Ky  = (-1+(Y./StdDev).^2).*normFactor(X,Y,StdDev,ORDER).*gfunc(X,Y,StdDev);
      Kxy = X.*Y.*normFactor(X,Y,StdDev,ORDER).*gfunc(X,Y,StdDev);
    end
    
end

% if (flgFourier)
%   KERNELX = gfunc(X,Y,StdDev); %KERNELX = KERNELX ./ max(KERNELX(:));
%   %KERNELY = gfunc(Y,StdDev); %KERNELY = KERNELY ./ max(KERNELY(:));
% else
%   KERNELX = gfunc(X,Y,StdDev);
% end
% % normalize
%KERNEL = KERNEL .* sum(KERNEL(:)).^-1;
clear X Y Z x y z hood c
end