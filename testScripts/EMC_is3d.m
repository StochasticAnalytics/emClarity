function [flg3d, ndim] = EMC_is3d(SIZE)
%
% Check if the SIZE vector is describing a 2d or 3d IMAGE.
% If the z dimension is 1, it is considered 2d.
%
if numel(SIZE) == 3
    if SIZE(3) == 1
        flg3d = false;
        ndim = 2;
    else
        flg3d = true;
        ndim = 3;
    end
elseif numel(SIZE) == 2
    flg3d = false;
    ndim = 2;
else
    error('Only 2D or 3D masks are supported, got %dD', numel(SIZE));
end

end