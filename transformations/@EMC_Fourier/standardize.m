function DFT = standardize(obj, DFT)
%
% DFT = obj.standardize(DFT)
% Set the real space mean to ~0 and real space variance to ~1 of the ARRAY.
% The mean and variance are changed in frequency space.
%
% Input:
%   DFT (numeric):  2d/3d Discrete Fourier Transform.
%
% Output:
%   DFT (numeric):  Standardize DFT.
%

if ~obj.half
    DFT(1) = 0;
   	DFT = DFT ./ (sqrt(sum(abs(DFT).^2, 'all')) / numel(DFT));
if obj.half
    % Capture every independant chunk (same for 2d or 3d)
    if obj.centered
        error('not finished')
        
    else
        DFT(1) = 0;
        cD = ceil(obj.size_real(1)/2);  % center of donor
        factor = sum(abs(DFT(1,:,:)).^2, 'all');  % unique row/plane
        factor = factor + 2*sum(abs(DFT(2:cD,:,:)).^2, 'all');  % common chunk
        if obj.isEven(1); factor = factor + sum(abs(DFT(cD+1,:,:)).^2, 'all'); end  % unique row/plane
    end

  	DFT = DFT ./ (sqrt(factor) / prod(obj.size_real, 'native'));
end

end