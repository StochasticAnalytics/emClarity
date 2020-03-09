function FULL = half2full(WRAP, HALF, SIZE, varargin)
%
% FULL = EMC_Fourier.half2full(WRAP, HALF, SIZE, varargin)
% Applied the Hermitian symetry to the HALF spectrum to reconstruct the FULL spectrum.
%
% Input:
%   WRAP (str|char):        Type of wrapping (always half -> full).
%                           'nc2nc': Calculate the linear indices to go from a half
%                                    non-centered spectrum to a full non-centered spectrum.
%                           'c2nc':  Calculate the linear indices to go from a half centered
%                                    spectrum to a full non-centered spectrum.
%                           'c2c':   Calculate the linear indices to go from a half centered
%                                    spectrum to a full centered spectrum.
%
%   HALF (numeric):         Half spectrum of size floor(SIZE/2)+1.
%                           NOTE: if real (not complex), the Hermitian symetry becomes a central symetry.
%
%   SIZE (vector):          Size (in pixel) of the full spectrum; [x, y, z] or [x, y].
%                           NOTE: [1,1], [N,1] or [1,N] are not allowed.
%
%   (optional)
%   varargin{1} (numeric):  Use this wrapper to reconstruct the full grid.
%                         	Its size and method should correspond to SIZE
%                         	and method of HALF, respectively.
%
% Output:
%   FULL (numeric):     Full spectrum of size obj.size_real.
%                       Method and precision are unchanged.
%
% Example:
%   - Compute the full (redundant) spectrum using the half (non-redundant) spectrum:
%     >> ft = EMC_Fourier([128,128], 'cpu', {});
%     >> dft_half = ft.fft(rand(128,128,'single'));
%     >> dft_full = EMC_Fourier.half2full('nc2nc', dft_half, [128,128]);
%
% Other EMC-files required:
%   EMC_is3d, EMC_getClass, EMC_shareMethod, EMC_Fourier.getIndex
%

% Created: 8Mar2020, R2019a
% Version: 1.0.
%

[~, SIZE] = EMC_is3d(SIZE);
[precision, isOnGpu, method] = EMC_getClass(HALF);
if isOnGpu
    FULL = zeros(SIZE, precision, 'gpuArray');
else
    FULL = zeros(SIZE, precision);
end

oX = floor(size(SIZE,1)/2)+1;

% Every wrapper is expecting the half grid in this position.
FULL(1:oX+1, :, :) = HALF;

if nargin == 4
    INDEX = varargin{1};
    if ~isequal(size(INDEX), SIZE) || ~EMC_shareMethod(INDEX, HALF)
        error('EMC:Fourier', 'INDEX should have correspond to SIZE and have the same method as HALF')
    end
else
    INDEX = EMC_Fourier.getIndex(WRAP, SIZE, method, {});
end

FULL = FULL(INDEX);  % wrap

% If complex, apply conjugate
if ~isreal(HALF)
    if strcmpi(WRAP, 'nc2nc') || strcmpi(WRAP, 'c2nc')
        FULL(oX(1)+1:end, :, :) = conj(FULL(oX(1)+1:end, :, :));
    elseif strcmpi(WRAP, 'c2c')  % half centered -> full centered
        if ~mod(SIZE(1),2)
            FULL(2:oX(1), :, :) = conj(FULL(2:oX(1), :, :));
        else
            FULL(1:oX(1), :, :) = conj(FULL(1:oX(1), :, :));
        end
    else
        error('EMC:Fourier', "WRAP should be 'n2c', 'c2nc' or 'c2c'")
    end
end

end  % half2full
