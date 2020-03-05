function INDEX = EMC_maskIndex(TYPE, SIZE, METHOD, OPTION)
%
% INDEX = EMC_maskIndex(TYPE, SIZE, METHOD, OPTION)
% Compute linear indexes for different TYPE of wrapping.
%
% Input:
%   TYPE (str):             Type of wrapping available.
%                           'fftshift':  Calculate the linear indices to go from a not-centered
%                                        (zero first) to a centered spectrum.
%
%                           'ifftshift': Calculate the linear indices to go from a centered
%                                        to a not-centered (zero first) spectrum.
%
%                           'half2full': Calculate the linear indices to go from a half grid to
%                                        full grid. This is for not-centered (zero first) grids.
%                                        These indexes are meant to be applied on the full grid
%                                        containing half of the spectrum.
%                                        See EMC_irfftn for more details.
%
%   SIZE (vector):          Size (in pixel) of the 3d/2d grid; [x, y, z] or [x, y].
%                           NOTE: [1,1], [N,1] or [1,N] are not allowed.
%                           NOTE: It must correspond to the full grid size even if OPTION.half = true.
%
%   METHOD (str):           Method used to compute the INDEX grid; 'cpu' or 'gpu'.
%
%   OPTION (cell|struct):   Optional parameters.
%                           If cell: {field,value ; ...}, note the ';' between parameters.
%                           NOTE: Can be empty.
%                           NOTE: Unknown fields will raise an error.
%
%     -> 'half' (bool):     Whether or not the output mask should correspond to the half (non-redondant)
%                           grid. If true, the first dimension of the output will have floor(SIZE(1)/2)+1
%                           pixels.
%                           NOTE: It must be false if TYPE = 'half2full'.
%                           default = false
%
%     -> 'precision' (str): Precision of the INDEX grid.
%                           float: 'single' or 'double'.
%                           integer: 'int', 'uint'.
%                           default = 'uint' (select the most appropriate unassigned integer)
%
% Output:
%   INDEX (numeric):        Linear indexes.
%                           NOTE: If OPTION.half=true, size(INDEX,1) == floor(SIZE(1)/2) + 1.
%
% Note:
% - 2d: sub2ind(SIZE, gX, gY) <=> gX + (gY-1)*SIZE(1)
%   3d: sub2ind(SIZE, gX, gY, gZ) <=> gX + (gY-1)*SIZE(1) + (gZ-1)*SIZE(1)*SIZE(2)
%   This is already much faster than sub2ind and it can be further simplified using vectors.
%   Moreover, it can be done directly with integers, whereas sub2ind requires floats.
%
% Other EMC-files required:
%   EMC_is3d, EMC_getOption, EMC_setPrecision, EMC_setMethod
%
% See also fftshift, ifftshift, EMC_rfftn, EMC_irfftn
%

% Created:  15Feb2020, R2019a
% Version:  v.1.0   Unittest (TF, 16Feb2020).
%

%% checkIN
[is3d, SIZE] = EMC_is3d(SIZE);
if any(SIZE == 1)
    error('EMC:SIZE', 'SIZE should corresponding to a 2d or 3d array, got %s', mat2str(SIZE))
end

OPTION = EMC_getOption(OPTION, {'half', 'precision'}, false);

% half
if isfield(OPTION, 'half')
    if ~islogical(OPTION.half) || ~isscalar(OPTION.half)
        error('EMC:half', 'OPTION.half should be a boolean, got %s', class(OPTION.half));
    elseif OPTION.half
        SIZE(1) = floor(SIZE(1)/2) + 1;  % switch to corresponding half size
    end
else
    OPTION.half = false;  % default
end

% precision
if isfield(OPTION, 'precision')
    if (ischar(OPTION.precision) || isstring(OPTION.precision))
        if strcmpi(OPTION.precision, 'int')
            if prod(SIZE) <= 2^16
                OPTION.precision = 'int16';
            elseif prod(SIZE) <= 2^32
                OPTION.precision = 'int32';
            else
                OPTION.precision = 'int64';
            end
        elseif strcmpi(OPTION.precision, 'uint')
            if prod(SIZE) <= 2^16
                OPTION.precision = 'uint16';
            elseif prod(SIZE) <= 2^32
                OPTION.precision = 'uint32';
            else
                OPTION.precision = 'uint64';
            end
        elseif ~strcmpi(OPTION.precision, 'single') && ~strcmpi(OPTION.precision, 'double')
            error('EMC:precision', "OPTION.precision should be 'single', 'double', 'int' or 'uint', got %s", OPTION.precision)
        end
    else
      	error('EMC:precision', "OPTION.precision should be a string|char vector, got %s", class(OPTION.precision))
    end
elseif prod(SIZE) <= 2^16
    OPTION.precision = 'uint16';  % default
elseif prod(SIZE) <= 2^32
    OPTION.precision = 'uint32';  % default
else
    OPTION.precision = 'uint64';  % default
end

%% Compute the indexes
% half2full
if strcmpi(TYPE, 'half2full')
    cR = floor(SIZE/2) + 1;  % center receiver
    cD = ceil(SIZE/2);  % center donor

    SIZE = EMC_setMethod(EMC_setPrecision(SIZE, OPTION.precision), METHOD);
    INDEX = reshape(1:prod(SIZE, 'native'), SIZE);  % linear indexes of a grid with desired size

    if is3d
        INDEX(cR(1)+1:end, 1,           1)           = INDEX(cD(1):-1:2, 1,              1);
        INDEX(cR(1)+1:end, 2:cR(2),     1)           = INDEX(cD(1):-1:2, end:-1:cD(2)+1, 1);
        INDEX(cR(1)+1:end, cR(2)+1:end, 1)           = INDEX(cD(1):-1:2, cD(2):-1:2,     1);
        INDEX(cR(1)+1:end, 1,           2:cR(3))     = INDEX(cD(1):-1:2, 1,              end:-1:cD(3)+1);
        INDEX(cR(1)+1:end, 1,           cR(3)+1:end) = INDEX(cD(1):-1:2, 1,              cD(3):-1:2);
        INDEX(cR(1)+1:end, 2:cR(2),     2:cR(3))     = INDEX(cD(1):-1:2, end:-1:cD(2)+1, end:-1:cD(3)+1);
        INDEX(cR(1)+1:end, 2:cR(2),     cR(3)+1:end) = INDEX(cD(1):-1:2, end:-1:cD(2)+1, cD(3):-1:2);
        INDEX(cR(1)+1:end, cR(2)+1:end, 2:cR(3))     = INDEX(cD(1):-1:2, cD(2):-1:2,     end:-1:cD(3)+1);
        INDEX(cR(1)+1:end, cR(2)+1:end, cR(3)+1:end) = INDEX(cD(1):-1:2, cD(2):-1:2,     cD(3):-1:2);
    else  % 2d
        INDEX(cR(1)+1:end, 1)           = INDEX(cD(1):-1:2, 1);
        INDEX(cR(1)+1:end, 2:cR(2))     = INDEX(cD(1):-1:2, end:-1:cD(2)+1);
        INDEX(cR(1)+1:end, cR(2)+1:end) = INDEX(cD(1):-1:2, cD(2):-1:2);
    end

% fftshift and ifftshift
else
    if strcmpi(TYPE, 'fftshift')
        half = ceil(SIZE/2);
    elseif strcmpi(TYPE, 'ifftshift')
        half = floor(SIZE/2);
    else
        error('EMC:TYPE', "TYPE should be 'fftshift', 'ifftshift' or 'half2full'")
    end

    SIZE = EMC_setMethod(EMC_setPrecision(SIZE, OPTION.precision), METHOD);  % convert after the division
    if OPTION.half; vX = (1:SIZE(1))'; else; vX = [half(1)+1:SIZE(1), 1:half(1)]'; end

    % Concatenation appears to be faster than fftshift and circshift.
    if is3d
        INDEX = vX + [half(2):SIZE(2)-1, 0:half(2)-1] * SIZE(1) + ...
                     reshape([half(3):SIZE(3)-1, 0:half(3)-1],1,1,[]) * (SIZE(1) * SIZE(2));
    else
        INDEX = vX + [half(2):SIZE(2)-1, 0:half(2)-1] * SIZE(1);
    end
end

end  % EMC_coordIndex
