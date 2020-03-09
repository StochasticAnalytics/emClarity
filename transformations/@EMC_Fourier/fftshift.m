function DFT = fftshift(DFT)
%
% DFT = fftshift(DFT)
% Shift zero-frequency component to the center of spectrum.
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
if strcmpi(TYPE, 'nc2nc')
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

elseif strcmpi(TYPE, 'c2c')
    cD = floor(SIZE/2) + 1;  % center donor
    cR = cD;  % center receiver
    if mod(SIZE(1),2); isEven = 1; else; isEven = 0; end 

    SIZE = EMC_setMethod(EMC_setPrecision(SIZE, OPTION.precision), METHOD);
    INDEX = reshape(1:prod(SIZE, 'native'), SIZE);  % linear indexes of a grid with desired size

    if is3d
        
        
    else
        INDEX(cR(1):end, :) = INDEX(1:cD(1)-isEven, :);
        INDEX(1:cR(1)-1, 1:cR(2)-1) = INDEX(cD(1)+1:-1:2, cD(2):end);
        INDEX(1:cR(1)-1, cR(2):end) = INDEX(cD(1)+1:-1:2, 1:cD(2)-1);
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
