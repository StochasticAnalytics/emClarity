function NUM = EMC_setPrecision(NUM, PRECISION)
%
% Cast the numerical array NUM to the desired PRECISION.
% This is similar to cast, but it's faster for single/double
% and allows raison EMC:precision.
%
switch PRECISION
    case 'single'
        NUM = single(NUM);
    case 'double'
        NUM = double(NUM);
    case 'int32'
        NUM = int32(NUM);
    case 'uint32'
        NUM = uint32(NUM);
    case 'logical'
        NUM = logical(NUM);
    case 'int16'
        NUM = int16(NUM);
    case 'uint16'
        NUM = uint16(NUM);
    case 'int64'
        NUM = int64(NUM);
    case 'uint64'
        NUM = uint64(NUM);
    otherwise
        error('EMC:precision', ["PRECISION should be 'single', 'double', 'logical', ", ...
                                "'int16', 'int32', 'int64', 'uint16', 'uint32' or 'uint64'"])
end

end
