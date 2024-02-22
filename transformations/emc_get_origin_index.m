function [ origin_index ] = emc_get_origin_index( obj )

    if isempty(obj)
        error('Input must not be empty');
    end
    if ~isnumeric(obj)
        error('Input must be a numeric object');
    end
    if any(obj <= 0)
        error('Input must be a positive number');
    end
    
    origin_index = floor(obj/2) + 1;
end