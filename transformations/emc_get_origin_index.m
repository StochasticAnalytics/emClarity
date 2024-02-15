function [ origin ] = emc_get_origin_index( obj )

    if isempty(obj)
        error('Input must not be empty');
    end
    if ~isnumeric(obj)
        error('Input must be a numeric object');
    end

    origin = floor(obj/2) + 1;

end