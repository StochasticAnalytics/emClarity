function [output_vol] =  emc_halfcast(input_vol)

    if isa(input_vol, 'uint16')
        output_vol = zeros(size(input_vol), 'single');
        mexFP16(output_vol, input_vol, false);
    elseif isa(input_vol, 'single')
        output_vol = zeros(size(input_vol), 'uint16');
        mexFP16(input_vol, output_vol, true);
    else
        error('Unknown precision');
    end

end