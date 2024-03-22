function [output_vol] =  emc_halfcast(input_vol, swap_host_device)

    if nargin < 2
        swap_host_device = false;
    end

    to_gpu = false;
    to_cpu = false;
    to_half = false;
    
    % Determin if we are going to or from half based on the input precision.
    switch underlyingType(input_vol)
        case 'uint16'
            to_half = false;
        case 'single'
            to_half = true;
        otherwise
            error('Unknown precision');
    end

    % By default, stay on cpu or gpu and just convert type.
    % If swap_host_device is true, then we will swap to the other device.
    if (swap_host_device)
        if isa(input_vol, 'gpuArray')
            to_cpu = true;
        else
            to_gpu = true;
        end
    else
        if isa(input_vol, 'gpuArray')
            to_gpu = true;
        else
            to_cpu = true;
        end
    end

    if (to_cpu && to_gpu)
        error('Cannot convert to and from GPU at the same time');
    end

    % fprintf("To half: %d\n", to_half);
    % fprintf("To GPU: %d\n", to_gpu);
    % fprintf("To CPU: %d\n", to_cpu);
    % fprintf("Swap: %d\n", swap_host_device);

    if (to_half)
        if (to_gpu)
            output_vol = zeros(size(input_vol), 'uint16', 'gpuArray');
        else
            output_vol = zeros(size(input_vol), 'uint16');
        end
        mexFP16(input_vol, output_vol, to_half, int64(numel(input_vol)));
    else 
        if (to_gpu)
            output_vol = zeros(size(input_vol), 'single', 'gpuArray')+2;
        else
            output_vol = zeros(size(input_vol), 'single') +1;
        end
        mexFP16(output_vol, input_vol, to_half, int64(numel(input_vol)));
    end

end