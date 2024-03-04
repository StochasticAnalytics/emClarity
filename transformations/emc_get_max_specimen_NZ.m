function [ max_specimen_nz ] = emc_get_max_specimen_NZ(subTomoMeta_coords, tomo_name_list, n_tomograms, samplingRate)
 
    max_z_value = -inf;
    min_z_value = inf;


    for iTomo = 1:n_tomograms    
        % 4 is the unbinned pixel size of the tomogram in Z
        % 6 is location of the specimen origin in Z relative to the origin of the tomogram
        if isa(subTomoMeta_coords,'cell')
            nZ = subTomoMeta_coords{iTomo}.NZ;
            dZ = subTomoMeta_coords{iTomo}.dZ_specimen_to_tomo;
        elseif isa(subTomoMeta_coords,'struct')
            nZ = subTomoMeta_coords.(tomo_name_list{iTomo}).NZ;
            dZ = subTomoMeta_coords.(tomo_name_list{iTomo}).dZ_specimen_to_tomo;
        else
            error('The field NZ or the field %s is not present in the subTomoMeta_coords', tomo_name_list{iTomo});
        end

        nZ = nZ ./ samplingRate;
        dZ = dZ ./ samplingRate;

        if (dZ + nZ / 2 > max_z_value)
            max_z_value = dZ + nZ / 2;
        end
        if (dZ - nZ / 2 < min_z_value)
            min_z_value = dZ - nZ / 2;
        end
    end

    if (max_z_value < min_z_value)
        error('The max z value is less than the min z value');
    end

    max_specimen_nz = ceil(max_z_value - min_z_value + (samplingRate * 2));

end