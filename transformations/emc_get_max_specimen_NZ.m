function [ max_specimen_nz, tomoNumber ] = emc_get_max_specimen_NZ(subTomoMeta_tomoName, subTomoMeta_coords, tomo_name_list, n_tomograms, samplingRate)
 
    max_z_value = -inf;
    min_z_value = inf;

    tomoNumber = zeros(n_tomograms,1);

    for iTomo = 1:n_tomograms
        % subTomoMeta_tomoName = subTomoMeta.mapBackGeometry.tomoName
        % subTomoMeta_coords = subTomoMeta.mapBackGeometry.(tiltName).coords
        if isa(subTomoMeta_tomoName,'struct')
            tomoNumber(iTomo) = subTomoMeta_tomoName.(tomo_name_list{iTomo}).tomoNumber;
        else
            tomoNumber(iTomo) = iTomo;
        end
    
        % 4 is the unbinned pixel size of the tomogram in Z
        % 6 is location of the specimen origin in Z relative to the origin of the tomogram
        nZ = subTomoMeta_coords(tomoNumber(iTomo),4) ./ samplingRate;
        oZ = subTomoMeta_coords(tomoNumber(iTomo),6) ./ samplingRate;
        
        % We need to consider the shift of the tomogram relative to the specimen origin,
        % which is -oZ
        dZ = -oZ;

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