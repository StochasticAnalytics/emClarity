function [ is_valid ] = emc_check_for_valid_image_file(wanted_filename)

    % Check if the file exists
    if isfile(wanted_filename)
        try 
            test_header = MRCImage(wanted_filename, 0);
            is_valid = true;
            return;
        catch
            system(['rm ' wanted_filename]);
            is_valid = false;
            return;
        end
    else
        is_valid = false;
        return;
    end

end