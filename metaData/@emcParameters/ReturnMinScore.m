function [ output_float ] = ReturnMaxScore(obj,exclude_negative_film_numbers)
%   Detailed explanation goes here



	output_float = floatmax('single');

	for line = 0:length(all_parameters)
	
		if (obj.ReturnImageIsActive(line) >= 0 || ~exclude_negative_film_numbers)
			temp_int = obj.ReturnScore(line);
			if (output_float > temp_int) 
        output_float = temp_int;
      end
    end
    
  end
  
end