function [ output_int ] = ReturnMaxPositionInStack(obj,exclude_negative_film_numbers)
%   Detailed explanation goes here


	output_int = intmin();

	for line = 0:length(all_parameters)
	
		if (obj.ReturnImageIsActive(line) >= 0 || ~exclude_negative_film_numbers)
			temp_int = obj.ReturnPositionInStack(line);
			if (output_int < temp_int) 
        output_int = temp_int;
      end
    end
    
  end
  
end