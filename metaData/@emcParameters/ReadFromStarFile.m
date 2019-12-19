function [ ] = ReadFromStarFile(obj,wanted_filename, program_name, exclude_negative_film_numbers )
%   Detailed explanation goes here


% 	temp_line = emcParameterLine();
	obj.all_parameters.Clear();
	star_reader = emcStarFileReader(wanted_filename, obj.all_parameters, exclude_negative_film_numbers);
	obj.parameters_that_were_read = star_reader.parameters_that_were_read;
  
end

