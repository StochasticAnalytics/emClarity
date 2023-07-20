function [ output_double ] = EMC_str2double( input_str )
    %Make a local copy of the default cluster, and modify the job storage
    %location
    
      output_double = str2double(input_str);
      if (isnan(output_double))
        output_double = str2num(input_str);
        isempty(output_double)
        if (isempty(output_double))   
        
          err_msg = sprintf('EMC_str2double: input string is not a number!\nReceived: %s', input_str);
          error(err_msg);
        end
      end
end