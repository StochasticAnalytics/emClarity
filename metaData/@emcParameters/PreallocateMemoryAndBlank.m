function [ ] = PreallocateMemoryAndBlank(obj,number_to_allocate)
%   Detailed explanation goes here

	obj.ClearAll();
	temp_line = cisTEMParameterLine();
	obj.all_parameters.Add(temp_line, number_to_allocate);
  
end