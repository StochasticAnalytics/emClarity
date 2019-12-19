function [ ] = SetAllReference3DFilename(obj,wanted_filename)
%   Detailed explanation goes here


	for  counter = 1:length(obj.all_parameters)
		obj.all_parameters{counter}.reference_3d_filename = wanted_filename;
  end

end