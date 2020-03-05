function NUM = EMC_setPrecision(NUM, PRECISION)
%
% Cast the numerical array NUM to the desired PRECISION.
%
if strcmp(PRECISION, 'single')
	NUM = single(NUM);
elseif strcmp(PRECISION, 'double')
	NUM = double(NUM);
else
	error('EMC:precision', "PRECISION should be 'single' or 'double'")
end

end
