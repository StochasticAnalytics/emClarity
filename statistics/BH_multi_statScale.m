function [IMG] = BH_multi_statScale(IMG,DATATYPE)
%Scale to zero mean, R variance, and cast as data type
%   Detailed explanation goes here
 

m = mean(IMG(:));

if abs(m) > 1e-4
  % set mean to zero
  IMG = IMG - m;
end

r = rms(IMG(:));

if abs(r-1) > 1e-4
  % set the rms to the desired val
  IMG = IMG ./ (r);
end

% number of standard deviations above which data values are compressed
% (rather than truncated)
maxNonCompressed = 4;
      
if ~strcmp(class(IMG), DATATYPE)

  % This should only be used with "raw" data like image stacks or tomograms
  % where no individual set of pixels should have extremely large values. 
  switch DATATYPE
    case 'int16'

      compressMask = IMG > maxNonCompressed;
      IMG(compressMask) = maxNonCompressed + ...
                            (sqrt(single(IMG(compressMask))) - sqrt(maxNonCompressed));
      
      compressMask = IMG < -1*maxNonCompressed;
      IMG(compressMask) = maxNonCompressed - ...
                            (sqrt(abs(single(IMG(compressMask)))) - sqrt(maxNonCompressed));                          
      
      IMG = int16(IMG .* (2^16 / (2*max(abs(IMG(:))))));
     
    case 'uint16'
 
      compressMask = IMG > maxNonCompressed;
      IMG(compressMask) = maxNonCompressed + ...
                            (sqrt(IMG(compressMask)) - sqrt(maxNonCompressed));
      
      compressMask = IMG < -1*maxNonCompressed;
      IMG(compressMask) = maxNonCompressed - ...
                            (sqrt(abs(IMG(compressMask))) - sqrt(maxNonCompressed));  
                          
      IMG = IMG - min(IMG(:));
      IMG = uint16(IMG .* (2^16 / (max(abs(IMG(:))))));
      
    otherwise
      error('Only set up for int16 and uint16');
  end
end

end

