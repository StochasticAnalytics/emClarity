function [ padVal ] = BH_multi_padVal( size1, size2 )
%Calc correct padding to keep origin at ceil(N+1)./2)
%   If the difference is even, pad the same on both ends. If it is odd, the
%   "extra" position should preced the origin for an odd dimension and follow
%   the origin for an even volume.


  sizeDiff = size2 - size1;
  oddDiff = mod(sizeDiff, 2);
  oddInput= mod(size1, 2);
  
  padVal = [ floor(sizeDiff./2) +  (oddInput) .* oddDiff ; ... 
             floor(sizeDiff./2) + ~(oddInput) .* oddDiff];
  
end

