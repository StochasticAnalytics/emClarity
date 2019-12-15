function [TAPER] = EMC_multi_taper(TYPE, START, END, LENGTH)
%
% Compute a taper of a given type and a given length.
%
% TYPE (str):       Type of the taper; 'cosine' or 'linear'.
%
% START (float):    First value of the taper.
%
% END (float):      Last value of the taper.
%
% LENGTH (int):     Length (in pixel) of the taper.
%
%-------
% TODO:             Add gaussian.

if strcmpi(TYPE, 'cosine')
    TAPER = cos(((0:LENGTH-1) * pi) ./ (LENGTH-1)) .* (START-END)/2 + abs(START-END)/2 + min(START, END);
elseif strcmpi(TYPE, 'linear')
    TAPER = linspace(START, END, LENGTH);
else
    error("TYPE should be 'cosine' or 'linear', got %s", TYPE)
end