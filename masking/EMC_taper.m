function [TAPER] = EMC_taper(TYPE, FIRST, LAST, SIZE)
%
% Compute a taper of a given type and a given length.
%
% TYPE (str):       Type of the taper; 'cosine', 'linear' or 'cosine'.
%
% START (float):    First value of the taper.
%
% END (float):      Last value of the taper.
%
% SIZE (int):       Length (in pixel) of the taper.
%
%-------
% NOTE:             The tapers do no include the first and last values.
%
%-------
% TODO:             Add gaussian?
%

if strcmpi(TYPE, 'cosine')
    TAPER = cos(((1:SIZE) * pi) ./ SIZE) .* (FIRST-LAST)/2 + abs(FIRST-LAST)/2 + min(FIRST, LAST);
elseif strcmpi(TYPE, 'linear')
    TAPER = linspace(FIRST, LAST, SIZE);
    TAPER = TAPER(2:end);
else
    error("TYPE should be 'cosine' or 'linear', got %s", TYPE)
end
