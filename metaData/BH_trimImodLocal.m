function [outputArg1,outputArg2] = BH_trimImodLocal(inputArg1,inputArg2)
%Remove tilts from an imod local alignment file
%   Detailed explanation goes here

% From the inputs
% Get a filename for the .local .rawtlt
% Get an array for tilts to remove

% Read count of starting tilt numbers using len()
% Use text scan to get one array of all alignments inlcuding header
% f = fopen('tilt1.local_patchTracking','r');
% p = textscan(f,'%f','TreatAsEmpty',{' '}, 'Delimiter','\n')

% First nine entries are header
%header = p{1}(1:9)
% Total number of blocks is header1 * header2


% 7 and 8 in the header are booleans, if xtilt and zfactor
% So remove header
% then there are nTilts * (1 + ifXtilt + ifZfactor) entries that are per
% tilt
% Then there are nTilts with 6 columsn for the XF

% loop through, create some kind of backup, and save the new alignments.



end

