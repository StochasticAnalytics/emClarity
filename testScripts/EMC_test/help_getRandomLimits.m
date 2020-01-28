function [LIMITS] = help_getRandomLimits(REPEAT, RANGE, PROFILE, AXIS, SIDE, DIMENSIONS)
%
% NUMBER (int):     Number of batch to repeat
% RANGE (numeric):  [min, max]
% PROFILE (cll):    'both', 'pad' or 'crop'
% AXIS (cell):       'all', 'x', 'y', 'z', 'xy', 'xz' or 'yz'
% SIDE (cell):       'both', 'left' or 'right'
%

if strcmpi(DIMENSIONS, '3d')
    ndim = 6;
elseif strcmpi(DIMENSIONS, '2d')
    ndim = 4;
else
    error('number of dimensions not supported')
end

numberOfRuns = REPEAT * length(PROFILE) * length(AXIS) * length(SIDE);
LIMITS = cell(numberOfRuns, 1);

increment = 1;
for iSide = 1:length(SIDE)
    if strcmpi(SIDE{iSide}, 'left')
        if ndim == 6; maskSide = [1,0,1,0,1,0]; else; maskSide = [1,0,1,0]; end
    elseif strcmpi(SIDE, 'right')
        if ndim == 6; maskSide = [0,1,0,1,0,1]; else; maskSide = [0,1,0,1]; end
    else
        if ndim == 6; maskSide = [1,1,1,1,1,1]; else; maskSide = [1,1,1,1]; end
    end
            
    for iAxis = 1:length(AXIS)
       	if strcmpi(AXIS{iAxis}, 'x')
          	if ndim == 6; maskAxis = [1,1,0,0,0,0]; else; maskAxis = [1,1,0,0]; end
        elseif strcmpi(AXIS{iAxis}, 'y')
            if ndim == 6; maskAxis = [0,0,1,1,0,0]; else; maskAxis = [0,0,1,1]; end
        elseif strcmpi(AXIS{iAxis}, 'z')
            if ndim == 6; maskAxis = [0,0,0,0,1,1]; else; maskAxis = [1,1,1,1]; end
        elseif strcmpi(AXIS{iAxis}, 'xy')
            if ndim == 6; maskAxis = [1,1,1,1,0,0]; else; maskAxis = [1,1,1,1]; end
        elseif strcmpi(AXIS{iAxis}, 'xz')
            if ndim == 6; maskAxis = [1,1,0,0,1,1]; else; maskAxis = [1,1,1,1]; end
        elseif strcmpi(AXIS{iAxis}, 'yz')
            if ndim == 6; maskAxis = [0,0,1,1,1,1]; else; maskAxis = [1,1,1,1]; end
        else
            if ndim == 6; maskAxis = [1,1,1,1,1,1]; else; maskAxis = [1,1,1,1]; end
        end
            
        for iProfile = 1:length(PROFILE)
            
            lim = floor((RANGE(2)-RANGE(1)) .* rand(REPEAT, ndim) + RANGE(1));
            
            % profile
            if strcmpi(PROFILE{iProfile}, 'pad')
                lim = abs(lim);
            elseif strcmpi(PROFILE{iProfile}, 'crop')
                lim = -1 .* abs(lim);
            end
            
            lim = lim .* (maskAxis .* maskSide);
            
            LIMITS(increment:increment+REPEAT-1) = num2cell(lim, 2);
            increment = increment + REPEAT;
        end
    end
end

end  % end help_getRandomLimits