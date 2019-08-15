function [ INDICES, PADVALUES, SHIFTS ] = ...
                            BH_isWindowValid( VOLUME_SIZE, WINDOW_SIZE, MASK_RADIUS, CENTER )
%Address out of bounds conditions.
%  
%
%   Input Variables:
%
%   VOLUMES_SIZE = Size of the volume to extract from.
% 
%   WINDOW_SIZE  = Size of the window to be extracted.
%
%   CENTER = Coordinates of the center of the window in the tomogram frame.
%
%
%   Output Variables:
%
%   INDICES = 2,3 array with [x1, y1, z1 ; x2, y2, z2]
%
%   PADVALUES = 2,3 array with pre and post padding values, same order.
% 
%   SHIFTS = fractional shifts
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%   Goals & Limitations:
%
%    ShiftVAL describes the shift from the window center to the particle
%    center, which is a result of the windowing process. XCF finds a shift
%    relative to the center (N+1)/2 of the window. This means any estimated
%    peak shifts should = + shiftVAL, and at the end of XCF this value
%    should be subtracted so the peak in the original volume coordinates is
%    returned.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%   TODO:
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  
minSizeMask = (max(MASK_RADIUS)+6).*[2,2,2];
winLowCorner = ceil((WINDOW_SIZE-1) ./ 2);
% if window size is odd then there should be as many pixels to the left and to
% the right of the origin. mod(ODD -1,2) = 0 otherwise it is 1 setting this
% condition.
winTopCorner = winLowCorner - mod(WINDOW_SIZE-1,2);


% Get the down to the nearest integer position, and save any fractional shift
winCenter = floor(CENTER);
deltaWinCenter = CENTER - winCenter;


% Test for out of bounds conditions
LOW = winCenter - winLowCorner ;
TOP = VOLUME_SIZE - (winCenter + winTopCorner);
% Select only negative or zero shifts
% Since there is no pixel zero, shift low pad by 1
padLOW = (abs(LOW)+1) .* ( LOW <= 0 );
padTOP =  (abs(TOP) )    .* ( TOP <= 0 );

INDICES = [LOW + padLOW  ; winCenter + winTopCorner - padTOP ];
PADVALUES = [padLOW ; padTOP];
SHIFTS =   deltaWinCenter;

     % The window is often much larger than the particle. This is to ensure that all of the
     % delocalized information is available for full CTF restoration. Ideally, this is
     % fully represented in every particle, but as long as the particle itself is there,
     % don't worry if the high resolution information is not. 

     % Should keep track of this and figure into the quality weight somehow.

     availableArea = WINDOW_SIZE - PADVALUES(1,:) - PADVALUES(2,:);
  %   if any(availableArea - minSizeMask < -2) 
     if any(availableArea ./ WINDOW_SIZE < 0.5) 
        fprintf(['\nvs %d %d %d\nws %d %d %d\nmr %2.1f %2.1f %2.1f\n',...
                'minArea %d %d %d\navailArea %d %d %d\nc %2.1f %2.1f %2.1f\nwindowCutoff %f\n'], ...
                VOLUME_SIZE, WINDOW_SIZE, MASK_RADIUS, minSizeMask, availableArea, CENTER); 
        INDICES = 'noUse';
        PADVALUES = [availableArea];  
     end
     if any(isnan(PADVALUES(:)))
      fprintf('center %f %f %f\n',CENTER);
      fprintf('min %f %f %f\n',minSizeMask);
      fprintf('winLowCorner %f %f %f\n', winLowCorner);
      fprintf('top %f %f %f\n',winTopCorner);
      fprintf('%f %f %f\n',winCenter);
      fprintf('del %f %f %f\n',deltaWinCenter);

      fprintf('%f %f %f\n',TOP);
      error('\n\nFound a NaN in the pad values. But Why ben why?\n\n');
      INDICES='noUse';
      PADVALUES = [availableArea];
     end   
     % padSUM = sum(PADVALUES(:,1)) .* WINDOW_SIZE(2) .* WINDOW_SIZE(3);
     % padSUM = padSUM + sum(PADVALUES(:,2)) .* WINDOW_SIZE(1) .* WINDOW_SIZE(3);
     % padSUM = padSUM + sum(PADVALUES(:,3)) .* WINDOW_SIZE(1) .* WINDOW_SIZE(2);
     % padSUM = padSUM ./ prod(WINDOW_SIZE);
      
      % minimal sampling
     % if (padSUM) > .30 || any(any((PADVALUES - [WINDOW_SIZE;WINDOW_SIZE]) >0))
     %   INDICES = 'noUse';
     %   PADVALUES = padSUM;
     % end



end

