function [ ROTATION_MATRIX ] = BH_defineMatrix( ANGLES, CONVENTION, DIRECTION)
%Create a rotation matrix.
%
%   Input variables:
%
%   ANGLES = [ first, second, third ] euler angles
%
%   CONVENTION = 'Protomo', 'Imod', 'Spider'
%     Protomo: ZXZ, passive, extrinsic
%     Imod   : ZYX, active,  intrinsic
%     Spider : ZYZ, active, extrinsic
%     Bah    : ZXZ, active,  extrinsic
%
%       For the time being, I am only updating Bah, and will revisit the other
%       conventions as time permits.
%
%   DIRECTION = forward : rotation from microscope frame to particle frame.
%               inverse : rotation from particle frame to microscope frame.
%   
%   Output variables:
%
%   ROTATION_MATRIX = 3d rotation matrix
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Goals & Restrictions
%
%   The individual rotation matrices defined below (cosine matrices) define an
%   active right handed transformation, that is a vector is rotated about the
%   given axis by the angle specified.
%
%   These are general, but in the scope of the BH_subTomo programs, they are
%   generally applied to an ndgrid which is transformed and used as the query to
%   an interpolation. 
%
%   Regardless of how they are used, the angles are interpreted to reflect an
%   active, intrinsic transformation on a particle, and the convention and 
%   direction are taken into account in order for this to work.
%
%	A good test is to create wedge masks of varying orientation because these
%	rely on this function for resampling. An aysmmetric wedge is more clear.
% 	e.g. [-50,70,0,90,0] [-50,70,60,90,0] [-50,70,0,90,60] [-50,70,90,90,-90]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   TODO
%     - update Protomo, Spider, Imod conversions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Convert from degrees to radians.
angles = ANGLES.*(pi./180);

%%%%%%%%%%%%%%%%%  Angle Definitions  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Creating the anonymous functions more than doubles the run time without.
% Using them also adds more time due to ~ 2x the trig calls.
% Leaving them here for a reference to the matrices used.
% Rx = @(t)[     1      0       0 ;...
%            0  cos(t)  -sin(t);...
%            0  sin(t)  cos(t) ];
% 
% Ry = @(t)[ cos(t)     0   sin(t);...
%            0      1        0;...
%           -sin(t)     0   cos(t) ];
% 
% Rz = @(t)[ cos(t)  -sin(t)      0;...
%           sin(t)  cos(t)      0;...
%             0       0      1 ];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
if strcmpi(DIRECTION, 'forward') || strcmpi(DIRECTION, 'invVector')
  angles =  -1.*angles;      
elseif strcmpi(DIRECTION, 'inv') || strcmpi(DIRECTION, 'forwardVector')
  % For interpolation the vectors are applied to a grid, so the sense must
    % be inverted to make the final transformation active.


    % In order to rotate the particle from a position defined by the input
    % angles, back to the proper reference frame, the sense is already
    % inverted, and just the order must be inverted.
    %
    % Think of this as taking an average in the proper frame, applying a given
    % rotation with 'forward', then this undoes that action.
    %
    % IMPORTANT NOTE: because the order is flipped, successive rotations by
    % multiple matrices must be right multplied for inverse operations. eg:
    % R1(e1,e2,e3) & R2(e4,e5,e6) then Rtot = R1 * R2 = e1*e2*e3*e4*e5*e6*Mat
    angles = flip(angles);
%   angles = [angles(3), angles(2), angles(1)];
      else
  error('Direction must be forward or inv, not %s', DIRECTION)
end

% Reduce number of trig functions
cosA = cos(angles);
sinA = sin(angles);
    
    
switch CONVENTION
  
  case 'Bah'
       
%     ROTATION_MATRIX = Rz(angles(3)) * Rx(angles(2)) * Rz(angles(1));
    ROTATION_MATRIX = [cosA(3),-sinA(3),0;...
                       sinA(3),cosA(3),0;...
                       0,0,1] * ...
                      [1,0,0; ...
                      0,cosA(2),-sinA(2);...
                      0,sinA(2),cosA(2)] * ...
                      [cosA(1),-sinA(1),0;...
                       sinA(1),cosA(1),0;...
                       0,0,1] ;
                     
  case 'TILT'
 
    ROTATION_MATRIX =  [cosA,0,sinA; ...
                             0,1,0;...
                         -sinA,0,cosA];

                      
  case 'SPIDER'
    
    
%     ROTATION_MATRIX = Rz(angles(3)) * Ry(angles(2)) * Rz(angles(1));

    ROTATION_MATRIX = [cosA(3),-sinA(3),0;...
                       sinA(3),cosA(3),0;...
                       0,0,1] * ...
                      [cosA(2),0,sinA(2); ...
                      0,1,0;...
                      -sinA(2),0,cosA(2)] * ...
                      [cosA(1),-sinA(1),0;...
                       sinA(1),cosA(1),0;...
                       0,0,1] ;

  case 'Helical'

    
%   ROTATION_MATRIX = Ry(angles(3)) * Rx(angles(2)) * Ry(angles(1));  

    ROTATION_MATRIX = [cosA(3),0,sinA(3); ...
                      0,1,0;...
                      -sinA(3),0,cosA(3)] * ...
                      [1,0,0; ...
                      0,cosA(2),-sinA(2);...
                      0,sinA(2),cosA(2)] * ...
                      [cosA(1),0,sinA(1); ...
                      0,1,0;...
                      -sinA(1),0,cosA(1)];

  case 'IMOD'

    cosA = cos(angles);
    sinA = sin(angles);
    
%   ROTATION_MATRIX = Rz(angles(3)) * Ry(angles(2)) * Rx(angles(1));  

    ROTATION_MATRIX = [cosA(3),-sinA(3),0;...
                       sinA(3),cosA(3),0;...
                       0,0,1] * ...
                      [cosA(2),0,sinA(2); ...
                      0,1,0;...
                      -sinA(2),0,cosA(2)] * ...
                      [cosA(1),-sinA(1),0;...
                       sinA(1),cosA(1),0;...
                       0,0,1] ;
  
  otherwise
    error('Convention must be Bah,SPI,Helical, not %s', CONVENTION)
end

end % end of function defineMatrix.

%     case 'Protomo'
%         % passive, intrinsic, Z X Z
%         % i3euler e1 e2 e3
% 
%         if strcmpi(dir,'forward')
% 
%         elseif strcmpi(dir, 'inv')
%             ang = -1.* [ang(3), ang(2), ang(1)];
% 
% 
%         elseif strcmpi(dir, 'i3')
%             ang = [ang(3), ang(2), ang(1)];

%         end
% 
%         RotMat = Rz(ang(3)) * Rx(ang(2)) * Rz(ang(1));
% 
% 
%     case 'Imod'
%         % active, extrinsic, Z Y X
% 
%         if strcmpi(dir, 'forward')
%             ang = -1 .* ang ;
%         end
% 
%         RotMat = Rx(ang(3))*Ry(ang(2))*Rz(ang(1)) ;
% 
%     case 'Spider'
%         % passive, extrinsic Z Y Z
%         % Note that in their documents they refer to the "object"
%         % rotating clockwise, which sounds active, but this is the
%         % same as the CS anti-clockwise, which is just a passive 
%         % (alias) rotation. I believe Frealign, and Relion also use.
% 
%         % Spider puts the origin at top left, with first z on top
% 
%         RotMat = Rz(-e3)*Ry(-e2)*Rz(-e1) ;
% 
%     case '2d'
%         % active rotation, second two euler angles are dummy var
% 
%        RotMat = Rz(e1);
%        RotMat = RotMat(1:2,1:2) ;
% 
%    case 'NegProtomo'
%         % passive, intrinsic, Z X Z
%         % i3euler e1 e2 e3
% 
% 
%        RotMat = Rz(-ang(1)) * Rx(-ang(2)) * Rz(-ang(3)) ;
% 
% end
                
            

