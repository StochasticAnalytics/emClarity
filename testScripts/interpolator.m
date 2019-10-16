classdef interpolator < handle
  %UNTITLED Summary of this class goes here
  %   Detailed explanation goes here
  
  properties
    texObject;
  end
  
  methods
    function [ obj, resampledVol ] = interpolator(inputVol, angles, shifts, convention, direction)
      %UNTITLED Construct an instance of this class
      %   Detailed explanation goes here
      

      [angles, shifts] = check_anglesAndShifts(obj,angles, shifts, convention, direction);
      boolDirection = check_inputs(obj,direction);
      
      [resampledVol, obj.texObject] = mexXform3d(inputVol, ...
                                                 angles, ...
                                                 shifts, ...
                                                 boolDirection);

    end
    
    function boolDirection = check_inputs(obj, direction)
      %METHOD1 Summary of this method goes here
      %   Detailed explanation goes here
      if strcmpi(direction, 'forward')
        boolDirection = true;
      elseif strcmpi(direction, 'inv')
        boolDirection = false;
      else
        error('Did not understand the transform direction %s\n', direction);
      end
      
    end

    
    function [ resampledVol ] = interp3d(obj, inputVol, angles, shifts, convention, direction)
      %METHOD1 Summary of this method goes here
      %   Detailed explanation goes here
      
      [angles, shifts] = check_anglesAndShifts(obj, angles, shifts, convention, direction);
      boolDirection = check_inputs(obj, direction);
      
      [resampledVol, ~] = mexXform3d(inputVol, ...
                                     angles, ...
                                     shifts, ...
                                     boolDirection, ...
                                     obj.texObject);
      
    end
    
   function [angles, shifts] = check_anglesAndShifts(obj,angles, shifts, convention, direction)
      %METHOD1 Summary of this method goes here
      %   Detailed explanation goes here
      if (numel(angles) == 3)
        angles = BH_defineMatrix(angles, convention, direction);
      end
      
      angles = gather(single(angles));
      
      if (numel(shifts) ~= 3)
        error('Shifts must have 3 elements');
      else
        shifts = gather(single(shifts));
      end
      
    end
    
  end
end

