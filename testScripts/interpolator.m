classdef interpolator < handle
  %UNTITLED Summary of this class goes here
  %   Detailed explanation goes here
  
  properties
    texObject;
    nSymMats;
    symmetry_matrices;
  end
  
  methods
    function [ obj, resampledVol ] = interpolator(inputVol, angles, shifts, convention, direction, symmetry)
      %UNTITLED Construct an instance of this class
      %   Detailed explanation goes here
      
      check_symmetry(obj, symmetry);
      [angles, shifts] = check_anglesAndShifts(obj,angles, shifts, convention, direction);
      boolDirection = check_inputs(obj,direction);
      
      [resampledVol, obj.texObject] = mexXform3d(inputVol, ...
                                                 angles, ...
                                                 shifts, ...
                                                 boolDirection);
      for iSym = 2:obj.nSymMats
       
        symAngles = gather(obj.symmetry_matrices{iSym} * angles);
        resampledVol = resampledVol + mexXform3d(inputVol, ...
                                                 symAngles, ...
                                                 shifts, ...
                                                 boolDirection, ...
                                                 obj.texObject);
        
      end
      
      if (obj.nSymMats > 1)
        resampledVol = resampledVol ./ obj.nSymMats;
      end

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

    
    function [ resampledVol ] = interp3d(obj, inputVol, angles, shifts, convention, direction, symmetry)
      %METHOD1 Summary of this method goes here
      %   Detailed explanation goes here
      
      check_symmetry(obj, symmetry);
      [angles, shifts] = check_anglesAndShifts(obj, angles, shifts, convention, direction);
      boolDirection = check_inputs(obj, direction);
 
      resampledVol = zeros(size(inputVol), 'single', 'gpuArray');
      for iSym = 1:obj.nSymMats
       
        symAngles = gather(obj.symmetry_matrices{iSym} * angles);
        resampledVol = resampledVol + mexXform3d(inputVol, ...
                                                 symAngles, ...
                                                 shifts, ...
                                                 boolDirection, ...
                                                 obj.texObject);
        
      end
      
      if (obj.nSymMats > 1)
        resampledVol = resampledVol ./ obj.nSymMats;
      end

      
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
    
   function [  ] = check_symmetry(obj, symmetry)
     
     if (strcmpi(symmetry, 'C1'))
       obj.symmetry_matrices = cell(1);
       symmetry_type = 0;
     elseif (strcmpi(symmetry, 'C2'))
       obj.symmetry_matrices = cell(2);
       symmetry_type = 0;
     elseif (strcmpi(symmetry, 'C3'))
        obj.symmetry_matrices = cell(3);
        symmetry_type = 0;
     elseif (strcmpi(symmetry, 'C4'))
        obj.symmetry_matrices = cell(4);
        symmetry_type = 0;
      elseif (strcmpi(symmetry, 'C5'))
        obj.symmetry_matrices = cell(5);
        symmetry_type = 0;
      elseif (strcmpi(symmetry, 'C6'))
        obj.symmetry_matrices = cell(6);
        symmetry_type = 0;
     else
       error('Only C1-6 have been implemented in the new interp class');
     end
     
     switch symmetry_type
       case 0
         obj.nSymMats = length(obj.symmetry_matrices);
         symInc = 360 / obj.nSymMats;
         for iSym = 0: obj.nSymMats - 1
          obj.symmetry_matrices{iSym + 1} = BH_defineMatrix([0,0,iSym.*symInc], 'Bah', 'forward');
         end
     end
      
     
   end
    
  end
end

