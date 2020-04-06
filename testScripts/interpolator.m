classdef interpolator < handle
  %UNTITLED Summary of this class goes here
  %   Detailed explanation goes here
  
  properties
    % Pointers to underlying cuda objects
    texObject = '';
    cuArray = '';
    % Symmetry matrices
    nSymMats;
    symmetry_matrices;
    symmetry_type = '';
    % The resampled volume
%     resampledVol;
  end
  
  methods
    function [ obj, resampledVol ] = interpolator(inputVol, angles, shifts, convention, direction, symmetry, varargin)
      %UNTITLED Construct an instance of this class
      %   Detailed explanation goes here
      
      if (nargin == 7)
        useOnlyOnce = varargin{1};
      else
        useOnlyOnce = false;
      end
      
      check_symmetry(obj, symmetry);
      [angles, shifts] = check_anglesAndShifts(obj,angles, shifts, convention, direction);
      boolDirection = check_inputs(obj,direction);
      
      [resampledVol, obj.texObject, obj.cuArray] = mexXform3d(inputVol, ...
                                                                 angles, ...
                                                                 shifts, ...
                                                                 boolDirection);
      for iSym = 2:obj.nSymMats
       
        symAngles = gather(angles*obj.symmetry_matrices{iSym});
        resampledVol = resampledVol + mexXform3d(inputVol, ...
                                                 symAngles, ...
                                                 shifts, ...
                                                 boolDirection, ...
                                                 obj.texObject);
        
      end
      
      if (obj.nSymMats > 1)
        resampledVol = resampledVol ./ obj.nSymMats;
      end
      
      if (useOnlyOnce)
        obj.deallocate();
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
 
%       resampledVol = zeros(size(inputVol), 'single', 'gpuArray');
      for iSym = 1:obj.nSymMats
       
        symAngles = gather(angles*obj.symmetry_matrices{iSym});
        if (iSym == 1)
        resampledVol =   mexXform3d(inputVol, ...
                                                 symAngles, ...
                                                 shifts, ...
                                                 boolDirection, ...
                                                 obj.texObject);  
        else
        resampledVol = resampledVol + mexXform3d(inputVol, ...
                                                 symAngles, ...
                                                 shifts, ...
                                                 boolDirection, ...
                                                 obj.texObject);
        end
        
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
     
     if isempty(obj.symmetry_type |  ~strcmpi(symmetry, obj.symmetry_type))
       
       if (strcmpi(symmetry, obj.symmetry_type))
         return;
       else
         obj.symmetry_type = symmetry;
%          warning('The requested symmetry (%s) is different from that initialized (%s)\n',symmetry,obj.symmetry_type);
       end
       
       if (strcmpi(symmetry, 'C1'))
         obj.symmetry_matrices = cell(1);
       elseif (strcmpi(symmetry, 'C2'))
         obj.symmetry_matrices = cell(2);
       elseif (strcmpi(symmetry, 'C3'))
          obj.symmetry_matrices = cell(3);
       elseif (strcmpi(symmetry, 'C4'))
          obj.symmetry_matrices = cell(4);
        elseif (strcmpi(symmetry, 'C5'))
          obj.symmetry_matrices = cell(5);
        elseif (strcmpi(symmetry, 'C6'))
          obj.symmetry_matrices = cell(6);
       else
         error('Only C1-6 have been implemented in the new interp class');
       end
     
%      switch symmetry_type
%        case 0
         obj.nSymMats = length(obj.symmetry_matrices);
         symInc = 360 / obj.nSymMats;
         for iSym = 0: obj.nSymMats - 1
          obj.symmetry_matrices{iSym + 1} = BH_defineMatrix([0,0,iSym.*symInc], 'Bah', 'forward');
         end
%      end
     end

      
     
   end
   
   function [  ] = delete(obj)
     if ~isempty(obj.texObject) && ~isempty(obj.cuArray)
      mexXform3d(obj.texObject, obj.cuArray);
     end
     obj.symmetry_matrices = [];
     obj.nSymMats = [];
     obj.cuArray = '';
     obj.texObject = '';
   end
   
  end
end

