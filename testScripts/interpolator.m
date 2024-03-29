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
    % 
    input_size = '';
    dummy_vol = '';
    
    % Keeping the texture alive works much of the time, but for some reason
    % deallocating the cuArray can lead to seg faults only some of the
    % time. It seems to be worse on 2080 ti vs titans or V100. For now,
    % keep it alive for symmetry ops, but delete after interp, and store
    % the input volume here.
    input_volume = '';
    make_tex_persistent = false;
  end
  
  methods
    function [ obj, resampledVol ] = interpolator(inputVol, angles, shifts, convention, direction, symmetry, varargin)
      %UNTITLED Construct an instance of this classlt
      
      %   Detailed explanation goes here
      
      if (nargin == 7)
        useOnlyOnce = varargin{1};
      else
        useOnlyOnce = false;
      end
      
      check_symmetry(obj, symmetry, convention);
      [angles, shifts] = check_anglesAndShifts(obj,angles, shifts, convention, direction);
      boolDirection = check_inputs(obj,direction);
      
      % Set the input size
      obj.input_size = uint64((size(inputVol)));
      obj.dummy_vol = inputVol(1:2,1:2,1:2);
      
      [resampledVol, obj.texObject, obj.cuArray] = mexXform3d(obj.input_size,...
                                                              inputVol, ...
                                                                 angles, ...
                                                                 shifts, ...
                                                                 boolDirection);
      for iSym = 2:obj.nSymMats
       
        symAngles = gather(angles*obj.symmetry_matrices{iSym});
        resampledVol = resampledVol + mexXform3d(obj.input_size, ...
                                                 obj.dummy_vol,...
                                                 symAngles, ...
                                                 shifts, ...
                                                 boolDirection, ...
                                                 obj.texObject);
        
      end
      
      if (obj.nSymMats > 1)
        resampledVol = resampledVol ./ obj.nSymMats;
      end
      
      if (useOnlyOnce)
        obj.delete();
      elseif (obj.make_tex_persistent)
        % The texture is kept alive, so no need to also store the volume
        % for re-initialization of the texture
        
        % Until I can figure out the segfaults, I'm not even putting in
        % somethign to enable this.
      else
        obj.input_volume = inputVol;
        obj.free_cuda_objects();
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

    
    function [ resampledVol ] = interp3d(obj, angles, shifts, convention, direction, symmetry)
      %METHOD1 Summary of this method goes here
      %   Detailed explanation goes here
      
      check_symmetry(obj, symmetry, convention);
      [angles, shifts] = check_anglesAndShifts(obj, angles, shifts, convention, direction);
      boolDirection = check_inputs(obj, direction);
 
%       resampledVol = zeros(size(inputVol), 'single', 'gpuArray');
      for iSym = 1:obj.nSymMats
       
        symAngles = gather(angles*obj.symmetry_matrices{iSym});
        if (iSym == 1)
  
          if (obj.make_tex_persistent)
           resampledVol =  mexXform3d(obj.input_size, ...
                                       obj.dummy_vol,...
                                       symAngles, ...
                                       shifts, ...
                                       boolDirection, ...
                                       obj.texObject);
          else
            % If we didn't make the tex obj persistent then we need to
            % return it here, for use at least with the symmetric vols
            [resampledVol, obj.texObject, obj.cuArray] =  ...
                              mexXform3d(obj.input_size, ...
                                         obj.input_volume,...
                                         symAngles, ...
                                         shifts, ...
                                         boolDirection);   

          end
 
        else
          
           % w or w/o persistent tex, this call is the same.
           resampledVol = resampledVol + mexXform3d(obj.input_size, ...
                                                     obj.dummy_vol,...
                                                     symAngles, ...
                                                     shifts, ...
                                                     boolDirection, ...
                                                     obj.texObject);
        end
        
      end
      
      if (obj.nSymMats > 1)
        resampledVol = resampledVol ./ obj.nSymMats;
      end
      
      if ~(obj.make_tex_persistent)
        obj.free_cuda_objects();
      end

      
    end
    
   function [angles, shifts] = check_anglesAndShifts(obj,angles, shifts, convention, direction)
      %METHOD1 Summary of this method goes here
      %   Detailed explanation goes here
      if (numel(angles) == 3)
        angles = BH_defineMatrix(angles, convention, direction);
      end
      
      angles = gather(single(reshape(angles,3,3)));
      
      if (numel(shifts) ~= 3)
        error('Shifts must have 3 elements');
      else
        shifts = gather(single(shifts));
      end
      
   end
    
   function [  ] = check_symmetry(obj, symmetry, convention)
     
     if isempty(obj.symmetry_type |  ~strcmpi(symmetry, obj.symmetry_type))
       
       if (strcmpi(symmetry, obj.symmetry_type))
         return;
       else
         obj.symmetry_type = symmetry;
%          warning('The requested symmetry (%s) is different from that initialized (%s)\n',symmetry,obj.symmetry_type);
       end
       
       if (symmetry(1) ~= 'C' && ~strcmp(convention, 'Bah'))
         error('Alternate conventions like Helical (%s encountered) only support CX symmetry\n', convention)
       end
       
       switch symmetry(1)
         case 'C'
           if (length(symmetry) < 2)
             error('Cyclic symmetry requires an int specifying CX');
           end
            obj.symmetry_matrices = cell(EMC_str2double(symmetry(2:end)),1); 
            obj.nSymMats = length(obj.symmetry_matrices);
            symInc = 360 / obj.nSymMats;
            for iSym = 0: obj.nSymMats - 1
              obj.symmetry_matrices{iSym + 1} = BH_defineMatrix([0,0,iSym.*symInc], convention, 'forward');
            end  
         case 'D'
           if (length(symmetry) < 2)
             error('D symmetry requires an int specifying DX');
           end
            n_inplane = EMC_str2double(symmetry(2:end));
            obj.symmetry_matrices = cell(n_inplane,1); 
            obj.nSymMats = length(obj.symmetry_matrices);
            symInc = 360 / obj.nSymMats;
            for iSym = 0: obj.nSymMats - 1
              obj.symmetry_matrices{iSym + 1} = BH_defineMatrix([0,0,iSym.*symInc], 'Bah', 'forward');
              obj.symmetry_matrices{iSym + 1 + n_inplane} =  obj.symmetry_matrices{iSym + 1} .* [1,-1,1;1,-1,1;1,1,-1];
            end             
         case 'O'
           if (length(symmetry) > 1)
             error('Octahedral symmetry requires no int');
           end
             obj.nSymMats = 24;
             obj.symmetry_matrices = cell(obj.nSymMats,1); 
            % Pulling this from symmetry matrix.cpp from cisTEM. If it
            % works, get rid of the addition and just change it to 1 based
            % indexing manually.
            obj.symmetry_matrices{ 0 + 1} = [ 1.000000, 0.000000, 0.000000; 0.000000, 1.000000, 0.000000; 0.000000, 0.000000, 1.000000]';
            obj.symmetry_matrices{ 1 + 1} = [ 0.000000, 0.000000,-1.000000;-1.000000, 0.000000, 0.000000; 0.000000, 1.000000, 0.000000]';
            obj.symmetry_matrices{ 2 + 1} = [ 0.000000, 0.000000,-1.000000; 0.000000,-1.000000, 0.000000;-1.000000, 0.000000, 0.000000]';
            obj.symmetry_matrices{ 3 + 1} = [ 0.000000, 0.000000,-1.000000; 1.000000, 0.000000, 0.000000; 0.000000,-1.000000, 0.000000]';
            obj.symmetry_matrices{ 4 + 1} = [ 0.000000, 0.000000,-1.000000; 0.000000, 1.000000, 0.000000; 1.000000, 0.000000, 0.000000]';
            obj.symmetry_matrices{ 5 + 1} = [ 0.000000, 0.000000, 1.000000; 0.000000, 1.000000, 0.000000;-1.000000, 0.000000, 0.000000]';
            obj.symmetry_matrices{ 6 + 1} = [ 0.000000, 0.000000, 1.000000;-1.000000, 0.000000, 0.000000; 0.000000,-1.000000, 0.000000]';
            obj.symmetry_matrices{ 7 + 1} = [ 0.000000, 0.000000, 1.000000; 0.000000,-1.000000, 0.000000; 1.000000, 0.000000, 0.000000]';
            obj.symmetry_matrices{ 8 + 1} = [ 0.000000, 0.000000, 1.000000; 1.000000, 0.000000, 0.000000; 0.000000, 1.000000, 0.000000]';
            obj.symmetry_matrices{ 9 + 1} = [ 0.000000, 1.000000, 0.000000; 0.000000, 0.000000,-1.000000;-1.000000, 0.000000, 0.000000]';
            obj.symmetry_matrices{10 + 1} = [-1.000000, 0.000000, 0.000000; 0.000000, 0.000000,-1.000000; 0.000000,-1.000000, 0.000000]';
            obj.symmetry_matrices{11 + 1} = [ 0.000000,-1.000000, 0.000000; 0.000000, 0.000000,-1.000000; 1.000000, 0.000000, 0.000000]';
            obj.symmetry_matrices{12 + 1} = [ 1.000000, 0.000000, 0.000000; 0.000000, 0.000000,-1.000000; 0.000000, 1.000000, 0.000000]';
            obj.symmetry_matrices{13 + 1} = [-1.000000, 0.000000, 0.000000; 0.000000, 0.000000, 1.000000; 0.000000, 1.000000, 0.000000]';
            obj.symmetry_matrices{14 + 1} = [ 0.000000,-1.000000, 0.000000; 0.000000, 0.000000, 1.000000;-1.000000, 0.000000, 0.000000]';
            obj.symmetry_matrices{15 + 1} = [ 1.000000, 0.000000, 0.000000; 0.000000, 0.000000, 1.000000; 0.000000,-1.000000, 0.000000]';
            obj.symmetry_matrices{16 + 1} = [ 0.000000, 1.000000, 0.000000; 0.000000, 0.000000, 1.000000; 1.000000, 0.000000, 0.000000]';
            obj.symmetry_matrices{17 + 1} = [-1.000000, 0.000000, 0.000000; 0.000000, 1.000000, 0.000000; 0.000000, 0.000000,-1.000000]';
            obj.symmetry_matrices{18 + 1} = [ 0.000000,-1.000000, 0.000000;-1.000000, 0.000000, 0.000000; 0.000000, 0.000000,-1.000000]';
            obj.symmetry_matrices{19 + 1} = [ 1.000000, 0.000000, 0.000000; 0.000000,-1.000000, 0.000000; 0.000000, 0.000000,-1.000000]';
            obj.symmetry_matrices{20 + 1} = [ 0.000000, 1.000000, 0.000000; 1.000000, 0.000000, 0.000000; 0.000000, 0.000000,-1.000000]';
            obj.symmetry_matrices{21 + 1} = [ 0.000000, 1.000000, 0.000000;-1.000000, 0.000000, 0.000000; 0.000000, 0.000000, 1.000000]';
            obj.symmetry_matrices{22 + 1} = [-1.000000, 0.000000, 0.000000; 0.000000,-1.000000, 0.000000; 0.000000, 0.000000, 1.000000]';
            obj.symmetry_matrices{23 + 1} = [ 0.000000,-1.000000, 0.000000; 1.000000, 0.000000, 0.000000; 0.000000, 0.000000, 1.000000]';
           
         case 'I'
           
          obj.nSymMats = 60;
          obj.symmetry_matrices = cell(obj.nSymMats,1); 
          if (length(symmetry) < 2)
        
            obj.symmetry_matrices{ 0 + 1} = [ 1.000000, 0.000000, 0.000000; 0.000000, 1.000000, 0.000000; 0.000000, 0.000000, 1.000000]';
            obj.symmetry_matrices{ 1 + 1} = [-0.500000,-0.309017,-0.809017;-0.309017,-0.809017, 0.500000;-0.809017, 0.500000, 0.309017]';
            obj.symmetry_matrices{ 2 + 1} = [-0.309017,-0.809017,-0.500000; 0.809017,-0.500000, 0.309017;-0.500000,-0.309017, 0.809017]';
            obj.symmetry_matrices{ 3 + 1} = [ 0.309017,-0.809017,-0.500000; 0.809017, 0.500000,-0.309017; 0.500000,-0.309017, 0.809017]';
            obj.symmetry_matrices{ 4 + 1} = [ 0.500000,-0.309017,-0.809017;-0.309017, 0.809017,-0.500000; 0.809017, 0.500000, 0.309017]';
            obj.symmetry_matrices{ 5 + 1} = [ 0.000000, 0.000000,-1.000000;-1.000000, 0.000000, 0.000000; 0.000000, 1.000000, 0.000000]';
            obj.symmetry_matrices{ 6 + 1} = [ 0.809017,-0.500000,-0.309017;-0.500000,-0.309017,-0.809017; 0.309017, 0.809017,-0.500000]';
            obj.symmetry_matrices{ 7 + 1} = [ 0.500000, 0.309017,-0.809017;-0.309017,-0.809017,-0.500000;-0.809017, 0.500000,-0.309017]';
            obj.symmetry_matrices{ 8 + 1} = [-0.500000, 0.309017,-0.809017; 0.309017,-0.809017,-0.500000;-0.809017,-0.500000, 0.309017]';
            obj.symmetry_matrices{ 9 + 1} = [-0.809017,-0.500000,-0.309017; 0.500000,-0.309017,-0.809017; 0.309017,-0.809017, 0.500000]';
            obj.symmetry_matrices{10 + 1} = [ 0.000000,-1.000000, 0.000000; 0.000000, 0.000000,-1.000000; 1.000000, 0.000000, 0.000000]';
            obj.symmetry_matrices{11 + 1} = [-0.309017,-0.809017, 0.500000; 0.809017,-0.500000,-0.309017; 0.500000, 0.309017, 0.809017]';
            obj.symmetry_matrices{12 + 1} = [ 0.809017,-0.500000, 0.309017; 0.500000, 0.309017,-0.809017; 0.309017, 0.809017, 0.500000]';
            obj.symmetry_matrices{13 + 1} = [ 0.809017, 0.500000,-0.309017;-0.500000, 0.309017,-0.809017;-0.309017, 0.809017, 0.500000]';
            obj.symmetry_matrices{14 + 1} = [-0.309017, 0.809017,-0.500000;-0.809017,-0.500000,-0.309017;-0.500000, 0.309017, 0.809017]';
            obj.symmetry_matrices{15 + 1} = [-1.000000, 0.000000, 0.000000; 0.000000,-1.000000, 0.000000; 0.000000, 0.000000, 1.000000]';
            obj.symmetry_matrices{16 + 1} = [-0.500000,-0.309017,-0.809017; 0.309017, 0.809017,-0.500000; 0.809017,-0.500000,-0.309017]';
            obj.symmetry_matrices{17 + 1} = [-0.309017,-0.809017,-0.500000;-0.809017, 0.500000,-0.309017; 0.500000, 0.309017,-0.809017]';
            obj.symmetry_matrices{18 + 1} = [ 0.309017,-0.809017,-0.500000;-0.809017,-0.500000, 0.309017;-0.500000, 0.309017,-0.809017]';
            obj.symmetry_matrices{19 + 1} = [ 0.500000,-0.309017,-0.809017; 0.309017,-0.809017, 0.500000;-0.809017,-0.500000,-0.309017]';
            obj.symmetry_matrices{20 + 1} = [ 0.000000, 0.000000,-1.000000; 1.000000, 0.000000, 0.000000; 0.000000,-1.000000, 0.000000]';
            obj.symmetry_matrices{21 + 1} = [ 0.809017,-0.500000,-0.309017; 0.500000, 0.309017, 0.809017;-0.309017,-0.809017, 0.500000]';
            obj.symmetry_matrices{22 + 1} = [ 0.500000, 0.309017,-0.809017; 0.309017, 0.809017, 0.500000; 0.809017,-0.500000, 0.309017]';
            obj.symmetry_matrices{23 + 1} = [-0.500000, 0.309017,-0.809017;-0.309017, 0.809017, 0.500000; 0.809017, 0.500000,-0.309017]';
            obj.symmetry_matrices{24 + 1} = [-0.809017,-0.500000,-0.309017;-0.500000, 0.309017, 0.809017;-0.309017, 0.809017,-0.500000]';
            obj.symmetry_matrices{25 + 1} = [ 0.000000,-1.000000, 0.000000; 0.000000, 0.000000, 1.000000;-1.000000, 0.000000, 0.000000]';
            obj.symmetry_matrices{26 + 1} = [-0.309017,-0.809017, 0.500000;-0.809017, 0.500000, 0.309017;-0.500000,-0.309017,-0.809017]';
            obj.symmetry_matrices{27 + 1} = [ 0.809017,-0.500000, 0.309017;-0.500000,-0.309017, 0.809017;-0.309017,-0.809017,-0.500000]';
            obj.symmetry_matrices{28 + 1} = [ 0.809017, 0.500000,-0.309017; 0.500000,-0.309017, 0.809017; 0.309017,-0.809017,-0.500000]';
            obj.symmetry_matrices{29 + 1} = [-0.309017, 0.809017,-0.500000; 0.809017, 0.500000, 0.309017; 0.500000,-0.309017,-0.809017]';
            obj.symmetry_matrices{30 + 1} = [-1.000000, 0.000000, 0.000000; 0.000000, 1.000000, 0.000000; 0.000000, 0.000000,-1.000000]';
            obj.symmetry_matrices{31 + 1} = [ 0.500000, 0.309017, 0.809017;-0.309017,-0.809017, 0.500000; 0.809017,-0.500000,-0.309017]';
            obj.symmetry_matrices{32 + 1} = [ 0.309017, 0.809017, 0.500000; 0.809017,-0.500000, 0.309017; 0.500000, 0.309017,-0.809017]';
            obj.symmetry_matrices{33 + 1} = [-0.309017, 0.809017, 0.500000; 0.809017, 0.500000,-0.309017;-0.500000, 0.309017,-0.809017]';
            obj.symmetry_matrices{34 + 1} = [-0.500000, 0.309017, 0.809017;-0.309017, 0.809017,-0.500000;-0.809017,-0.500000,-0.309017]';
            obj.symmetry_matrices{35 + 1} = [ 0.000000, 0.000000, 1.000000;-1.000000, 0.000000, 0.000000; 0.000000,-1.000000, 0.000000]';
            obj.symmetry_matrices{36 + 1} = [-0.809017, 0.500000, 0.309017;-0.500000,-0.309017,-0.809017;-0.309017,-0.809017, 0.500000]';
            obj.symmetry_matrices{37 + 1} = [-0.500000,-0.309017, 0.809017;-0.309017,-0.809017,-0.500000; 0.809017,-0.500000, 0.309017]';
            obj.symmetry_matrices{38 + 1} = [ 0.500000,-0.309017, 0.809017; 0.309017,-0.809017,-0.500000; 0.809017, 0.500000,-0.309017]';
            obj.symmetry_matrices{39 + 1} = [ 0.809017, 0.500000, 0.309017; 0.500000,-0.309017,-0.809017;-0.309017, 0.809017,-0.500000]';
            obj.symmetry_matrices{40 + 1} = [ 0.000000, 1.000000, 0.000000; 0.000000, 0.000000,-1.000000;-1.000000, 0.000000, 0.000000]';
            obj.symmetry_matrices{41 + 1} = [ 0.309017, 0.809017,-0.500000; 0.809017,-0.500000,-0.309017;-0.500000,-0.309017,-0.809017]';
            obj.symmetry_matrices{42 + 1} = [-0.809017, 0.500000,-0.309017; 0.500000, 0.309017,-0.809017;-0.309017,-0.809017,-0.500000]';
            obj.symmetry_matrices{43 + 1} = [-0.809017,-0.500000, 0.309017;-0.500000, 0.309017,-0.809017; 0.309017,-0.809017,-0.500000]';
            obj.symmetry_matrices{44 + 1} = [ 0.309017,-0.809017, 0.500000;-0.809017,-0.500000,-0.309017; 0.500000,-0.309017,-0.809017]';
            obj.symmetry_matrices{45 + 1} = [ 1.000000, 0.000000, 0.000000; 0.000000,-1.000000, 0.000000; 0.000000, 0.000000,-1.000000]';
            obj.symmetry_matrices{46 + 1} = [ 0.500000, 0.309017, 0.809017; 0.309017, 0.809017,-0.500000;-0.809017, 0.500000, 0.309017]';
            obj.symmetry_matrices{47 + 1} = [ 0.309017, 0.809017, 0.500000;-0.809017, 0.500000,-0.309017;-0.500000,-0.309017, 0.809017]';
            obj.symmetry_matrices{48 + 1} = [-0.309017, 0.809017, 0.500000;-0.809017,-0.500000, 0.309017; 0.500000,-0.309017, 0.809017]';
            obj.symmetry_matrices{49 + 1} = [-0.500000, 0.309017, 0.809017; 0.309017,-0.809017, 0.500000; 0.809017, 0.500000, 0.309017]';
            obj.symmetry_matrices{50 + 1} = [ 0.000000, 0.000000, 1.000000; 1.000000, 0.000000, 0.000000; 0.000000, 1.000000, 0.000000]';
            obj.symmetry_matrices{51 + 1} = [-0.809017, 0.500000, 0.309017; 0.500000, 0.309017, 0.809017; 0.309017, 0.809017,-0.500000]';
            obj.symmetry_matrices{52 + 1} = [-0.500000,-0.309017, 0.809017; 0.309017, 0.809017, 0.500000;-0.809017, 0.500000,-0.309017]';
            obj.symmetry_matrices{53 + 1} = [ 0.500000,-0.309017, 0.809017;-0.309017, 0.809017, 0.500000;-0.809017,-0.500000, 0.309017]';
            obj.symmetry_matrices{54 + 1} = [ 0.809017, 0.500000, 0.309017;-0.500000, 0.309017, 0.809017; 0.309017,-0.809017, 0.500000]';
            obj.symmetry_matrices{55 + 1} = [ 0.000000, 1.000000, 0.000000; 0.000000, 0.000000, 1.000000; 1.000000, 0.000000, 0.000000]';
            obj.symmetry_matrices{56 + 1} = [ 0.309017, 0.809017,-0.500000;-0.809017, 0.500000, 0.309017; 0.500000, 0.309017, 0.809017]';
            obj.symmetry_matrices{57 + 1} = [-0.809017, 0.500000,-0.309017;-0.500000,-0.309017, 0.809017; 0.309017, 0.809017, 0.500000]';
            obj.symmetry_matrices{58 + 1} = [-0.809017,-0.500000, 0.309017; 0.500000,-0.309017, 0.809017;-0.309017, 0.809017, 0.500000]';
            obj.symmetry_matrices{59 + 1} = [ 0.309017,-0.809017, 0.500000; 0.809017, 0.500000, 0.309017;-0.500000, 0.309017, 0.809017]';
          elseif strcmp(symmetry(2),'2')
            obj.symmetry_matrices{ 0 + 1} = [ 1.000000, 0.000000, 0.000000; 0.000000, 1.000000, 0.000000; 0.000000, 0.000000, 1.000000]';
            obj.symmetry_matrices{ 1 + 1} = [ 0.500000,-0.809017,-0.309017;-0.809017,-0.309017,-0.500000; 0.309017, 0.500000,-0.809017]';
            obj.symmetry_matrices{ 2 + 1} = [ 0.309017,-0.500000, 0.809017;-0.500000,-0.809017,-0.309017; 0.809017,-0.309017,-0.500000]';
            obj.symmetry_matrices{ 3 + 1} = [-0.309017, 0.500000, 0.809017;-0.500000,-0.809017, 0.309017; 0.809017,-0.309017, 0.500000]';
            obj.symmetry_matrices{ 4 + 1} = [-0.500000, 0.809017,-0.309017;-0.809017,-0.309017, 0.500000; 0.309017, 0.500000, 0.809017]';
            obj.symmetry_matrices{ 5 + 1} = [ 0.000000, 0.000000,-1.000000;-1.000000, 0.000000, 0.000000; 0.000000, 1.000000, 0.000000]';
            obj.symmetry_matrices{ 6 + 1} = [-0.309017,-0.500000, 0.809017; 0.500000,-0.809017,-0.309017; 0.809017, 0.309017, 0.500000]';
            obj.symmetry_matrices{ 7 + 1} = [-0.809017, 0.309017, 0.500000; 0.309017,-0.500000, 0.809017; 0.500000, 0.809017, 0.309017]';
            obj.symmetry_matrices{ 8 + 1} = [-0.809017, 0.309017,-0.500000;-0.309017, 0.500000, 0.809017; 0.500000, 0.809017,-0.309017]';
            obj.symmetry_matrices{ 9 + 1} = [-0.309017,-0.500000,-0.809017;-0.500000, 0.809017,-0.309017; 0.809017, 0.309017,-0.500000]';
            obj.symmetry_matrices{10 + 1} = [ 0.000000,-1.000000, 0.000000; 0.000000, 0.000000,-1.000000; 1.000000, 0.000000, 0.000000]';
            obj.symmetry_matrices{11 + 1} = [-0.809017,-0.309017,-0.500000;-0.309017,-0.500000, 0.809017;-0.500000, 0.809017, 0.309017]';
            obj.symmetry_matrices{12 + 1} = [-0.500000,-0.809017,-0.309017;-0.809017, 0.309017, 0.500000;-0.309017, 0.500000,-0.809017]';
            obj.symmetry_matrices{13 + 1} = [-0.500000,-0.809017, 0.309017;-0.809017, 0.309017,-0.500000; 0.309017,-0.500000,-0.809017]';
            obj.symmetry_matrices{14 + 1} = [-0.809017,-0.309017, 0.500000;-0.309017,-0.500000,-0.809017; 0.500000,-0.809017, 0.309017]';
            obj.symmetry_matrices{15 + 1} = [-1.000000, 0.000000, 0.000000; 0.000000,-1.000000, 0.000000; 0.000000, 0.000000, 1.000000]';
            obj.symmetry_matrices{16 + 1} = [ 0.500000,-0.809017,-0.309017; 0.809017, 0.309017, 0.500000;-0.309017,-0.500000, 0.809017]';
            obj.symmetry_matrices{17 + 1} = [ 0.309017,-0.500000, 0.809017; 0.500000, 0.809017, 0.309017;-0.809017, 0.309017, 0.500000]';
            obj.symmetry_matrices{18 + 1} = [-0.309017, 0.500000, 0.809017; 0.500000, 0.809017,-0.309017;-0.809017, 0.309017,-0.500000]';
            obj.symmetry_matrices{19 + 1} = [-0.500000, 0.809017,-0.309017; 0.809017, 0.309017,-0.500000;-0.309017,-0.500000,-0.809017]';
            obj.symmetry_matrices{20 + 1} = [ 0.000000, 0.000000,-1.000000; 1.000000, 0.000000, 0.000000; 0.000000,-1.000000, 0.000000]';
            obj.symmetry_matrices{21 + 1} = [-0.309017,-0.500000, 0.809017;-0.500000, 0.809017, 0.309017;-0.809017,-0.309017,-0.500000]';
            obj.symmetry_matrices{22 + 1} = [-0.809017, 0.309017, 0.500000;-0.309017, 0.500000,-0.809017;-0.500000,-0.809017,-0.309017]';
            obj.symmetry_matrices{23 + 1} = [-0.809017, 0.309017,-0.500000; 0.309017,-0.500000,-0.809017;-0.500000,-0.809017, 0.309017]';
            obj.symmetry_matrices{24 + 1} = [-0.309017,-0.500000,-0.809017; 0.500000,-0.809017, 0.309017;-0.809017,-0.309017, 0.500000]';
            obj.symmetry_matrices{25 + 1} = [ 0.000000,-1.000000, 0.000000; 0.000000, 0.000000, 1.000000;-1.000000, 0.000000, 0.000000]';
            obj.symmetry_matrices{26 + 1} = [-0.809017,-0.309017,-0.500000; 0.309017, 0.500000,-0.809017; 0.500000,-0.809017,-0.309017]';
            obj.symmetry_matrices{27 + 1} = [-0.500000,-0.809017,-0.309017; 0.809017,-0.309017,-0.500000; 0.309017,-0.500000, 0.809017]';
            obj.symmetry_matrices{28 + 1} = [-0.500000,-0.809017, 0.309017; 0.809017,-0.309017, 0.500000;-0.309017, 0.500000, 0.809017]';
            obj.symmetry_matrices{29 + 1} = [-0.809017,-0.309017, 0.500000; 0.309017, 0.500000, 0.809017;-0.500000, 0.809017,-0.309017]';
            obj.symmetry_matrices{30 + 1} = [-1.000000, 0.000000, 0.000000; 0.000000, 1.000000, 0.000000; 0.000000, 0.000000,-1.000000]';
            obj.symmetry_matrices{31 + 1} = [-0.500000, 0.809017, 0.309017;-0.809017,-0.309017,-0.500000;-0.309017,-0.500000, 0.809017]';
            obj.symmetry_matrices{32 + 1} = [-0.309017, 0.500000,-0.809017;-0.500000,-0.809017,-0.309017;-0.809017, 0.309017, 0.500000]';
            obj.symmetry_matrices{33 + 1} = [ 0.309017,-0.500000,-0.809017;-0.500000,-0.809017, 0.309017;-0.809017, 0.309017,-0.500000]';
            obj.symmetry_matrices{34 + 1} = [ 0.500000,-0.809017, 0.309017;-0.809017,-0.309017, 0.500000;-0.309017,-0.500000,-0.809017]';
            obj.symmetry_matrices{35 + 1} = [ 0.000000, 0.000000, 1.000000;-1.000000, 0.000000, 0.000000; 0.000000,-1.000000, 0.000000]';
            obj.symmetry_matrices{36 + 1} = [ 0.309017, 0.500000,-0.809017; 0.500000,-0.809017,-0.309017;-0.809017,-0.309017,-0.500000]';
            obj.symmetry_matrices{37 + 1} = [ 0.809017,-0.309017,-0.500000; 0.309017,-0.500000, 0.809017;-0.500000,-0.809017,-0.309017]';
            obj.symmetry_matrices{38 + 1} = [ 0.809017,-0.309017, 0.500000;-0.309017, 0.500000, 0.809017;-0.500000,-0.809017, 0.309017]';
            obj.symmetry_matrices{39 + 1} = [ 0.309017, 0.500000, 0.809017;-0.500000, 0.809017,-0.309017;-0.809017,-0.309017, 0.500000]';
            obj.symmetry_matrices{40 + 1} = [ 0.000000, 1.000000, 0.000000; 0.000000, 0.000000,-1.000000;-1.000000, 0.000000, 0.000000]';
            obj.symmetry_matrices{41 + 1} = [ 0.809017, 0.309017, 0.500000;-0.309017,-0.500000, 0.809017; 0.500000,-0.809017,-0.309017]';
            obj.symmetry_matrices{42 + 1} = [ 0.500000, 0.809017, 0.309017;-0.809017, 0.309017, 0.500000; 0.309017,-0.500000, 0.809017]';
            obj.symmetry_matrices{43 + 1} = [ 0.500000, 0.809017,-0.309017;-0.809017, 0.309017,-0.500000;-0.309017, 0.500000, 0.809017]';
            obj.symmetry_matrices{44 + 1} = [ 0.809017, 0.309017,-0.500000;-0.309017,-0.500000,-0.809017;-0.500000, 0.809017,-0.309017]';
            obj.symmetry_matrices{45 + 1} = [ 1.000000, 0.000000, 0.000000; 0.000000,-1.000000, 0.000000; 0.000000, 0.000000,-1.000000]';
            obj.symmetry_matrices{46 + 1} = [-0.500000, 0.809017, 0.309017; 0.809017, 0.309017, 0.500000; 0.309017, 0.500000,-0.809017]';
            obj.symmetry_matrices{47 + 1} = [-0.309017, 0.500000,-0.809017; 0.500000, 0.809017, 0.309017; 0.809017,-0.309017,-0.500000]';
            obj.symmetry_matrices{48 + 1} = [ 0.309017,-0.500000,-0.809017; 0.500000, 0.809017,-0.309017; 0.809017,-0.309017, 0.500000]';
            obj.symmetry_matrices{49 + 1} = [ 0.500000,-0.809017, 0.309017; 0.809017, 0.309017,-0.500000; 0.309017, 0.500000, 0.809017]';
            obj.symmetry_matrices{50 + 1} = [ 0.000000, 0.000000, 1.000000; 1.000000, 0.000000, 0.000000; 0.000000, 1.000000, 0.000000]';
            obj.symmetry_matrices{51 + 1} = [ 0.309017, 0.500000,-0.809017;-0.500000, 0.809017, 0.309017; 0.809017, 0.309017, 0.500000]';
            obj.symmetry_matrices{52 + 1} = [ 0.809017,-0.309017,-0.500000;-0.309017, 0.500000,-0.809017; 0.500000, 0.809017, 0.309017]';
            obj.symmetry_matrices{53 + 1} = [ 0.809017,-0.309017, 0.500000; 0.309017,-0.500000,-0.809017; 0.500000, 0.809017,-0.309017]';
            obj.symmetry_matrices{54 + 1} = [ 0.309017, 0.500000, 0.809017; 0.500000,-0.809017, 0.309017; 0.809017, 0.309017,-0.500000]';
            obj.symmetry_matrices{55 + 1} = [ 0.000000, 1.000000, 0.000000; 0.000000, 0.000000, 1.000000; 1.000000, 0.000000, 0.000000]';
            obj.symmetry_matrices{56 + 1} = [ 0.809017, 0.309017, 0.500000; 0.309017, 0.500000,-0.809017;-0.500000, 0.809017, 0.309017]';
            obj.symmetry_matrices{57 + 1} = [ 0.500000, 0.809017, 0.309017; 0.809017,-0.309017,-0.500000;-0.309017, 0.500000,-0.809017]';
            obj.symmetry_matrices{58 + 1} = [ 0.500000, 0.809017,-0.309017; 0.809017,-0.309017, 0.500000; 0.309017,-0.500000,-0.809017]';
            obj.symmetry_matrices{59 + 1} = [ 0.809017, 0.309017,-0.500000; 0.309017, 0.500000, 0.809017; 0.500000,-0.809017, 0.309017]';           
          else
            error('Icosahderal symmetry may be I or I2');
          end
         
         otherwise
           error('Only CX, DX, O, I(2) symmetry is implemented');
       end

       
   
        
  
     
%      switch symmetry_type
%        case 0

%      end
     end

      
     
   end
   
   function [  ] = delete(obj)
     
      free_cuda_objects(obj);
      
   end
   
   function [ ] = free_cuda_objects(obj)
      if ~isempty(obj.texObject) && ~isempty(obj.cuArray)
%       fprintf('Deleting the texture obj, it is (%d) a uint64 (%d)\n',isa(obj.texObject,'uint64'),isa(obj.cuArray,'uint64'));
        mexXform3d(obj.texObject, obj.cuArray);  
        obj.texObject = '';
        obj.cuArray = '';
      end
   end
   
  end
end

