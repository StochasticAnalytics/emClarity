classdef fourierTransformer < handle
  %Create and store cufft plans, do R2C and C2R with optional normalization.

  
  properties (Access = 'public')

    bandpass; % for now create one internally and just chop off half.
    inputSize = '';
    is2d;
    halfDim;
    halfDimSize;
    padValIn;
    padValOut;
    
%   end
%   
%   properties (Access = 'protected')
 
    normalization_factor;

    plan_FWD = '';
    plan_INV = ''; 
    bpVals = [0,0,0,0];
    bpDoesNotExist = true;
    invTrim; % boolean for the input size reduced dimensino
    phaseCenter = '';
    indexCenterFWD = '';  
    indexCenterINV = ''
    
    useFwdSwapForInverse; % works for even sized images
    
  end
  
  methods
    
    function [obj, ft] = fourierTransformer(inputVol)
      
      % Must be single and on gpu
      if ~isa(inputVol(1),'gpuArray')
        inputVol = gpuArray(single(inputVol));
      else
        inputVol = single(inputVol);
      end
      
      % Assuming a forward FFT if the object is not initialized.
      obj.inputSize = size(inputVol);
      obj.halfDim = 1;
      obj.halfDimSize = floor(obj.inputSize(1)/2)+1;
      

      
      if (mod(obj.inputSize(obj.halfDim),2))
        obj.invTrim = int16(1);
      else
        obj.invTrim = int16(2);
      end
      % Second condition checks for even dimensions which require an
      % explicit ifft mask
      if (sum(mod(obj.inputSize,2)))
        obj.useFwdSwapForInverse = false;
      else
        % All dimensions even, only use a fwd mask for fwd/inv swap
        obj.useFwdSwapForInverse = true;
      end
        
      if (numel(obj.inputSize) > 2)
        obj.is2d = false;
      else
        obj.is2d = true;
      end
      
      [ft] = fwdFFT(obj,inputVol);

      obj.normalization_factor = 1./sqrt(numel(inputVol));
      
    end

    function [ft] = fwdFFT(obj, inputVol, varargin)
      % Vararginr = 
      % 1 - Bandpass filter ( or a 1 )
      % 2 - bool center and standardize to 1
      % 3 - normalize scaling, 0 for none, 1 if only fourier comp, 2 for
      % complete.
      doBandpass = false;
      doCenter = false;
      if nargin > 2
        doNorm   = varargin{1};
        doCenter = varargin{2};       
        if length(varargin) > 2 && ~isempty(varargin{3})
          
          doBandpass = true;     
          obj.makeBandPass(size(inputVol),varargin{3})
          
        end      

      else
        doNorm = 0;
      end
      
      if isempty([obj.plan_FWD,obj.plan_INV])
        [ ft, obj.plan_FWD, obj.plan_INV ] = mexFFT(inputVol,obj.invTrim);
      else
        [ ft ] = mexFFT(inputVol,obj.invTrim,obj.plan_FWD, obj.plan_INV);    
      end
      
      if (doNorm && ~doCenter)
        ft = ft .* (obj.normalization_factor^doNorm);
      end
      if (doBandpass)
 
        ft = ft .* obj.bandpass;
      end
      if (doCenter)
        ft(1) = 0;
        ft = ft ./ sqrt(2.*sum(sum(sum(abs(ft(1:end-obj.invTrim,:,:)).^2))));
      end
      
    end

    function [ft] = invFFT(obj, inputVol, varargin)
      if nargin > 2
        doNorm = varargin{1};
      else
        doNorm = 0;
      end

      [ ft ] = mexFFT(inputVol,obj.invTrim,obj.plan_FWD, obj.plan_INV);
      
      if (doNorm)
        ft = ft .* (obj.normalization_factor^doNorm);
      end
    end
    

    function to_cpu(obj)
      % There are a number of places where I reset the gpuDevice. Use this
      % to protect the properties.
      obj.bandpass = gather(obj.bandpass);
      
    end
    
    function to_GPU(obj)
      % There are a number of places where I reset the gpuDevice. Use this
      % to un-protect the properties.
      obj.bandpass = gpuArray(obj.bandpass);
      
    end
    
    function delete(obj)
      % Destroy the cuda plans. mexFFT knows to do this because the
      % first argument has one element;
      if ~(isempty(obj.plan_FWD) && isempty(obj.plan_INV))
        
          mexFFT(gpuArray(1),obj.invTrim,obj.plan_FWD,obj.plan_INV);
         
      end
    end 
    
    function [inputVol] = shiftStretch(obj, inputVol, shiftXY, Mag, isCentered)
      
        [ dU, dV ] = BH_multi_gridCoordinates(obj.inputSize,'Cartesian','GPU', ...
                                    {'none'},1,isCentered,0,{'halfgrid'});
         inputVol = inputVol .* (Mag.^-2.*exp(-2i.*pi.*(dU.*shiftXY(1)+dV.*shiftXY(2))));
         
         clear dU dV
                                
    end
    
    function [inputVol] = swapPhase(obj, inputVol, direction)
      

      % Create the phase swap indices if needed
      if isempty(obj.phaseCenter) 
        if obj.is2d
          [ obj.phaseCenter, dV ] = BH_multi_gridCoordinates(obj.inputSize,'Cartesian','GPU', ...
                                    {'none'},1,1,0,{'halfgrid'});
          if (obj.inputSize(1) == obj.inputSize(2))
            obj.phaseCenter = exp(-2i.*pi.*obj.halfDimSize.*(obj.phaseCenter+dV));
            clear dU dV
          else
            hX = floor(obj.inputSize(1)/2) + 1;
            hY = floor(obj.inputSize(2)/2) + 1;
            obj.phaseCenter = exp(-2i.*pi.*(hX.*obj.phaseCenter+hY.*dV));
          end
        else
          [ obj.phaseCenter, dV, dW] = BH_multi_gridCoordinates(obj.inputSize,'Cartesian','GPU', ...
                                    {'none'},1,1,0,{'halfgrid'});
                                  obj.inputSize
          if (obj.inputSize(1) == obj.inputSize(2) == obj.inputSize(3))
            obj.phaseCenter = exp(-2i.*pi.*obj.halfDimSize.*(obj.phaseCenter+dV+dW));
            clear dU dV dW        
          else
            hX = floor(obj.inputSize(1)/2) + 1;
            hY = floor(obj.inputSize(2)/2) + 1;
            hZ = floor(obj.inputSize(3)/2) + 1;
            obj.phaseCenter = exp(-2i.*pi.*(hX.*obj.phaseCenter+hY.*dV+hZ.*dW));
            clear dU dV dW
          end
        end
      end
      


      if strcmp(direction,'fwd')
        inputVol = inputVol .* obj.phaseCenter;
      elseif strcmp(direction,'inv')
        inputVol = inputVol ./ obj.phaseCenter;
      else
        error('Direction not fwd or inv: %s',direction);
      end
     
      
    end
    
    function [inputVol] = swapIndexFWD(obj, inputVol)
      
      % TODO setup to take a window to reduce ifftshift on CCC maps where
      % the peak is expected near the center.
      window = 0;
      
      % Create the phase swap indices if needed
      if isempty(obj.indexCenterFWD)         
        obj.indexCenterFWD = BH_fftShift(window,obj.inputSize,1,'halfgrid');
      end
      
      inputVol = inputVol(obj.indexCenterFWD);

    end
    
    function [inputVol] = swapIndexINV(obj, inputVol)
      
      % TODO setup to take a window to reduce ifftshift on CCC maps where
      % the peak is expected near the center.
      window = 0;
      
      if (obj.useFwdSwapForInverse)
        inputVol = obj.swapIndexFWD(inputVol);
      else
        % Create the phase swap indices if needed
        if isempty(obj.indexCenterINV)         
          obj.indexCenterINV = BH_fftShift(window,-1.*obj.inputSize,1,'halfgrid');
        end
        inputVol = inputVol(obj.indexCenterINV);
      end

    end
    
    function [inputVol] = fwdSwap(obj,inputVol)
      
      % For transformations in Fourier space. 
      % Assumed that the input vol is not index swapped. TODO add a check
      inputVol = obj.swapPhase(obj.swapIndexFWD(inputVol),'fwd');
      
    end
    
    function [inputVol] = invSwap(obj,inputVol)
      
      % For transformations in Fourier space

      inputVol = obj.swapIndexINV(obj.swapPhase(inputVol,'inv'));
      
    end
    
   
    
    function makeBandPass(obj, sizeInput, bpValsNew)
      
      doCalc = false;
      % First see if the filter exists
      if (obj.bpDoesNotExist)
        doCalc = true;
      else
        if abs(sum( obj.bpVals - bpValsNew )) > 1e-4      
          doCalc = true;
        end        
      end
        
      
      if (doCalc)        
        obj.bandpass = BH_bandpass3d(sizeInput,bpValsNew(1),...
                                                    bpValsNew(2),...
                                                    bpValsNew(3),...
                                                    'GPU', ...
                                                    bpValsNew(4));
        switch ndims(obj.bandpass)
          case 3
            obj.bandpass = obj.bandpass(1:obj.halfDimSize,:,:);
          case 2
            obj.bandpass = obj.bandpass(1:obj.halfDimSize,:);
          case 1
            obj.bandpass = obj.bandpass(1:obj.halfDimSize);
        end
        
        % Update the properties
        obj.bpDoesNotExist = false;
        obj.bpVals = bpValsNew;
      end
      
      
    end
    
  end
end

