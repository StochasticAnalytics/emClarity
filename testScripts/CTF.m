classdef CTF < handle
  %UNTITLED Summary of this class goes here
  %   Detailed explanation goes here
  
  properties
    img_nX;
    img_nY;
    img_nZ;
    fou_nX;
    pixelSize;
    coord_grids = cell(3,1);
    useGPU;
    ctf_img = '';
    ctf_abs = '';
    ctf_sq = '';
    
  end
  
  methods
    function obj = CTF(SIZE,pixelSize,useGPU)
      %UNTITLED Construct an instance of this class
      %   Detailed explanation goes here
      obj.pixelSize = pixelSize;
      obj.useGPU = useGPU;
      obj.img_nX = SIZE(1); obj.img_nY = SIZE(2);
      if numel(SIZE) == 3
        obj.img_nZ = SIZE(3);
      else
        obj.img_nZ = 1;
      end
      obj.fou_nX = floor(obj.img_nX/2) + 1;
      
      if strcmp(useGPU,'GPU')
        % Maybe add flag for shifted coordinates
        [ obj.coord_grids{1},obj.coord_grids{3},~,~,~,~ ] = ...
          BH_multi_gridCoordinates( ...
          [obj.img_nX,obj.img_nY,obj.img_nZ],...
          'Cylindrical','GPU', {'none'},1,0,0);
      elseif strcmp(useGPU,'cpu')
        % Maybe add flag for shifted coordinates
        [ obj.coord_grids{1},obj.coord_grids{3},~,~,~,~ ] = ...
          BH_multi_gridCoordinates( ...
          [obj.img_nX,obj.img_nY,obj.img_nZ],...
          'Cylindrical','cpu', {'none'},1,0,0);
      else
        error('useGPU must be GPU or cpu string');
      end
      obj.coord_grids{2} = 0;
      % When half transforms are integrated, update multi_gridCoords to do
      % this directly rather than making the full grid then trimming.
      obj.coord_grids{1} = obj.coord_grids{1}(1:obj.fou_nX,:,:)./pixelSize;
      obj.coord_grids{3} = obj.coord_grids{3}(1:obj.fou_nX,:,:);
      
    end
    
    function obj = new_img(obj,defocusVector,CS,WL,AMP,DC_SCALE,varargin)
      %METHOD1 Summary of this method goes here
      %   Detailed explanation goes here
      if nargin > 6
        obj.ctf_img = BH_ctfCalc(obj.coord_grids,CS,WL,defocusVector,...
          [obj.img_nX,obj.img_nY,obj.img_nZ],...
          AMP,DC_SCALE,varargin{1});
      else
        obj.ctf_img = BH_ctfCalc(obj.coord_grids,CS,WL,defocusVector,...
          [obj.img_nX,obj.img_nY,obj.img_nZ],...
          AMP,DC_SCALE);
      end
      
      
    end
    
    function ctf_product = multiply(obj,otherImg,varargin)
      
      % make sure there is an image
      if isempty(obj.ctf_img)
        error('There is no ctf_img to take the absolute value of!');
      end
      
      if nargin > 2
        form = varargin{1};
      else
        form = 'none';
      end
      
      switch form
        case 'none'
          ctf_product = obj.ctf_img .* otherImg;
        case 'abs'
          if isempty(obj.ctf_abs)
            obj.ctf_abs = abs(obj.ctf_img);
          end
          ctf_product = obj.ctf_abs .* otherImg;
        case 'sq'
          if isempty(obj.ctf_sq)
            obj.ctf_sq = abs(obj.ctf_img).^2;
          end
          ctf_product = obj.ctf_sq .* otherImg;
          
          
        otherwise
          error('form must be none, abs, or sq, not %s\n',form)
      end
      
    end
    
  end
end

