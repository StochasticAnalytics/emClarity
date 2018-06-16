%MRCImage        MRCImage Constructor 
%
%   mRCImage = MRCImage
%   mRCImage = MRCImage(filename)
%   mRCImage = MRCImage(filename, flgLoad)
%   mRCImage = MRCImage(header)
%   mRCImage = MRCImage(header, filename)
%   mRCImage = MRCImage(header, filename, flgLoad)
%   mRCImage = MRCImage(volume)
%   mRCImage = MRCImage(MRCImage)
%   mRCImage = MRCImage(MRCImage, fileName)
%
%   mRCImage    The constructed MRCImage object
%
%   fileName    The name of an MRC image file to use in initializing the
%               object.
%
%   flgLoad     A flag specifying whether to load the volume into memory.
%               (default: 1, load volume).
%
%   header      A header for creating an empty (zeroed) volume.
%
%
%   MRCImage constructs an MRCImage object and optionally initializes it
%   with the specified MRC image file.
%
%   Bugs: none known
%
% This file is part of PEET (Particle Estimation for Electron Tomography).
% Copyright 2000-2014 The Regents of the University of Colorado & BL3DEMC:
%           The Boulder Laboratory For 3D Electron Microscopy of Cells.
% See PEETCopyright.txt for more details.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  $Author: John Heumann $
%
%  $Date: 2014/01/13 20:00:38 $
%
%  $Revision: 6b413b88334c $
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mRCImage = MRCImage(varargin)
  
  % This sometimes creates weird errors, just put the default function inside this 
  % function. 20170621 BAH
  %%% Create a default MRCImage image structure
  %%%mRCImage = default;
  [ mRCImage ] = init_default_struct();
  mRCImage = class(mRCImage, 'MRCImage');    

  % Default constructor
  if length(varargin) < 1
    return;
  end

  % Decide which type of constructor is being called
  % - If the first argument is a string it is the name of the MRCImage file.
  % - if the first argument is a struct it is (partial) header; create a
  %   zero'd volume.
  % - If numeric, it is the volume to use.
  % - If the first argument is another MRCImage then do a copy.
  if isa(varargin{1}, 'char')
    if length(varargin) < 2
      
      mRCImage = open(mRCImage, varargin{1});
    else
      %fprintf('opening\n')
      
      mRCImage = open(mRCImage, varargin{1}, varargin{2});
     % fprintf('done opening\n');
    end

  elseif isa(varargin{1}, 'struct')
    if nargin > 2
      flgLoad = varargin{3};
      filename = varargin{2};
    else
      flgLoad = 0; % % The PEET default is to laod the volume into the MRC Object
                       % There is probably some good reason for this, but I
                       % don't see it. Double memory and slower. BAH
                       % 20171203
      filename = '';
    end
    mRCImage = emptyVolume(mRCImage, varargin{1}, filename, flgLoad);
    
  elseif isa(varargin{1}, 'numeric')
    
    if nargin > 1
      % Scale to new data type
      varargin{1} = BH_multi_statScale(varargin{1},varargin{2});
    end
    
    mRCImage = setVolumeAndHeaderFromVolume(mRCImage, varargin{1});
    
  else
    if nargin < 2
      % Direct copy
      mRCImage = varargin{1};
    else
      % Copy an existing MRCImage (file and object) and give it a new filename
      destFilename = varargin{2};
      if destFilename(1) ~= '/'
        workingDir = cd;
        destFilename = [workingDir '/' destFilename];
      end
      srcFilename = varargin{1}.filename;
      [stat message] = copyfile(srcFilename,  destFilename);
      if ~ stat
        disp(message);
        PEETError('Unable to copy file!');
      end
      mRCImage = varargin{1};

      % Reset the file descriptor
      mRCImage.fid = [];

      % Open the copied file in the same form as the source MRCImage
      mRCImage = open(mRCImage, destFilename, mRCImage.flgVolume);   
    end
  end

end 

function [ mRCImage ] = init_default_struct()

  mRCImage.fid = [];
  mRCImage.filename = [];
  mRCImage.endianFormat = 'ieee-le';
  mRCImage.type = 'BL3DFS';
  mRCImage.version = '1.0';

  mRCImage.dataIndex = -Inf;
  mRCImage.volume = [];
  mRCImage.flgVolume = 0;

  mRCImage.header.nX = -Inf;
  mRCImage.header.nY = -Inf;
  mRCImage.header.nZ = -Inf;
  mRCImage.header.mode = -Inf;
  mRCImage.header.nXStart = -Inf;
  mRCImage.header.nYStart = -Inf;
  mRCImage.header.nZStart = -Inf;
  mRCImage.header.mX = -Inf;
  mRCImage.header.mY = -Inf;
  mRCImage.header.mZ = -Inf;
  mRCImage.header.cellDimensionX = -Inf;
  mRCImage.header.cellDimensionY = -Inf;
  mRCImage.header.cellDimensionZ = -Inf;
  mRCImage.header.cellAngleX = -Inf;
  mRCImage.header.cellAngleY = -Inf;
  mRCImage.header.cellAngleZ = -Inf;
  mRCImage.header.mapColumns = 1;
  mRCImage.header.mapRows = 2;
  mRCImage.header.mapSections = 3;
  mRCImage.header.minDensity = -Inf;
  mRCImage.header.maxDensity = -Inf;
  mRCImage.header.meanDensity = -Inf;
  mRCImage.header.spaceGroup = -Inf;
  mRCImage.header.nSymmetryBytes = -Inf;
  mRCImage.header.nBytesExtended = -Inf;
  mRCImage.header.creatorID = -Inf;
  mRCImage.header.nBytesPerSection = -Inf;
  mRCImage.header.serialEMType = -Inf;
  mRCImage.header.imodStamp = -Inf;
  mRCImage.header.imodFlags = -Inf;
  mRCImage.header.idtype = -Inf;
  mRCImage.header.lens = -Inf;
  mRCImage.header.ndl = -Inf;
  mRCImage.header.nd2 = -Inf;
  mRCImage.header.vdl = -Inf;
  mRCImage.header.vd2 = -Inf;
  mRCImage.header.tiltAngles = [];
  mRCImage.header.extra = [];
  mRCImage.header.xOrigin = -Inf;
  mRCImage.header.yOrigin = -Inf;
  mRCImage.header.zOrigin = -Inf;
  mRCImage.header.map = '';
  mRCImage.header.machineStamp = '';
  mRCImage.header.densityRMS = -Inf;
  mRCImage.header.nLabels = -Inf;
  mRCImage.header.labels = [blanks(80); blanks(80); 
                      blanks(80); blanks(80); 
                      blanks(80); blanks(80);
                      blanks(80); blanks(80);
                      blanks(80); blanks(80) ];

  mRCImage.extended = [];

  mRCImage.forceWriteByteMode = [];

end
