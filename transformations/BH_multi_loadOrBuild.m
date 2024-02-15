function [ IMG_OUT, reconGeom ] = BH_multi_loadOrBuild( tomoName,  ...
                                                        rCoords, ...
                                                        mapBackIter, ...
                                                        SAMPLING, ...
                                                        gpuIDX,...
                                                        reconScaling, ...
                                                        varargin)
%Check to see if a cached binned image exists, either load or reconstruct
%   Switched to using imod's newstack and binvol to create binning and
%   removed inline binning from my workflow.

IMG_OUT = '';

% If gpuIDX is negative, this is from template matching, so allow
% reconstruction of non-CTF stack. Otherwise throw error if existing
% reconstruction is not available.
if gpuIDX < 0
  doRecon = 1;
  gpuIDX = abs(gpuIDX);
else
  doRecon = 0;
end


if nargin > 6
  flgLoad = varargin{1};
end

recon = '';
if nargin > 7
  recon = varargin{2};
end

ctf = '';
ali = 'ali';
super_sample = '';
expand_lines = '';
if nargin > 8
  if ~isempty(varargin{3})
    super_sample = varargin{3};
  end
  
  %   expand_lines = varargin{4}
end

!mkdir -p cache

nameSplit = strsplit(tomoName,'_');
tomoName = strjoin(nameSplit(1:end-1),'_')
tomoNumber = EMC_str2double(nameSplit{end})

rCoords = rCoords ./ SAMPLING;
% fix is like floor but rounds towards zero, not sure why I'm doing this here anymore.
rCoords(1:4) = fix(rCoords(1:4));



checkStack = sprintf('%sStacks/%s_ali%d%s.fixed',ali,tomoName,mapBackIter+1,ctf);

if isempty(recon)
  % Otherwise name is varargin 8 from mapBack
  recon = sprintf('cache/%s_%d_bin%d.rec',tomoName,tomoNumber,SAMPLING);
end


if SAMPLING > 1
  
  stack = sprintf('cache/%s_ali%d%s_bin%d.fixed',tomoName,mapBackIter+1,ctf,SAMPLING);
  if ~exist(stack, 'file')
    BH_multi_loadOrBin(checkStack,SAMPLING, 2); %%%%% med filt flag
  end
  
else
  
  stack = sprintf('%sStacks/%s_ali%d%s.fixed',ali,tomoName,mapBackIter+1,ctf);
  
end


if exist(recon,'file') || ~doRecon
  header = getHeader(MRCImage(stack,0));
  [ reconGeom ] = calc_rg(  header, rCoords );
elseif (doRecon)
  
  if (mapBackIter)
    TLT = sprintf('mapBack%d/%s_ali%d_ctf.tlt',mapBackIter,tomoName,...
      mapBackIter);
    LOCAL = sprintf('mapBack%d/%s_ali%d_ctf.local',mapBackIter,tomoName, ...
      mapBackIter);
  else
    TLT = sprintf('fixedStacks/%s.tlt',tomoName);
    LOCAL = sprintf('fixedStacks/%s.local',tomoName);
  end
  
  if exist(LOCAL,'file')
    flgLocal = 1;
  else
    fprintf('Did not find local alignment information at %s\n',LOCAL);
    flgLocal = 0;
  end
  
  % check to see if the binned stack exists and is readable
  [initialCheckFail,~] = system(sprintf('header %s',stack));
  
  if (initialCheckFail)
    % See if the file exists but is being written by another process.
    if exist( stack, 'file')
      % It is there but possibly being written, run imod wait which throws an error if not growing.
      BH_imodWait(stack)
    else
      
      error('Did not find the full aligned stack at %s\n',stack);
      
    end
    
  else
    
    fprintf('Reconstructing from cached stack %s\n', stack);
    
  end
  
  header = getHeader(MRCImage(stack));
  
  if exist(recon,'file')
    fprintf('Using cached file %s\n', recon);
    [ reconGeom ] = calc_rg(  header, rCoords );
  else
    
    
    % Check that no out of bounds occur on slices
    
    if (rCoords(2) == 0)
      fprintf('shifting slices up 1\n')
      rCoords(2) = 1;
    end
    if (rCoords(3) > header.nY)
      if (header.nY - rCoords(3) > 2)
        fprintf('Slice index is too high.\n')
      else
        rCoords(3) = header.nY-1;
      end
    end
    
    [ reconGeom ] = calc_rg(  header, rCoords );
    
    rCMD = sprintf(['-input %s -output %s -TILTFILE %s -UseGPU %d ', ...
      '-WIDTH %d -SLICE %d,%d -THICKNESS %d -SHIFT %f,%f '],...
      stack, recon, TLT, gpuIDX, rCoords(1:6));
    
    
    % Explicitly set Radial to Nyquist
    if (flgLocal)
      rCMD = [rCMD sprintf('-LOCALFILE %s -RotateBy90 -RADIAL 0.5,.05 -MODE 2 -SCALE 0,%f',LOCAL,reconScaling)];
    else
      rCMD = [rCMD sprintf('-RotateBy90 -RADIAL 0.5,.05 -MODE 2 -SCALE 0,%f',reconScaling)];
    end
    
    if system('which tilt')
      error('Did not find IMOD tilt funciton on path')
    else
      fprintf('Reconstructing from newly cached stack %s\n', stack);
      fprintf('tilt %s %s %s\n',rCMD,super_sample,expand_lines)
      system(sprintf('tilt %s %s %s',rCMD,super_sample,expand_lines));
    end
    
  end
else
  error('An appropriate reconstruction was not found for %s\n',recon)
  
end

if (strcmpi(recon,'tomoCPR'))
  m = '';
else
  failedLoads = 0;
  while failedLoads < 6
    try
      % fprintf('pwd is %s\n', pwd);
      % fprintf(...
      %  'attempting to load %s\n', recon);
      m = MRCImage(sprintf('%s', recon),0);
      % fprintf('Loaded the MRCImage\n');
      if ( flgLoad )
        IMG_OUT = single(getVolume(m));
        % fprintf('Loaded the volume\n');
      else
        IMG_OUT = m;
        % fprintf('Did not load the full Volume\n');
      end
      failedLoads = 6;
    catch
      failedLoads = failedLoads + 1
      if failedLoads < 4
        pause(failedLoads.^3); % 1,8,27 second pauses
      else
        pause(randi(100));
      end
    end
  end
end



  function [ reconGeom ] = calc_rg(  header, rCoords )
    
    % Origin in the binned tilt series
    oY = 1 + floor(header.nY ./ 2);
    % Size in the binned tilt series
    nY = rCoords(3) - rCoords(2) + 1;
    % Origin in the reconstructed area (active shift from origin in tilt series)
    dY = floor(rCoords(2) + nY/2) - oY;
    reconGeom = zeros(2,3);
    % FIXME index 4 and 5 should be 3 and 4
    reconGeom(1,1:3) = [rCoords(1), nY, rCoords(4)];
    % value specify location of origin, but SHIFT in IMOD's tilt takes the
    % location to shift the origin too, so multiply oX by -1. The notion for Z is
    % flipped since imod does reconstruction on a rotated frame. I.e. a positive
    % number shifts the recon "up" which when rotated to the microscope frame is
    % actually "down" (in Z) so no need to multipliy oZ by -1.
    reconGeom(2,1:3)   = round([-1*rCoords(5),dY,rCoords(6)]);
  end
end
