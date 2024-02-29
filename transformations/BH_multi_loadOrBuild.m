function [ IMG_OUT ] = BH_multi_loadOrBuild(tomoName,  ...
                                            mapBackIter, ...
                                            SAMPLING, ...
                                            gpuIDX,...
                                            varargin)
%Check to see if a cached binned image exists, either load or reconstruct
%   Switched to using imod's newstack and binvol to create binning and
%   removed inline binning from my workflow.

IMG_OUT = '';

if nargin > 4
  flgLoad = varargin{1};
end

recon = '';
if nargin > 5
  recon = varargin{2};
end

ali = 'ali';
super_sample = ''; % not used
expand_lines = '';
if nargin > 6
  if ~isempty(varargin{3})
    super_sample = varargin{3};
  end
end

!mkdir -p cache

nameSplit = strsplit(tomoName,'_');
tomoName = strjoin(nameSplit(1:end-1),'_');
tomoIdx = EMC_str2double(nameSplit{end});

checkStack = sprintf('%sStacks/%s_ali%d.fixed',ali,tomoName,mapBackIter+1);

if isempty(recon)
  % Otherwise name is varargin 8 from mapBack
  recon = sprintf('cache/%s_%d_bin%d.rec',tomoName,tomoIdx,SAMPLING);
end


if SAMPLING > 1
  stack = sprintf('cache/%s_ali%d_bin%d.fixed',tomoName,mapBackIter+1,SAMPLING);
  if ~exist(stack, 'file')
    BH_multi_loadOrBin(checkStack, SAMPLING, 2, true); %%%%% med filt flag
  end
else
  stack = sprintf('%sStacks/%s_ali%d.fixed',ali,tomoName,mapBackIter+1);
end

if exist(recon,'file')
  header = getHeader(MRCImage(stack,0));
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
        IMG_OUT = OPEN_IMG('single', m);
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

end
