function [tmpCache, flgCleanCache, CWD] = EMC_setup_tmp_cache(emc_fastScratchDisk, existing_tmpCache, caller_name, add_trailing_slash)

    % After parseParameter file
    % emc.fastScratchDisk= '' || getenv('MCR_CACHE_ROOT');

    tmpCache = emc_fastScratchDisk;
    flgCleanCache = false;
    if isempty(existing_tmpCache)
        % Default condition we need to setup the cache which is now coming
        % from the parser and can be empty or a path to the local MCR on ramdisk
        if isempty(tmpCache)
            % Local project directory, we don't wan to clean the cache as this will conflict
            % with other processes
            tmpCache='cache';
        end
    else        
        % We are using an existing cache directory, we need to make sure it exists
        if isfolder(existing_tmpCache)
            tmpCache = existing_tmpCache;
        else
            error('The existing tmpCache %s does not exist',existing_tmpCache);
        end
    end

    [ tmpCache ] = strip_trailing_slash(tmpCache);
    [filepath,name,ext] = fileparts(tmpCache);
    % Check if the path is local or remote, this logic only works if the trailing slash is gone,
    % eg. /to/path, '', '' = fileparts(/to/path/)
    if (isempty(name) && isempty(filepath))
        error('The tmpCache (%s) should not be empty at this point in the code.', tmpCache)
    end

    if isempty(name)
        % We must be in the local project directory, add this little check
        if ~isfolder('fixedStacks')
            % Should work for soft links too
            error('The fixedStacks directory does not exist in the current directory, %s',pwd);
        end
        tmpCache = 'cache';
        CWD = '';
    else
        if ~strcmp(name,'cache')
           tmpCache = fullfile(filepath,name, 'cache'); 
        end
        CWD = sprintf('%s/',pwd);
    end

    system(sprintf('mkdir -p %s', tmpCache));

    % Check for legal names
    switch caller_name
        case 'ctf3d'
            tmpCache = fullfile(tmpCache, 'ctf3d');
            system(sprintf('mkdir -p %s', tmpCache));
            if ~isfolder(tmpCache)
                error('The tmpCache %s does not exist',tmpCache);
            end
        case 'tomoCPR'
            % Nothing to do here
            
        case 'cisTEM'
            tmpCache = fullfile(tmpCache, 'to_cisTEM');
            system(sprintf('mkdir -p %s', tmpCache));
        otherwise
            error('This function is not allowed to be called by %s',caller_name);
    end

    % Finally check if the program asks for a trailing slash.
    if (add_trailing_slash)
        tmpCache = sprintf('%s/',tmpCache);
    end
    fprintf('tmpCache is %s\n',tmpCache);

end

function [ cleaned ] = strip_trailing_slash(input)

  % Check for a trailing slash
  slashCheck = strsplit(input,'/');
  if isempty(slashCheck{end})
    cleaned = strjoin(slashCheck(1:end-1),'/');
  else
    cleaned = input;
  end

end