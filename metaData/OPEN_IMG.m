function [ vol ] = OPEN_IMG(precision, filename, varargin)

    if isa(filename, 'MRCImage')
        vol = getVolume(filename, varargin{:});
    else 
        vol = getVolume(MRCImage(filename), varargin{:});
    end

    if strcmp(precision, 'single')
        vol = single(vol);
    elseif strcmp(precision, 'double')
        vol = double(vol);
    elseif strcmp(precision, 'half')
        vol = half(vol);
    else
        error('Unknown precision');
    end

end