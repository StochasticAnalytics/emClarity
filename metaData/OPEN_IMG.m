function [ vol ] = OPEN_IMG(precision, filename, varargin)

    vol = getVolume(MRCImage(filename, precision), varargin{:});
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