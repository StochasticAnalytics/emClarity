classdef alignmentVol < handle
    properties (Access = 'private')
        data;
        input_size;
        pad_val;
    end
    properties (Access = 'public')
    end
    
    methods
        
        function obj = alignmentVol(data, pad_val)
            pause(3)
            fprintf('alignmentVol constructor\n')
            pause(3)
            obj.data = data;
            obj.input_size = size(obj.data);
            obj.pad_val = pad_val;
        end
        
        function [] = isAllocated(obj)
            if isempty(obj.data)
                error('data is not allocated')
            else
                fprintf('data is allocated\n')
            end
        end
        
        function [] = zero_data(obj)
            obj.data = obj.data .* 0;
        end
    end
    
    
    
end

% function [ IMAGE ] = BH_padZerosSimple3d(IMAGE)

%PADLOW, PADTOP, ...
%METHOD, PRECISION, varargin )
%Pad an image volume with zeros.
%
%
%   Input variables:
%
%   IMAGE = 3d image to be padded
%
%   PADLOW = size of padding pre
%
%   PADTOP = size of padding post
%
%   PRECISION = 'single' or 'double'
%
%   METHOD = case sensitive 'GPU' otherwise cpu
%
%   Output variables:BH_bandpass3d.m
%
%   PADDED_IMG = the padded image.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Goals & limitations:
%
%   When doing both pre and post padding with the build in padarray function,
%   two steps are required, this slows things down. Here the whole matrix is
%   allocated once, and much faster than zeros.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   TODO:
%
%   - Test for expected function, and return value (template matching depends)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% pause(3)
% IMAGE = IMAGE .* 0;
% IMAGE(1:7,1:7,1:7) = 1;
% pause(3)
% if isnumeric(PADLOW)
% padLOW = PADLOW ;
% padTOP = PADTOP;
% elseif strcmpi(PADLOW,'fwd')
% padLOW = PADTOP(1,:);
% padTOP = PADTOP(2,:);
% elseif strcmpi(PADLOW,'inv')
% padLOW = PADTOP(3,:);
% padTOP = PADTOP(4,:);
% else
% error('padLOW should be numeric, fwd or inv');
% end


% twoD = 0;
% doRand = false;
% if numel(padLOW) == 2
% twoD = 1;
% padLOW= [padLOW, 0];
% padTOP= [padTOP, 0];
% end

% if nargin > 5
% if isnumeric(varargin{1})
% extrapVal = varargin{1};
% else
% extrapVal = std(IMAGE(:));
% extrapMean = mean(IMAGE(:));
% doRand = 1;
% end
% else
% extrapVal = 0;
% end

% fourierOverSample = 0;
% if nargin > 5 && length(varargin) > 1
% extrapVal = 0;
% % It is assumed that a bandpass is already applied to that the value at
% % Nyquist is ~0, only applies for 2d right now.
% fourierOverSample = 1;
% end

% trimLOW = abs(padLOW .* (padLOW < 0));
% trimTOP = abs(padTOP .* (padTOP < 0));
% padLOW  = padLOW .* (padLOW >= 0);
% padTOP  = padTOP .* (padTOP >= 0);

% % If any pad values are negative, first trim the image.

% try
% IMAGE = IMAGE(1+trimLOW(1):end-trimTOP(1),...
% 1+trimLOW(2):end-trimTOP(2),...
% 1+trimLOW(3):end-trimTOP(3));
% catch
% fprintf('%f %f %f\n,%f %f %f\n',PADLOW,PADTOP);
% fprintf('%f %f %f\n, %f %f %f\n',trimLOW,trimTOP);
% end

% if ismatrix(IMAGE)
% imgSize = [size(IMAGE),1];
% else
% imgSize = size(IMAGE);
% end

% padSize = imgSize + padLOW + padTOP;




% % Optionally taper the edges, rolling over 6 pixels with the first zero
% % outside the original image, in the padding
% % 1.0000    0.9505    0.8117    0.6113    0.3887    0.1883    0.0495

% taper=false;
% if strcmpi(PRECISION, 'singleTaper')
% taper = 0.5+0.5.*cos((((1:7)-1).*pi)./(length((1:7))));
% PRECISION = 'single';
% elseif strcmpi(PRECISION, 'doubleTaper')
% taper= 0.5+0.5.*cos((((1:7)-1).*pi)./(length((1:7))));
% PRECISION = 'double';
% end

% if strcmp(METHOD, 'GPU')
% if strcmpi(PRECISION, 'single')
% if (doRand)
% PADDED_IMG = randn(padSize,'single','gpuArray').*extrapVal+extrapMean;
% else
% PADDED_IMG = zeros(padSize,'single','gpuArray');
% end
% elseif strcmpi(PRECISION, 'double')
% if (doRand)
% PADDED_IMG = randn(padSize,'double','gpuArray').*extrapVal+extrapMean;
% else
% PADDED_IMG = zeros(padSize,'double','gpuArray');
% end
% else
% error('PRECISION must be single or double, not %s', PRECISION)
% end
% else
% if strcmpi(PRECISION, 'single')
% if (doRand)
% PADDED_IMG = randn(padSize,'single').*extrapVal+extrapMean;
% else
% PADDED_IMG = zeros(padSize,'single');
% end
% elseif strcmpi(PRECISION, 'double')
% if (doRand)
% PADDED_IMG = randn(padSize,'double').*extrapVal+extrapMean;
% else
% PADDED_IMG = zeros(padSize,'double');
% end
% else
% error('PRECISION must be single or double, not %s', PRECISION)
% end
% end

% if ( gather(extrapVal) )
% PADDED_IMG = PADDED_IMG + extrapVal;
% end

% if (twoD)
% if (taper)
% [d1,d2,d3] = size(IMAGE);

% IMAGE(:,1:7) = IMAGE(:,1:7) .* repmat(flip(taper),d1,1,d3) + ...
% repmat(flip(extrapVal.*(1-taper)),d1,1,d3);
% IMAGE(1:7,:) = IMAGE(1:7,:) .* repmat(flip(taper)',1,d2,d3) + ...
% repmat(flip(extrapVal.*(1-taper))',1,d2,d3);


% IMAGE(:,end-6:end) = IMAGE(:,end-6:end) .* repmat(taper,d1,1,d3) + ...
%       repmat(extrapVal.*(1-taper),d1,1,d3);
% IMAGE(end-6:end,:) = IMAGE(end-6:end,:) .* repmat(taper',1,d2,d3) + ...
%       repmat(extrapVal.*(1-taper'),1,d2,d3);

% end

% if ( fourierOverSample )
% sX1 = ceil((size(IMAGE,1)+1)./2);
% sX2 = size(IMAGE,1)-sX1-1;
% sY1 = ceil((size(IMAGE,2)+1)./2);
% sY2 = size(IMAGE,2)-sY1-1;

% PADDED_IMG(1:sX1,1:sY1) = IMAGE(1:sX1,1:sY1);
% PADDED_IMG(end-sX2:end,1:sY1) = IMAGE(end-sX2:end,1:sY1);
% PADDED_IMG(1:sX1,end-sY2:end) = IMAGE(1:sX1,end-sY2:end);
% PADDED_IMG(end-sX2:end,end-sY2:end) = IMAGE(end-sX2:end,end-sY2:end);
% else
% PADDED_IMG(padLOW(1)+1: end - padTOP(1), ...
% padLOW(2)+1: end - padTOP(2)) = IMAGE;
% end
% else
% if (taper)
% [d1,d2,d3] = size(IMAGE);

% IMAGE(:,1:7,:) = IMAGE(:,1:7,:) .* repmat(flip(taper),d1,1,d3) + ...
% repmat(flip(extrapVal.*(1-taper)),d1,1,d3);
% IMAGE(1:7,:,:) = IMAGE(1:7,:,:) .* repmat(flip(taper)',1,d2,d3) + ...
% repmat(flip(extrapVal.*(1-taper')),1,d2,d3);
% IMAGE(:,:,1:7) = IMAGE(:,:,1:7) .* repmat(permute(flip(taper),[3,1,2]),d1,d2,1) + ...
% repmat(permute(flip(extrapVal.*(1-taper)),[3,1,2]),d1,d2,1);

% IMAGE(:,end-6:end,:) = IMAGE(:,end-6:end,:) .* repmat(taper,d1,1,d3) + ...
%           repmat(extrapVal.*(1-taper),d1,1,d3);
% IMAGE(end-6:end,:,:) = IMAGE(end-6:end,:,:) .* repmat(taper',1,d2,d3) + ...
%           repmat(extrapVal.*(1-taper'),1,d2,d3);
% IMAGE(:,:,end-6:end) = IMAGE(:,:,end-6:end) .* repmat(permute(taper,[3,1,2]),d1,d2,1) + ...
%           repmat(permute(extrapVal.*(1-taper),[3,1,2]),d1,d2,1);
% end

% if ( fourierOverSample )
% sX1 = ceil((size(IMAGE,1)+1)./2);
% sX2 = size(IMAGE,1)-sX1-1;
% sY1 = ceil((size(IMAGE,2)+1)./2);
% sY2 = size(IMAGE,2)-sY1-1;
% sZ1 = ceil((size(IMAGE,3)+1)./2);
% sZ2 = size(IMAGE,3)-sZ1-1;

% PADDED_IMG(1:sX1,1:sY1,1:sZ1) = IMAGE(1:sX1,1:sY1,1:sZ1);
% PADDED_IMG(end-sX2:end,1:sY1,1:sZ1) = IMAGE(end-sX2:end,1:sY1,1:sZ1);
% PADDED_IMG(1:sX1,end-sY2:end,1:sZ1) = IMAGE(1:sX1,end-sY2:end,1:sZ1);
% PADDED_IMG(end-sX2:end,end-sY2:end,1:sZ1) = IMAGE(end-sX2:end,end-sY2:end,1:sZ1);

% PADDED_IMG(1:sX1,1:sY1,end-sZ2:end) = IMAGE(1:sX1,1:sY1,end-sZ2:end);
% PADDED_IMG(end-sX2:end,1:sY1,end-sZ2:end) = IMAGE(end-sX2:end,1:sY1,end-sZ2:end);
% PADDED_IMG(1:sX1,end-sY2:end,end-sZ2:end) = IMAGE(1:sX1,end-sY2:end,end-sZ2:end);
% PADDED_IMG(end-sX2:end,end-sY2:end,end-sZ2:end) = IMAGE(end-sX2:end,end-sY2:end,end-sZ2:end);
% else

% PADDED_IMG(padLOW(1)+1: end - padTOP(1), ...
% padLOW(2)+1: end - padTOP(2), ...
% padLOW(3)+1: end - padTOP(3)) = IMAGE;
% end

% end

% clear IMAGE
% end of the padZeros3d functions



