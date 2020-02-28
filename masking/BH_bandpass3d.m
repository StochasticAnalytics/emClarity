function [ BANDPASS ] = BH_bandpass3d( SIZE, HIGH_THRESH, HIGH_CUT, LOW_CUT, ...
                                       METHOD, PIXEL_SIZE )
%Create a bandpass filter, to apply to fft of real space 3d images.
%   
%   Input variables:

%   SIZE = size of the image filter is to be applied to : vector, float
%
%   HIGH_THRESH = Percent attenuation of low frequency : float
%
%   HIGH_CUT  = Spatial frequency high-pass back to 100% (A^-1) : float
%
%   LOW_CUT= Spatial frequency low-pass starts to roll off :  float
%
%   METHOD = 'GPU' case specific, create mask on GPU, otherwise on CPU
%
%   PIXEL_SIZE = Sampling frequency : Angstrom/Pixel
%   
%
%   Output variables:
%
%   BANDPASS  = 3d MRC image file, single precision float.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Goals & Limitations:
%
%   Creates a bandpass filter, that falls off smoothly enough to avoid
%   artifacts, while also restricting information to the expected ranges.
%   This is accomplished by oversampling the image to be filtered, so that
%   the fall off can be over a sufficient number of pixels, while still
%   happening over a small range of frequency.
%   
%   Because this is frequency space, the fall off depends on the resolution
%   where it begins. For sampling of 3A/pixel with a cutoff starting at
%   20A^-1 the spatial frequency drops to ~ 17.3 over six pixels if the 
%   frequency rectangle is 256 pixel sq. For a cutoff starting at 10A the 
%   drop over 6 pixels is only to ~9.3
% 
%   For the monolayer work, the largest dimension is ~140 pixels, s.t. 256
%   provides a substantial padding for cross-correlation, and also a
%   reasonably tight window for filtering. The memory requirements are ~
%   67/134 mb for single/double precision, compared to another order of
%   magnitude for 512. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   TODO:
%     -test with gpu flag
%   
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

       
if numel(SIZE) == 2
  SIZE = [SIZE,1];
end

% The following will adjust the apodization window, and is not intended for
% regular users to adjust. The value of 2.0 makes for a nice (soft) fall
% off over ~ 7 pixels. Larger value results in a steeper cutoff.
% Corresponds to the standard deviation in the gaussian roll.
if isnumeric(PIXEL_SIZE) 
  [bSize, highRoll, lowRoll, highCut, lowCut] = calc_frequencies( ...
                        SIZE, HIGH_THRESH, HIGH_CUT, LOW_CUT, PIXEL_SIZE );
else
  bSize = SIZE;
  if strcmpi(PIXEL_SIZE,'nyquistHigh')
      highCut = 7/min(bSize(bSize>1));
      highThresh = 1e-6;
      highRoll = sqrt((-1.*highCut.^2)./(2.*log(highThresh)));
  else
      highCut = 0; highThresh = 0; highRoll = 0;
  end
  lowRoll = 1.5 .* (1.0./min(bSize(bSize>1)));
  lowCut = 0.485+LOW_CUT;
end   

gaussian = @(x,m,s) exp( -1.*(x-m).^2 ./ (2.*s.^2) );
 


% initialize window of appropriate size
if strcmp(METHOD, 'GPU')
  mWindow(bSize(1),bSize(2),bSize(3)) = gpuArray(single(0)); 
else
  mWindow(bSize(1),bSize(2),bSize(3)) = single(0);
end
mWindow = mWindow + 1;

%%%% This is ~ 120x faster than = gpuArray(ones(bSize, 'single'));
% initialize nd grids of appropriate size


[ radius,~,~,~,~,~] = BH_multi_gridCoordinates( bSize, 'Cartesian', METHOD, ...
                                                  {'single',...
                                                  [1,0,0;0,1,0;0,0,1],...
                                                  [0,0,0]','forward',1,1}, ...
                                                  1, 0, 1 );


% Calc lowpass filter
mWindow = ...
   ( (radius < lowCut) .* mWindow + ...
     (radius >= lowCut) .* gaussian(radius, lowCut, lowRoll) );

% Breaks Hermitian symmetry
% % % Don't randomize the high-pass since the signal it modulates is very
% % % strong.
% % mWindow = BH_multi_randomizeTaper(mWindow);
% Add in high pass if required

if highCut ~= 0
  mWindow = (radius <= highCut) .* gaussian(radius, highCut, highRoll) + ...
         (radius > highCut) .* mWindow;
end   
   
%mWindow((mWindow<= 10^-8)) = 0;   
BANDPASS = mWindow;
clearvars -except BANDPASS
                           
end % end of BH_mask3d function.


function [bSize, highRoll, lowRoll, highCut, lowCut] = calc_frequencies(...
                         SIZE, HIGH_THRESH, HIGH_CUT, LOW_CUT, PIXEL_SIZE )

bSize = SIZE;

% Check that the value makes sense.
if ( 0 > HIGH_THRESH || HIGH_THRESH >= 1)
  error('HIGH_THRESH must be between 0 and 1, not %f', HIGH_THRESH)
end



% Translate boundries from A^-1 to cycles/pixel. A value of zero means no 
% high pass filter.
if HIGH_CUT ~= 0
  highCut = PIXEL_SIZE ./ HIGH_CUT;
else
  highCut = 0;
end
lowCut  = PIXEL_SIZE ./ LOW_CUT;

% fixed lowpass roll off, cycles/pix depends on dimension of image
% if lowCut is negative, indicates a "SIRT like" lowpass and filter rolls
% from 1 at abs(lowCut) to 10^-8 at 20A

if (lowCut > 0)
  if (bSize(3) == 1)
    lowRoll = 2.0 .* (1.0./min(bSize(1:2)));
  else
    lowRoll = 2.0 .* (1.0./min(bSize));
  end
else
  lowCut = abs(lowCut);
  lowEND = 0.5;%PIXEL_SIZE ./ 20;
  lowRoll = sqrt((-1.*(lowEND-lowCut).^2)./(2.*log(10^-3)));
end

% calc the highpass roll off
if HIGH_CUT ~= 0
  highRoll = sqrt((-1.*highCut.^2)./(2.*log(HIGH_THRESH)));
else
  highRoll = 0;
end

                      
                      
end % end of calc_frequencies function.



