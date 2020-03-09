function [ ] = SAVE_IMG( vol, varargin)
%Overloaded MRCImage method SAVE_IMG 
%   This exists so that the call SAVE_IMG(MRCImage( vol)) can be avoided
%   (using the MRCImage.SAVE_IMG). This also permits passing in a gpuArray
%   which in turn allows much faster calcs on the stats.

  imgMin = gather(min(vol(:)));
  imgMax = gather(max(vol(:)));
  imgMean = gather(mean(vol(:)));
  imgRMS = gather(rms(vol(:)));
  
  mRCImage = MRCImage(gather(vol));
   
  mRCImage.header.minDensity = imgMin;
  mRCImage.header.maxDensity = imgMax;

  mRCImage.header.meanDensity = imgMean;
  mRCImage.header.densityRMS = imgRMS;
%   mRCImage.writeHeader(mRCImage);   Changed SAVE_IMG to write out the header no matter what.   

  switch length(varargin)
    case 0
      SAVE_IMG(mRCImage);
    case 1
      SAVE_IMG(mRCImage, varargin{1});
    case 2
      SAVE_IMG(mRCImage, varargin{1},varargin{2});
    case 3
      SAVE_IMG(mRCImage, varargin{1},varargin{2},varargin{3});
    case 4
      SAVE_IMG(mRCImage, varargin{1},varargin{2},varargin{3},varargin{4});
    case 5
      SAVE_IMG(mRCImage, varargin{1},varargin{2},varargin{3},varargin{4},varargin{5});
    otherwise
      error('More than 5 additional args in SAVE_IMG');
  end
  
  
end

