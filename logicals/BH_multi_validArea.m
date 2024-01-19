function [ sizeWINDOW, sizeCALC, sizeMASK, padWINDOW, padCALC ] = ...
  BH_multi_validArea( MASK_SIZE, MASK_RADIUS, scaleCalcSize )
%Calculate the padding and sizes needed to work with a given dimension.
%
%   Input:
%
%   MASK_RADIUS: radius in pixels defining the size of valid area
%
%   SIZE_REF: size of the reference being used, check that it is large enough to
%             be masked properly.
%   Output:
%
%   sizeWINDOW: size to cut out of tomogram to allow generic rotation without
%               extrapolation.
%
%   sizeCALC: padded size to avoid wraparound effects in circular convolution.
%
%   sizeMASK: smallest size to use for mask creation that includes enough room
%             for an appropriate taper. This would be passed to BH_mask3d as the
%             size, or trim to this size and then pass volume to BH_padZeros3d
%             with the taper option.
%
%   padWINDOW: excess dimensions to remove from transformed window to get down
%              to sizeMASK.
%
%   padCALC: size to zero pad for wraparound
%
%   padREF: size cut ref to match sizeMASK
%
%           padding is [ pre X Y Z ; post X Y Z ]
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




sizeMASK = ceil(scaleCalcSize.* max(MASK_SIZE).*[1,1,1]);

% Ensure that the reference being used is large enough to be masked
% appropriately



% I don't remember what I was going for here. This should be simpler.
% % % padWINDOW = ceil((sqrt(mPlusTaper*mPlusTaper')-mPlusTaper));
% % % padWINDOW = [ padWINDOW ; padWINDOW ];
% % % sizeWINDOW = 2.* padWINDOW(1,:) + sizeMASK;
sizeWINDOW = max(sizeMASK).*[1,1,1];
padWINDOW = BH_multi_padVal(sizeMASK,sizeWINDOW);


% Calc padding, each image masked to dimesion 2*maskRad+20
% Also make cubic so that weight mask can be calculated once for each tomo and
% then interpolated. In the future, it would be better to calculate a mask for
% each particle, then this could just be 2.*sizeMask.
%sizeCALC = 2.*sizeMASK;

% sizeCALC = [2,2,2].*max(sizeMASK);

%%%%%%%%%%%%%%%%%%%%%%%%%
% Test not padding by a factor of 2
sizeCALC = scaleCalcSize.*[2,2,2].*max(MASK_RADIUS);

% Find the next largest size for fft
[ sizeCALC ] = BH_multi_iterator( sizeCALC, 'fourier' );

if any((sizeCALC - sizeWINDOW) < 0);
  sizeCALC = sizeWINDOW;
  % Find the next largest size for fft
  [ sizeCALC ] = BH_multi_iterator( sizeCALC, 'fourier' );
end

[ padCALC ] = BH_multi_padVal( sizeMASK, sizeCALC );



end

