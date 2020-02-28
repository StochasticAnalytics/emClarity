function [ STACK ] = BH_multi_loadAndMaskStack(STACK,TLT,mapBackIter,THICKNESS,PIXEL_SIZE)
%Relative weighting of a tilt-series
%   Weight projections after normalizing stats (mean and unit var)
%   1) For best expected signal, which is assumed to be the fraction of
%   inelestics. Assume 100nm thickness unless an estimate is known from
%   tomoCPR. Use meanfree path for in-elastics of 400nm
%   2) Normal schemes aim to have ~ the same dose at the detector. If this
%   is not true, the projections must be weighted appropriately. For the
%   initial implementation, this will require a "custom" order file which
%   has for each tilt:
%     1 - number in stack
%     2 - order in data collection (from 1) -9999 to skip
%     3 - total dose on *this* tilt (not cummulative)
%     4 - defocus 1
%     5 - defocus 2
%     6 - astigmatism angle
%
%     The final columns can be set to zero to calculate in ctf_estimate,
%     otherwise the program will now just set up the tlt info, resample the
%     stack, and exit.
%
%     % Store the dose per-tilt in col 14 (cummul. dose is in 11)
%
%
%     % STACK - 3d vol stack, or PATH to file
%     % TLT - loaded tilt geometry
%     % THICKNESS - estimate for thickness
%     % PIXEL_SIZE - angstrom
%     %
%     % STACK - Modified stack in mem, or on disk

EDGE_PAD = 64; % min dist to edge for stats calc

if (THICKNESS < 10 || THICKNESS > 1000)
  error('Your sample thickness is likely incorrect, it should be between 10 and 1000 nm, not %d \n',THICKNESS);
end

% Override the thickness for now, if the sample isn't flat, the method to
% calculate it based on subTomogram coordinates in 3d is not correct
THICKNESS = 75;
if ~isnumeric(STACK)

  if samplingRate > 1
    fullStack = sprintf('aliStacks/%s_ali%d.fixed', STACK,mapBackIter+1);
    inputStack = sprintf('cache/%s_ali%d%s_bin%d.fixed',STACK,mapBackIter+1,samplingRate);

    if ~exist(inputStack, 'file')
      BH_multi_loadOrBin(fullStack,-1.*samplingRate,2);
    end

  else
    inputStack = sprintf('aliStacks/%s_ali%d%s.fixed',STACK,mapBackIter+1,suffix);
  end


  STACK = single(getVolume(MRCImage(inputStack)));

end


[d1,d2,d3] = size(STACK);

if isa(STACK,'gpuArray')
  flgOnDevice = 1;
else
  flgOnDevice = 0;
end

bandNyquist = BH_bandpass3d([d1,d2,1],0,0,1,'GPU','nyquistHigh');
meanVariance = 0;
for iPrj = 1:d3
  
  if (flgOnDevice)
    iProjection = STACK(:,:,TLT(iPrj,1));
  else
    iProjection = gpuArray(STACK(:,:,TLT(iPrj,1)));
  end
  
  % Remove any residual outliers
  mask = (abs(iProjection(:)) > mean(iProjection(:))+5*std(iProjection(:)));
  
  iProjection(mask) = randn([sum(mask(:)),1],'single','gpuArray').*(0.5*std(std(iProjection(~mask))));
 
  maxEval = cosd(TLT(iPrj,4)).*(d1/2) + (THICKNESS*10./(PIXEL_SIZE))./2*abs(sind(TLT(iPrj,4)));
  oX = ceil((d1+1)./2);
  iEvalMask = max(EDGE_PAD,floor(oX-maxEval)):min(d1-EDGE_PAD,ceil(oX+maxEval));

  % Remove gradients
  iProjection  = real(ifftn(fftn(iProjection).*bandNyquist));
  % Center under mask
  iProjection = iProjection - mean2(iProjection(iEvalMask,EDGE_PAD:end-EDGE_PAD));
  % Unit variance
  iProjection = iProjection ./ rms(rms(iProjection(iEvalMask,EDGE_PAD:end-EDGE_PAD)));
  
  fractionOfDose = TLT(iPrj,14)/mean(TLT(:,14));
  fractionOfElastics = exp(-1.*THICKNESS/( cosd(TLT(iPrj,4))*400 ));
  fractionOfElastics = fractionOfElastics ./ max(fractionOfElastics(:));
  
  iProjection = iProjection .* (fractionOfDose*fractionOfElastics);

  meanVariance = meanVariance + rms(rms(iProjection(iEvalMask,EDGE_PAD:end-EDGE_PAD)));

  
  
  if (flgOnDevice)
    STACK(:,:,TLT(iPrj,1)) = iProjection;
  else
    STACK(:,:,TLT(iPrj,1)) = gather(iProjection);
  end
  
end % end loop over projections
meanVariance = meanVariance ./ d3;
STACK = STACK ./ gather( meanVariance);

end % end of function

