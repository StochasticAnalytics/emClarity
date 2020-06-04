function [ output_pos, nBeads ] = BH_fitBeads(pixelSize,beadSize,sampling,img_name,input_name,output_name)

pixel_radius = beadSize ./ pixelSize ./ sampling ./ 2;

% Oversample the cross correlation map by this factor to get sub-pixel
% fits using sinc interp.
padBy = 5;

% Factor descrbing 3 gaussians that are summed to make the reference.
  g1 = 2.5; % number of stdDev the central gaussian should fall by the bead radius defined below
  %s1 = pixelRadius. Set in loop so fractions of pixel radius can be searched b/c beads are not uniform diameter (but right now only looking at 1)
  s2 = 1.25; % radius for the edge gaussian, I suspect this should be closer to 1 or even 0.5?
  g2 = s2./2; % center the edge gaussian at pixel_radius + g2  
  g3 = 6; % give the central gaussian a more flat top to match what a bead looks like. Abs is taken so odd values okay

show_ref_profile = false;
show_mip = false;
show_results = false;

bgVal = 0;
edgeVal = 4;
beadVal = -4;
  
input_ts = gpuArray(getVolume(MRCImage(img_name)));
system(sprintf('model2point -contour %s %s.txt',input_name,input_name));
input_pos= gpuArray(load(sprintf('%s.txt',input_name)));
output_pos = zeros(size(input_pos),'single');

[d1,d2,nPrjs] = size(input_ts);

meanVect = [0,0,0];
firstLoop = false;

  KERNEL = EMC_gaussianKernel([1,3], 0.5, 'gpu', {});
  
for iPrj = 1:nPrjs
  iPrj
  img = input_ts(:,:,iPrj);
  img = medfilt2(img,[3,3]);

  lCoord = input_pos(:,4) == iPrj -1;

  
  [w] = BH_multi_gridCoordinates([d1,d2],'Cartesian','GPU',{'none'},1,0,1);
  img_derivative = real(ifftn(fftn(img).*w));
  
%       img_derivative = EMC_convn(single(img_derivative), KERNEL);

  img_derivative = img_derivative - mean(img_derivative(:));
  img_derivative = img_derivative ./ rms(img_derivative(:));
  

  r = ceil(pixel_radius.*3);
  avg_size = 2.*r+1 .*[1,1];
  x = input_pos(lCoord,2);
  y = input_pos(lCoord,3);
  xi = floor(x);
  yi = floor(y);
  xf = x - xi - 0.5; % in imod, model point 7.5 is in the "middle" of pixel 7 in emClarity, a pixel is 6.5 - 7.5 with 7 in the middle
  yf = y - yi - 0.5;


  % Note if you put this into parallel, this will not work.
  if (firstLoop)
    % Get an avgerage bead for this projection to determine the intensities in
    % the bead, on the edge and in the background. Then build a reference using
    % these.
    avg_bead = zeros(avg_size,'single','gpuArray');
    n = 0;
    for i = 1:length(x)
     % add edge checking
       avg_bead = avg_bead + img_derivative(xi(i)-r:xi(i)+r,yi(i)-r:yi(i)+r);
       n = n + 1;
    end
    avg_bead = avg_bead ./ n ;

    % Get the mean values for the three main pixel values
    gFit = fitgmdist(gather(avg_bead(:)),3,'RegularizationValue',0.1,'Replicates',20,'Options',statset('MaxIter',500));
    meanVect = gFit.mu;
    [beadVal,beadCoord] = min(meanVect);
    [edgeVal,edgeCoord] = max(meanVect);

    switch (beadCoord+edgeCoord)
      case 3
        bgVal = meanVect(3);
      case 4
        bgVal = meanVect(2);
      case 5 
        bgVal = meanVect(1);
      otherwise
        error('failed to find the correct index for the background mean');
    end
    firstLoop = false;
  end



  % The edge radius should be around 1 pix unless you blur the derivative
  % first.
  if (show_ref_profile)
    v = -6:0.01:6;
    figure, plot(v, beadVal.* exp(-abs((g1.*v./pixel_radius).^g3)) + ....
                    edgeVal.*(exp(-(g1.*(v-g2-pixel_radius)./s2).^2) + exp(-(g1.*(v+g2 + pixel_radius)./s2).^2)) + ...
                    bgVal);
  end


  [t] = BH_multi_gridCoordinates(avg_size,'Cartesian','GPU',{'none'},0,1,1);



  % figure, imshow3D(cat(3,t,avg_bead))  

  mip = zeros([d1,d2],'single','gpuArray');
  padVal = BH_multi_padVal(size(t),[d1,d2]);
  imgFT = fftn(img_derivative);

  % Loop over references of different fractions of the particle radius. For
  % now, just using 1.
  for iRef = [0.92,1.0,1.08]
    s1 = pixel_radius * iRef;

    ref = beadVal.* exp(-abs((g1.*t./s1).^g3)) + ....
          edgeVal.*(exp(-(g1.*(t-g2-s1)./s2).^2)) + ...
                    bgVal;
    ref = BH_padZeros3d(ref,'fwd',padVal,'GPU','single');

    ccf = real(fftshift(ifftn(conj(fftn(ref)).*imgFT)));
    lMip = ccf > mip;
    mip(lMip) = ccf(lMip);
  end

  if (show_mip)
    figure, imshow3D(gather(real(mip)));
  end
  % Get the updated peak positions, with optional over-sampling of the bead
  % an option to pad zeros in half transforms would be nice here. Just use
  % native matlab FFT for now.


  r = ceil(pixel_radius.*1.0);
  rp = padBy.*(2.*r+1).*[1,1];
  padVal = BH_multi_padVal( (2.*r+1).*[1,1] , rp );

  
  maskRadius = 0.5.*(1- 2.5./sampling).*r.*[1,1];
  
  peakMask = BH_mask3d('sphere',(2.*r+1).*[1,1],maskRadius,[0,0],'2d');
  xo = x;
  yo = y;

  ro = floor(rp/2) + 1;

  nBeads = length(x);
  nSkipped = 0;
  for i = 1:nBeads
   % add edge checking

     xl = xi(i)-r;
     xh = xi(i)+r;
     yl = yi(i)-r;
     yh = yi(i)+r;
     
     if (xl < 1 || yl < 1 || xh > d1 || yh > d2)
       fprintf('skipping bead %d/%d\n', i,nBeads);
       nSkipped = nSkipped + 1; 
       xo(i) = x(i);
       yo(i) = y(i);
     else
       ccf = mip(xl:xh,yl:yh);
       ccf = fftshift(fftn(ccf.*peakMask));
       ccf = BH_padZeros3d(ccf,'fwd',padVal,'GPU','single');
       ccf = real(ifftn(ifftshift(ccf)));

       [~,maxCoord] = max(ccf(:));
       [mi,mj] = ind2sub(rp,maxCoord);
       mi = (mi- ro(1))./ padBy;
       mj = (mj- ro(2))./ padBy;
       xo(i) = (mi  - xf(i) + x(i));
       yo(i) = (mj  - yf(i) + y(i));
     end
  end

  fprintf('Updated the fit for %d/%d beads\n', nBeads - nSkipped, nBeads);
  
  if (show_results)
    figure, imshow3D(img); hold on
    plot(y,x,'ro','MarkerSize',7); 
    plot(yo,xo, 'b+','MarkerSize', 7);
  end
  
  % We need to add back the 0.5 for imod model coords
  output_pos(lCoord,2) = gather(xo + 0.0);
  output_pos(lCoord,3) = gather(yo + 0.0);
  output_pos(lCoord,[1,4]) = gather(input_pos(lCoord,[1,4]));

end

f = fopen(sprintf('%s.txt',output_name),'w');
fprintf(f,'%d %f %f %f\n',output_pos');
fclose(f);


system(sprintf('point2model -circle %d -zero -image %s %s.txt %s',sampling,img_name,output_name,output_name));

end
