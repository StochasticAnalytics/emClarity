% Test fourierTransformer class

testSize = { [257,129,65], [256,128,64]} ;

testFFT = false;
testInterp3d = true;

if (testFFT)
  for iTest = 1:2
    testImg = randn(testSize{iTest}, 'single', 'gpuArray');

    BH_out = real(ifftn(fftn(testImg)));

    f = fourierTransformer(testImg);
    MEX_out = real(f.invFFT(f.fwdFFT(testImg,2,0)));

    diff = BH_returnImage2ImageRMSD(BH_out, MEX_out);
    
    fprintf("For fft test on size %d, the mean RMSD is %3.3e\n", ...
            iTest, diff);
  end
end

if (testInterp3d)
  for iTest = 1:2
    
    testImg = rand(testSize{iTest}, 'single', 'gpuArray');
    shift1 = [30,-50,0];
    rotMat = [10,-20,-30];
    direction = {'forward','inv'};
    firstRun = true;
    
    for iDir = 1:2
      for iShift = [-1,1]
        for iRot = [-1,1]
          
          dir = direction{iDir};
          BH_out = BH_resample3d(testImg,iRot.*rotMat, ...
                                         iShift.*shift1,'Bah','GPU',dir);
    
     if (firstRun)
        [resampler, MEX_out] = interpolator(testImg, iRot.*rotMat, iShift.*shift1, 'Bah', dir);
        firstRun = false;
     else
      MEX_out = resampler.interp3d(testImg, iRot.*rotMat, iShift.*shift1, 'Bah', dir);
     end
     
    diff = BH_returnImage2ImageRMSD(BH_out, MEX_out);
    
    fprintf("For interp test on size %d rot shift dir %d %d %s, the mean RMSD is %3.3e\n", ...
            iTest, iRot, iShift, dir, diff);
        end
      end
    end
    
    
  
  
  end
end



