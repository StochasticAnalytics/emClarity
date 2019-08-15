function [ subTomoMeta ] = BH_updateAli(subTomoMeta,cycle,mapBackIter, pixelSize,...
                                        samplingRate,particleRadius)

% % % mapBackIter =0;
% % % pixelSize = 2.17;
% % % samplingRate = 4;
% % % particleRadius = [160,160,160];

% Build a list of very large local shifts

% % % % % maskRadius = max(particleRadius)./(pixelSize*samplingRate);
% % % % % boxSize = floor(maskRadius*4); % Is a bigger box better?
% % % % % mask = BH_mask3d('sphere',boxSize.*[1,1,1],maskRadius.*[2,2,2],[0,0,0]);

% % % % % [U,V,W] = BH_multi_gridCoordinates(boxSize.*[1,1,1],'Cartesian','GPU',{'none'},1,0,0);
% % % % % U = U .* (2i*pi);
% % % % % V = V .* (2i*pi);
% % % % % W = W .* (2i*pi);

% % % % % tmpTomoBin = ceil(6/(pixelSize*samplingRate));

checkResults = false;
% % % doMedianFilter = true;
% % % 
% % % load('../cycle003_releaseTest_backup.mat')
% % % cycle = 2;% subTomoMeta.currentCycle;
cycleNumber = sprintf('cycle%0.3u',cycle);

tiltNameList = fieldnames(subTomoMeta.mapBackGeometry);
tiltNameList = tiltNameList(~ismember(tiltNameList,{'tomoName','viewGroups'}));
nTiltSeries = length(tiltNameList);
 % TODO run this in parallel
for iTilt = 1:nTiltSeries
  
  fprintf('Estimating subtomogram shifts after tomoCPR for tilt %d/ %d\n',iTilt,nTiltSeries);
  tiltName = tiltNameList{iTilt};  
  fullName = sprintf('%s_ali%d_ctf',tiltName,mapBackIter+1);
  
% % % 
% % %   if (doMedianFilter)
% % %     fprintf('Applying median filter to stack\n'); %#ok<*UNRCH>
% % %     medFiltStack = getVolume(MRCImage(sprintf('%s_ali%d_shifted.st',tiltName,mapBackIter+1)));
% % %     for iPrj = 1:size(medFiltStack,3)
% % %       medFiltStack(:,:,iPrj) = medfilt2(medFiltStack(:,:,iPrj),[3,3]);
% % %     end
% % %     SAVE_IMG(MRCImage(medFiltStack),sprintf('%s_MedFilt.st',fullName),pixelSize.*samplingRate);
% % %     medFiltStack = [];
% % % 
% % %     % Resample the synthetic stack - alignments need to be binned first
% % %     fprintf('resampling the tomoCPR ref to new globalAli\n');
% % %     system(sprintf('awk -v B=%d ''{print $1,$2,$3,$4,$5/B,$6/B}'' %s.tltxf > %s_bin.xf',tmpTomoBin*samplingRate, fullName, fullName));
% % %     system(sprintf('newstack -meansd 0,3 -float 2 -xf %s_bin.xf -bin %d %s_MedFilt.st %s_resamp.st >> ./.tomoCPR_log.txt', ...
% % %                     fullName,tmpTomoBin,fullName,fullName));
% % %   else
    
     % Resample the synthetic stack - alignments need to be binned first,
     % TODO just do this in mem prior to saving in tomoCPR
    fprintf('resampling the tomoCPR ref to new globalAli\n');
    system(sprintf('awk -v B=%d ''{print $1,$2,$3,$4,$5/B,$6/B}'' %s.tltxf > %s_bin.xf',samplingRate, fullName, fullName));
    system(sprintf('newstack  -xf %s_bin.xf  %s_ali%d_1_mapBack.st %s_resamp.st >> ./.tomoCPR_log.txt', ...
                    fullName,tiltName,mapBackIter+1,fullName)); 
                  
                  
% % %   end
              
              
  % Get size info (TODO get the bin number from tomoCPR)
  iMrcObj = MRCImage(sprintf('%s_ali%d.tmpTomo1',tiltName,mapBackIter + 1));
  iHeader = getHeader(iMrcObj);

  % Reconstruct the warped synthetic vol
  fprintf('Reconstructing the new positions\n');
  recCmd = sprintf(['tilt -input %s_resamp.st -output %s_resamp.rec ', ...
                    '-TILTFILE %s.tlt -RADIAL 0.5,0.05 -UseGPU 0 ', ...
                    '-THICKNESS %d -COSINTERP 0 -RotateBy90 ', ...
                    '-LOCALFILE %s.local -FakeSIRTiterations 30 >> ./.tomoCPR_log.txt'], ...
                    fullName, fullName, fullName, iHeader.nZ ,fullName);
                 
  system(recCmd);
  
  % Load the tomos into main mem
  fprintf('Loading the reference and target tomos\n');
  tomoPre = getVolume(MRCImage(sprintf('%s_ali%d.tmpTomo1',tiltName,mapBackIter+1)));
  tomoPost= getVolume(MRCImage(sprintf('%s_resamp.rec',fullName)));

  if (size(tomoPre) ~= size(tomoPost))
    error('tomoPre and tomoPost sizes do not match!\n');
  else
    [d1,d2,d3] = size(tomoPre);
  end
  
  % Get the local areas
  system(sprintf('head -n 1 %s.local > %s_header.txt',fullName,fullName));
  localHeader = load(sprintf('%s_header.txt',fullName));
  nX = localHeader(1);
  nY = localHeader(2);
  oX = floor(localHeader(3)./samplingRate);
  oY = floor(localHeader(4)./samplingRate);
  oZ = floor(d3/2) + 1;
  incX= floor(localHeader(5)./samplingRate);
  incY= floor(localHeader(6)./samplingRate);
  nLocalAreas = nX*nY;
  
  boxSize = floor([oX,oY,oZ].*1.9);
  shiftsOUT = zeros(5,nLocalAreas,'gpuArray');
  
  

% % % % %   % % Get the estimated shifts due to the global alignment
% % % % %   % % tomoNumber, subTomoIDX, xOUT,yOUT,zOUT,xIN,yIN,zIN ( Some positions
% % % % %   % may be entirely ignored
% % % % %   awkCmd=sprintf('awk ''FNR==NR{x[$6]=$2; y[$6]=$3; z[$6]=$4; next}{if(x[$1]) {print $2,$3,x[$1],y[$1],z[$1],$4,$5,$6} else {print $2,$3,-9999,-9999,-9999,$4,$5,$6}}'' ');
% % % % %   awkCmd=sprintf('%s %s.xyz %s_ali%d.coord_start  > %s.globalDeltaXYZ', awkCmd,fullName, tiltName, mapBackIter+1, fullName); 
% % % % %   system(awkCmd);
% % % % %   dXYZ = load(sprintf('%s.globalDeltaXYZ',fullName));
% % % % %   nonSampledMask = (dXYZ(:,3) == -9999 & dXYZ(:,4) == -9999 & dXYZ(:,5) == -9999);
% % % % %   dXYZ(~nonSampledMask,5) = dXYZ(~nonSampledMask,5) + (floor(((d3+samplingRate)*samplingRate)/2)+1);
% % % % %   dXYZ(nonSampledMask,3:5) = dXYZ(nonSampledMask,6:8);
% % % % %   % Z is in centered coordinates, and also inverted due to IMOD recon geom
% % % % %   dXYZ(:,3:8) = dXYZ(:,3:8) ./ (tmpTomoBin * samplingRate);
% % % % % 
% % % % %   % dXYZ(:,[5,8]) = (d3 - (dXYZ(:,[5,8]) + (floor(d3/2)+1)));
% % % % %   nSubTomos = size(dXYZ,1);
% % % % %   shiftsOUT = zeros(nSubTomos,5,'gpuArray');
% % % % % 
% % % % %   shiftsOUT(:,1:2) = dXYZ(:,1:2);
% % % % % 
% % % % %   for iSubTomo = 1:nSubTomos
% % % % % 
% % % % % 
% % % % %     img_XYZ = dXYZ(iSubTomo,3:5);
% % % % %     ref_XYZ = dXYZ(iSubTomo,6:8);
    
% % % % %     
% % % % %     if (abs(sum(img_XYZ - ref_XYZ))>1e-6)
iSubRegion = 1;
for iY = 1:nY
  for iX = 1:nX
    
      rXYZ = [oX + incX*(iX-1), oY + incY*(iY-1), oZ];
      shiftsOUT(iSubRegion,1:2) = rXYZ(1:2) .* samplingRate;
     
 
      % get the reference, i.e. the model in the original position
      [ indVAL, refPadVAL, refShiftVAL ] = ...
                              BH_isWindowValid([d1,d2,d3], ...
                                                boxSize, boxSize.*0.3, ...
                                                rXYZ);

      if ~isnumeric(indVAL)
        error('out of bounds');
      end
      
      iRef =  gpuArray(tomoPre(indVAL(1,1):indVAL(2,1), ...
                                indVAL(1,2):indVAL(2,2), ...
                                indVAL(1,3):indVAL(2,3)));

      if any(refPadVAL(:))
        [ iRef ] = BH_padZeros3d(iRef,  refPadVAL(1,1:3), ...
                                        refPadVAL(2,1:3), 'GPU', 'single');
      end

      % get the target img, i.e. the model in the shifted position
      [ indVAL, imgPadVAL, imgShiftVAL ] = ...
                              BH_isWindowValid([d1,d2,d3], ...
                                                boxSize, boxSize.*0.3, ...
                                                rXYZ);
      if ~isnumeric(indVAL)
        error('out of bounds');
      end

      iImg =  gpuArray(tomoPost(indVAL(1,1):indVAL(2,1), ...
                               indVAL(1,2):indVAL(2,2), ...
                               indVAL(1,3):indVAL(2,3)));      

      if any(imgPadVAL(:))
        [ iImg ] = BH_padZeros3d(iImg,  imgPadVAL(1,1:3), ...
                                        imgPadVAL(2,1:3), 'GPU', 'single');
      end  


      % TODO, apply this shifts as phase shifts during the cross correlation.
      if (checkResults)
        figure, imshow3D(gather(iRef));
        figure, imshow3D(gather(iImg));
        SAVE_IMG(MRCImage(gather(iRef)),'tmpRef.mrc');
        SAVE_IMG(MRCImage(gather(iImg)),'tmpImg.mrc');
      end


      % The boxing out and subsequent sub pixel shifts give the estimated
      % displacement (imgXYZ-refXYZ). Any residual error is in the peakShift
      [peakShift] = BH_multi_xcf_Translational(fftn(iImg),...
                                          conj(fftn(iRef)),'',[3,3,3]);
% % % %       [peakShift] = BH_multi_xcf_Translational(fftn(iImg.*mask).*exp(U.*imgShiftVAL(1)+V.*imgShiftVAL(2)+W.*imgShiftVAL(3)),...
% % % %                                           conj(fftn(iRef.*mask).*exp(U.*refShiftVAL(1)+V.*refShiftVAL(2)+W.*refShiftVAL(3))),'',[3,3,3]);


      totalShift = peakShift ;
      if (checkResults)
        peakShift
% 
%         refShiftVAL - imgShiftVAL
%         ref_XYZ - img_XYZ
%         totalShift = peakShift + (imgShiftVAL - refShiftVAL) + ( (img_XYZ - ref_XYZ))
        iImg = BH_resample3d(iRef,[0,0,0], peakShift ,'Bah','GPU','forward');    
        figure, imshow3D(gather(iImg));
        error('sdf')

      end
% % % % %     else
% % % % %       totalShift = [0,0,0];
% % % % %     end
    totalShift
    shiftsOUT(iSubRegion,3:5) = totalShift; % (samplingRate*tmpTomoBin).*totalShift;
    iSubRegion = iSubRegion + 1;
  end

end



  fName = fieldnames(subTomoMeta.mapBackGeometry.tomoName);
  n = 0;
    for iTomo = 1:length(fName)
      fName{iTomo};
      tomoNum = subTomoMeta.mapBackGeometry.tomoName.(fName{iTomo}).tomoNumber;
      nSubTomos = size(subTomoMeta.(cycleNumber).RawAlign.(fName{iTomo}),1);
        for iSubTomo = 1:nSubTomos
  
          subTomoMeta.(cycleNumber).RawAlign.(fName{iTomo})(iSubTomo,11:12);
          minDist = shiftsOUT(:,1:2) - subTomoMeta.(cycleNumber).RawAlign.(fName{iTomo})(iSubTomo,11:12);
          minDist = sqrt(sum(minDist.*minDist,2));
          [~,minCoord] = min(minDist);
              subTomoMeta.(cycleNumber).RawAlign.(fName{iTomo})(iSubTomo,11:13) = ...
                          subTomoMeta.(cycleNumber).RawAlign.(fName{iTomo})(iSubTomo,11:13) + ...
                          gather(shiftsOUT(minCoord,3:5));   
                        
                        shiftsOUT(minCoord,3:5);
        end
    end

% % % %   for iTomo = 1:length(fName)
% % % %     fName{iTomo}
% % % %     if strcmp(subTomoMeta.mapBackGeometry.tomoName.(fName{iTomo}).tiltName, tiltName)
% % % %       tomoNum = subTomoMeta.mapBackGeometry.tomoName.(fName{iTomo}).tomoNumber
% % % %         for iSubTomo = 1:nSubTomos
% % % %           if (shiftsOUT(iSubTomo,1) == tomoNum)
% % % %             idx=find(subTomoMeta.(cycleNumber).RawAlign.(fName{iTomo})(:,4) == shiftsOUT(iSubTomo,2),1,'first');
% % % %             idx
% % % %             if (shiftsOUT(iSubTomo,3) == -9999)
% % % %               subTomoMeta.(cycleNumber).RawAlign.(fName{iTomo})(idx,26) = -9999; % TODO add a counter   
% % % %               n = n + 1;
% % % %             else
% % % %               subTomoMeta.(cycleNumber).RawAlign.(fName{iTomo})(idx,11:13) = ...
% % % %                           subTomoMeta.(cycleNumber).RawAlign.(fName{iTomo})(idx,11:13) + ...
% % % %                           gather(shiftsOUT(iSubTomo,3:5));
% % % %             end
% % % %           end
% % % %         end
% % % %     end
% % % %   end

  % Remove the binned tomo
%   system(sprintf('rm %s_resamp.rec',fullName));
end % end loop over tilts

end
