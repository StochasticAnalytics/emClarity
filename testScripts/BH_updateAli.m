% % % % % function [ subTomoMeta ] = BH_updateAli(subTomoMeta,cycle,mapBackIter, pixelSize,...
% % % % %                                         samplingRate,particleRadius)

mapBackIter = 1;
pixelSize = 2.17;
samplingRate = 3;
particleRadius = [160];
cycle=5
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

delete(gcp('nocreate'));
parpool(nTiltSeries);
parResults = cell(nTiltSeries,1);

  try
    geometry = subTomoMeta.(cycleNumber).RawAlign;
    geometry = 'RawAlign';
    fprintf('Using Alignment geometry %s\n',cycleNumber);
  catch
    fprintf('Using Average geometry %s\n',cycleNumber);
    geometry = 'Avg_geometry';
  end
  
for i = 1:nTiltSeries
  parResults{i} = struct();
  parResults{i}.('RawAlign') = subTomoMeta.(cycleNumber).(geometry);
  parResults{i}.('fName') = fieldnames(subTomoMeta.mapBackGeometry.tomoName);
  parResults{i}.('tomoName') = subTomoMeta.mapBackGeometry.tomoName;
  parResults{i}.('coords') = subTomoMeta.mapBackGeometry.(tiltNameList{i}).coords;
end

  
parfor iTilt = 1:nTiltSeries
  
  fprintf('Estimating subtomogram shifts after tomoCPR for tilt %d/ %d\n',iTilt,nTiltSeries);
  tiltName = tiltNameList{iTilt};  
  fullName = sprintf('%s_ali%d_ctf',tiltName,mapBackIter+1);
  
  [~, reconThickness]=system(sprintf('awk ''{if(/^THICKNESS/) print $2}'' %s_ali%d_1_reMod.sh', tiltName, mapBackIter+1));
  reconThickness = str2double(reconThickness);
  [~, reconTiltFile]=system(sprintf('awk ''{if(/^TILTFILE/) print $2}'' %s_ali%d_1_reMod.sh', tiltName, mapBackIter+1));
  [~,baseName,extension]=fileparts(reconTiltFile);
  [~, reconLocalFile]=system(sprintf('awk ''{if(/^LOCALFILE/) print $2}'' %s_ali%d_1_reMod.sh', tiltName, mapBackIter+1));
  
  if (isempty(reconLocalFile))
    reconLocalFile = ' ';
  else
    reconLocalFile = sprintf('-LOCALFILE %s ',reconLocalFile);
  end

  
  
     % Resample the synthetic stack - alignments need to be binned first,
     % TODO just do this in mem prior to saving in tomoCPR
    fprintf('resampling the tomoCPR ref to new globalAli\n');
    system(sprintf('awk -v B=%d ''{print $1,$2,$3,$4,$5/B,$6/B}'' %s.tltxf > %s_bin.xf',samplingRate, fullName, fullName));
    system(sprintf('newstack  -xf %s_bin.xf  %s_ali%d_1_mapBack.st %s_resamp.st >> ./.tomoCPR_log.txt', ...
                    fullName,tiltName,mapBackIter+1,fullName)); 
                  
                  
% % %   end
              
              
  % Get size info (TODO get the bin number from tomoCPR)

  % Reconstruct the original and warped synthetic vol
  fprintf('Reconstructing the synthetic tomos positions\n');
  recCmd = sprintf(['tilt -input %s_resamp.st -output %s_resamp.rec ', ...
                    '-TILTFILE %s.tlt -RADIAL 0.5,0.05 -UseGPU 0 ', ...
                    '-THICKNESS %d -COSINTERP 0 -RotateBy90 ', ...
                    '-LOCALFILE %s.local -FakeSIRTiterations 30 >> ./.tomoCPR_log.txt'], ...
                    fullName, fullName, fullName, reconThickness ,fullName);
                 
  system(recCmd);
  
  recCmd = sprintf(['tilt -input %s_ali%d_1_mapBack.st -output %s_orig.rec ', ...
                    '-TILTFILE %s.rawtlt -RADIAL 0.5,0.05 -UseGPU 0 ', ...
                    '-THICKNESS %d -COSINTERP 0 -RotateBy90 ', ...
                    '%s -FakeSIRTiterations 30 >> ./.tomoCPR_log.txt'], ...
                    tiltName,mapBackIter+1, fullName, baseName, reconThickness , reconLocalFile);
                 
  system(recCmd);
  
  

  % Load the tomos into main mem
  fprintf('Loading the reference and target tomos\n');
  tomoPre = getVolume(MRCImage(sprintf('%s_orig.rec',fullName)));
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
  
  

iSubRegion = 1;
for iY = 1:nY
  for iX = 1:nX
      iSubRegion
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
        iImg = BH_resample3d(iRef,[0,0,0], peakShift ,'Bah','GPU','forward');    
        figure, imshow3D(gather(iImg));

      end

    totalShift;
    shiftsOUT(iSubRegion,3:5) = totalShift; % (samplingRate*tmpTomoBin).*totalShift;
    iSubRegion = iSubRegion + 1;
  end

end


boxSize = gather(ceil(3.*max(particleRadius)./(pixelSize.*samplingRate)).*[1,1,1]);
iMask = gather(BH_mask3d('sphere',boxSize, boxSize.*0.33, [0,0,0]));

  n = 0;

  for iTomo = 1:length(parResults{iTilt}.fName)
    randName = randi([1000,9999],1);
    if strcmp(parResults{iTilt}.tomoName.(parResults{iTilt}.fName{iTomo}).tiltName, tiltNameList{iTilt})
      
      writeOutPyAli(randName);
      
      tomoNum = parResults{iTilt}.tomoName.(parResults{iTilt}.fName{iTomo}).tomoNumber;
      rCoords = parResults{iTilt}.coords(tomoNum,:);
      tomoReconGeom = BH_offsets(rCoords, sprintf('%s_ali%d_1_mapBack.st',tiltName,mapBackIter+1), samplingRate);
      
      originRealTomo = floor(tomoReconGeom(1,:)./2) + 1;
      originFakeTomo = floor([d1,d2,d3]./2) + 1;
      reconShift = tomoReconGeom(2,1:3);

      lowerLeftVol = originFakeTomo + reconShift - originRealTomo + 1;

      nSubTomos = size(parResults{iTilt}.RawAlign.(parResults{iTilt}.fName{iTomo}),1);
        for iSubTomo = 1:nSubTomos
  
          minDist = shiftsOUT(:,1:2) - parResults{iTilt}.RawAlign.(parResults{iTilt}.fName{iTomo})(iSubTomo,11:12);
          minDist = sqrt(sum(minDist.*minDist,2));
          [~,minCoord] = min(minDist);
          
          center =  parResults{iTilt}.RawAlign.(parResults{iTilt}.fName{iTomo})(iSubTomo,11:13) ./ samplingRate + lowerLeftVol;
          angles =  reshape(parResults{iTilt}.RawAlign.(parResults{iTilt}.fName{iTomo})(iSubTomo,17:25),3,3) ;
          
 
                % get the reference, i.e. the model in the original position
          [ indVAL, refPadVAL, refShiftVAL ] = ...
                              BH_isWindowValid([d1,d2,d3], ...
                                                boxSize, boxSize.*0.3, ...
                                                center);

          iRef = (tomoPre(indVAL(1,1):indVAL(2,1), ...
                         indVAL(1,2):indVAL(2,2), ...
                         indVAL(1,3):indVAL(2,3)));

          if any(refPadVAL(:))
            [ iRef ] = BH_padZeros3d(iRef,  refPadVAL(1,1:3), ...
                                            refPadVAL(2,1:3), 'cpu', 'single');
          end
                                              
          [ indVAL, imgPadVAL, imgShiftVAL ] = ...
                              BH_isWindowValid([d1,d2,d3], ...
                                                boxSize, boxSize.*0.3, ...
                                                center+gather(shiftsOUT(minCoord,3:5)));
                                              
           
          iImg =  (tomoPost(indVAL(1,1):indVAL(2,1), ...
                               indVAL(1,2):indVAL(2,2), ...
                               indVAL(1,3):indVAL(2,3)));      

          if any(imgPadVAL(:))
            [ iImg ] = BH_padZeros3d(iImg,  imgPadVAL(1,1:3), ...
                                            imgPadVAL(2,1:3), 'cpu', 'single');
          end  
          

          % For now, we'll assum there is a writable /tmp dir and that the
          % python script to run chimera already exists in the FSC dir.

          
          refName = sprintf('/tmp/ref%d.mrc', randName);
          imgName = sprintf('/tmp/img%d.mrc', randName);
          SAVE_IMG(iRef.*iMask, refName, pixelSize.*samplingRate);
          SAVE_IMG(iImg.*iMask, imgName, pixelSize.*samplingRate);
          
          
        system(sprintf('chimera --nogui --script "/tmp/fitInMap_%d.py %s %s /tmp/fitInMap_%d.txt > /dev/null" ',...
                      randName, refName,imgName,randName));

      % Read in the results
      % 1:9 rotation matrix 10:12 = dXYZ
      bestAnglesAndShifts = load(sprintf('/tmp/fitInMap_%d.txt',randName));
      % The results from fit in map are the transpose of Bah, forward rotmat
      bestAngles = reshape(bestAnglesAndShifts([1,4,7,2,5,8,3,6,9]),3,3);
      
      % The results are in Angstrom, convert to unbinned pixels
      bestShifts = bestAnglesAndShifts(10:12)'./pixelSize;

%       resampleCheck = (BH_resample3d(iRef, ...
%                                   bestAngles, ...
%                                   1.*bestShifts'./samplingRate, ...
%                                   'Bah', 'GPU', 'forward'));
%                                 
%                 SAVE_IMG(resampleCheck, '/tmp/check.mrc', pixelSize.*samplingRate);
%                 error('sdf');

%       system(sprintf(' rm %s %s /tmp/fitInMap_%d.txt',refName,imgName,randName));%s %s 
      

      iShift = (gather(shiftsOUT(minCoord,3:5)).*samplingRate + bestShifts);
      
      parResults{iTilt}.RawAlign.(parResults{iTilt}.fName{iTomo})(iSubTomo,11:13) = ...
                          parResults{iTilt}.RawAlign.(parResults{iTilt}.fName{iTomo})(iSubTomo,11:13) + ...
                          gather(iShift);   
                        
        end
    end % if clause on tomo in this tilt

  system(sprintf(' rm /tmp/fitInMap_%d.py', randName));
  
  end % end of loop over tomos
end % end loop over tilts


% Gather up the results
for iTilt = 1:nTiltSeries
  for iTomo = 1:length(parResults{iTilt}.fName)
     if strcmp(parResults{iTilt}.tomoName.(parResults{iTilt}.fName{iTomo}).tiltName, tiltNameList{iTilt})
       subTomoMeta.(cycleNumber).(geometry).(parResults{iTilt}.fName{iTomo}) = ...
                          parResults{iTilt}.RawAlign.(parResults{iTilt}.fName{iTomo});
                            
       
     end
    
  end
end
% % % % % end

function writeOutPyAli(randName)

fOUT = fopen(sprintf('/tmp/fitInMap_%d.py', randName),'w');

fprintf(fOUT,['\n'...
'from sys import argv\n\n',...
'def fit_map_in_map(map1_path, map2_path, xformName,\n',...
'                   initial_map1_transform = None,\n',...
'                   map1_threshold = 1.5,\n',...
'                   ijk_step_size_min = 0.01,\n',...    
'                   ijk_step_size_max = 1.5,\n',...     
'                   max_steps = 5000,\n',...
'                   optimize_translation = True,\n',...
'                   optimize_rotation = False):\n',...
'  from VolumeViewer import open_volume_file\n',...
'  map1 = open_volume_file(map1_path)[0]\n',... 
'  map2 = open_volume_file(map2_path)[0]\n\n',...
'  if initial_map1_transform:\n',...
'    from Matrix import chimera_xform\n',...
'    xf = chimera_xform(initial_map1_transform)\n',...
'    map1.surface_model().openState.globalXform(xf)\n\n',...    
'  use_threshold = (map1_threshold != None)\n\n',...  
'  from FitMap.fitmap import map_points_and_weights, motion_to_maximum\n',...
'  points, point_weights = map_points_and_weights(map1, use_threshold)\n\n',...
'  move_tf, stats = motion_to_maximum(points, point_weights, map2, max_steps,\n',...
'                                     ijk_step_size_min, ijk_step_size_max,\n',...
'                                     optimize_translation, optimize_rotation)\n\n',...
'  import Matrix\n',...
'  if initial_map1_transform:\n',...
'    move_tf = Matrix.multiply_matrices(move_tf, initial_map1_transform)\n\n',...
'  f = open(xformName,''w'')\n',...
'  for i in range(4):\n',...
'    for j in range(3):\n',...
'      f.write(''{:f}''.format(move_tf[j][i]))\n',...
'      f.write("\\n")\n\n',...
'  f.close()\n',...
'  tfs = Matrix.transformation_description(move_tf)\n',...
't = fit_map_in_map(argv[1],argv[2],argv[3])\n']);

fclose(fOUT);


end
