function [ ] = BH_multi_loadOrCalcWeight(masterTM, ...
                                                            ctfGroupList, ...
                                                            tomoName, ...
                                                            samplingRate, ...
                                                            sizeCalc,...
                                                            geometry_tmp,...
                                                            cutPrecision,...
                                                            gpuIDX)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    iTiltName = masterTM.mapBackGeometry.tomoName.(tomoName).tiltName;
    wgtName = sprintf('cache/%s_bin%d.wgt',iTiltName,samplingRate); 
    
    try
      calcWeight = 0;
      iTiltMRCobj = MRCImage(wgtName,0);
      iTiltHeader = getHeader(iTiltMRCobj);
      
      % Explicitly assuming the weights are cubic, st the z dimension of
      % the montage corresponds to each weight dimension
      if any(sizeCalc - iTiltHeader.nZ)
        calcWeight = 1;
        fprintf('the current weight with z = %d is not equal %d %d %d\n',...
                 iTiltHeader.NZ,sizeCalc);
      else
        fprintf('the wgt %s exists and will be used as is\n',wgtName);
      end
    catch
      calcWeight = 1;
      gDev = gpuDevice(gpuIDX);
   
   
      nCtfGroups = ctfGroupList.(tomoName)(1);

      if ~(exist('./cache','dir'))
        system('mkdir cache');
      end
      fprintf('Reweighting.\n\n')
      % Make a wedge mask that can be interpolated with no extrapolation for
      % calculating wedge weighting in class average alignment. 
       kVal = 0;
      [ OUTPUT ] = BH_multi_iterator( [sizeCalc;kVal.*[1,1,1]], 'extrapolate' );
      
      if isempty(OUTPUT)
        fprintf('effing OUTPUT is empty\n');
        OUTPUT = sizeCalc;
      else
        fprintf('OUTPUT %f %f %f\n', OUTPUT(1,:));
      end
    end
                                            


        if (calcWeight)
          wedgeMask = [];

          
           fprintf('sizeCalc %d %d %d\n sizeWgt %d %d %d\n', sizeCalc, size(wedgeMask));
           % Keep tomo name as cell, so switching to explicit calc makes
           % sense to me later.
          [ wedgeMask ] = BH_weightMask_dp(masterTM, OUTPUT(1,:), samplingRate,...
                                           {{tomoName},geometry_tmp}, ...
                                           cutPrecision, 'GPU');           

          
          for iCtfGroup = 1:nCtfGroups
            wedgeMask{iCtfGroup} = single(wedgeMask{iCtfGroup});
          end
          [montOUT, ~] = BH_montage4d(wedgeMask, '');
          SAVE_IMG(MRCImage(montOUT),wgtName);
          montOUT = [];
           for iCtfGroup = 1:nCtfGroups
            wedgeMask{iCtfGroup} = [];
          end         
        end  
        
end

