function  BH_geometryAnalysis( PARAMETER_FILE, CYCLE, ...
                               STAGEofALIGNMENT, OPERATION, ...
                               VECTOR_OP, HALF_SET)
%Check the alignment at various stages.
%   Detailed explanation goes here



if (nargin ~= 6)
  error('args = (PARAMETER_FILE, CYCLE, STAGEofALIGNMENT, OPERATION, REMOVE_CLASS, HALF_SET)')
end

CYCLE = EMC_str2double(CYCLE);
cycleNumber = sprintf('cycle%0.3u', CYCLE);
undoOP = 0;
listTomos = 0;
if strcmpi(VECTOR_OP, 'undo')
  undoOP = 1;
  fprintf('Undoing %s action', OPERATION);
elseif (strcmpi(VECTOR_OP, 'listTomos'))
    listTomos = 1;
    fprintf('Listing tomos, delete those that should be removed, and run again replacing "listTomos" with "tomoList.txt"\n');
else
  try
    % should be a text file of only x,y,z (model2points imodModel classes.txt)
    % or
    % a list of names of tomograms to remove
    if ~(strcmp(VECTOR_OP,'tomoList.txt'))
      [~,modNAME,~] = fileparts(VECTOR_OP);
      system(sprintf('model2point %s %s.txt > /dev/null',VECTOR_OP,modNAME));

      VECTOR_OP = importdata(sprintf('%s.txt',modNAME));
    end

  catch
    % or a 3-vector containing shifts to apply to all volumes, or range for
    % randomization of 3 euler angles

      VECTOR_OP = EMC_str2double(VECTOR_OP);
      if size(VECTOR_OP) == [1,3]
        shiftXYZ = VECTOR_OP'
      elseif size(VECTOR_OP) == [3,1]
        shiftXYZ = VECTOR_OP
      elseif ~any(size(VECTOR_OP))
        % size of EMC_str2double('undo') = [0 0]

      else
        error('shift values must be a three vector or string "undo"\n');
      end
  end
end
  
startTime =  datetime("now");

cycleNumber = sprintf('cycle%0.3u', CYCLE);

emc = BH_parseParameterFile(PARAMETER_FILE);
try
  conserveDiskSpace = emc.('conserveDiskSpace');
catch
  conserveDiskSpace = 0;
end
try
  percentCut   = emc.('removeBottomPercent');
catch
  percentCut = 0.0;
end


if percentCut < 0 
  error('removeBottomPercent must be a fraction between 0 1, or the top N volumes to keep');
end


halfSet = HALF_SET
if strcmpi(halfSet, 'ODD')
  halfNUM = 1;
elseif strcmpi(halfSet, 'EVE')
  halfNUM = 2;
elseif strcmpi(halfSet, 'STD')
  halfNUM = [1,2];
else
  error('HALF_SET should be EVE or ODD or STD')
end



switch STAGEofALIGNMENT
  case 'TiltAlignment'
    fieldPrefix = 'Ref'
  case 'RawAlignment'
    fieldPrefix = 'Ref';
  case 'Cluster_cls'
    fieldPrefix = 'Cls';
    STAGEofALIGNMENT = 'Cluster';
  otherwise
    error('STAGEofALIGNMENT incorrect, should be Cluster_cls or [Tilt,Raw]Alignment, not %s', ...
                                                        STAGEofALIGNMENT);
end


samplingRate = emc.(sprintf('%s_samplingRate','Ali'));

className    = emc.(sprintf('%s_className',fieldPrefix));


outputPrefix = sprintf('%s_%s', cycleNumber, emc.('subTomoMeta'));


%try
  load(sprintf('%s.mat', emc.('subTomoMeta')), 'subTomoMeta');
  switch STAGEofALIGNMENT
    case 'TiltAlignment'
     geometry = subTomoMeta.tiltGeometry;
     subTomoMeta.(cycleNumber).(sprintf('Pre_%s_tiltGeometry', OPERATION)) = geometry;
    case 'RawAlignment'
      if (undoOP)
        subTomoMeta.(cycleNumber).RawAlign = ...
        subTomoMeta.(cycleNumber).(sprintf('Pre_%s_RawAlign', OPERATION));
        save(emc.('subTomoMeta'), 'subTomoMeta');
        error('No Error, just exiting.\n')
      else      
        geometry = subTomoMeta.(cycleNumber).RawAlign;
        subTomoMeta.(cycleNumber).(sprintf('Pre_%s_RawAlign', OPERATION)) = geometry;
      end
    case 'Cluster'
  
        try
          classVector{1}  = emc.(sprintf('%s_classes_odd',fieldPrefix));
        catch
          classVector{1}  = emc.(sprintf('%s_classes',fieldPrefix));
        end
        
        classVector{2}  = emc.(sprintf('%s_classes_eve',fieldPrefix));

        try
          classCoeffs{1} =  emc.('Pca_coeffs_odd');
          classCoeffs{2} =  emc.('Pca_coeffs_eve');
        catch
          classCoeffs{1} = emc.('Pca_coeffs');
        end

      cN = sprintf('%s_%d_%d_nClass_%d_%s',outputPrefix,classCoeffs{halfNUM(1)}(1,1), ...
                                    classCoeffs{halfNUM(1)}(1,end), className, halfSet)
      imgNAME = sprintf('class_%d_Locations_%s_%s_NoWgt', className, fieldPrefix, halfSet);
      imgNAMEStd = sprintf('class_%d_Locations_%s_%s_NoWgt', className, fieldPrefix, 'ODD');
     
      
      if (undoOP)
        subTomoMeta.(cycleNumber).ClusterResults.(cN) = ...
        subTomoMeta.(cycleNumber).(sprintf('Pre_%s_ClusterResults', OPERATION)).(cN);
        save(emc.('subTomoMeta'), 'subTomoMeta');
        error('No Error, just exiting.\n')
      else     

        if strcmpi(fieldPrefix, 'Cls')
         geometry = subTomoMeta.(cycleNumber).('ClusterClsGeom');
         clusterGeom = 'ClusterClsGeom';
        elseif strcmpi(fieldPrefix, 'Ref')
         geometry = subTomoMeta.(cycleNumber).('ClusterRefGeom');
         clusterGeom = 'ClusterRefGeom';
        end
        
        %geometry = subTomoMeta.(cycleNumber).ClusterResults.(cN);
        subTomoMeta.(cycleNumber).(sprintf('Pre_%s_ClusterResults', OPERATION)).(cN) = geometry;
        try
          locations= subTomoMeta.(cycleNumber).(imgNAME){2};
        catch
          locations= subTomoMeta.(cycleNumber).(imgNAMEStd){2};
        end
      end
    otherwise
      error('STAGE_ALIGNMENTS: [Tilt,Class,No,Raw]Alignment, not %s', ...
                                                            STAGEofALIGNMENT);
  end
% % %     pathList= subTomoMeta.mapPath;
% % %     extList = subTomoMeta.mapExt;
mapBackIter = subTomoMeta.currentTomoCPR; 

    masterTM = subTomoMeta; clear subTomoMeta
%   
% catch 
%   error('failed to load geometry')
% end

tomoList = fieldnames(geometry);
nTomograms = length(tomoList);


switch OPERATION
  case 'SwitchCurrentCycle'
    masterTM.currentCycle = VECTOR_OP(1);
  case 'SwitchCurrentTomoCpr'
    masterTM.currentTomoCPR = VECTOR_OP(1);    
  case 'SwitchExposureDose'
    % While transitioning , edit the dose column in the tilt geometry.
    
    for iTomo = 1:nTomograms
      TLT = geometry.(tomoList{iTomo});
      nPrjs = size(TLT,1);
      % Assuming a bi-directional tilt scheme with smallest abs value as first tilt,
      % which could be wrong
      dosePerTilt = [TLT(:,1), TLT(:,4), zeros(nPrjs,1)];
      % Make sure arranged from negative to positive
      dosePerTilt = sortrows(dosePerTilt,2);
      [~,firstTilt] = min(abs(dosePerTilt(:,2)));
      dosePerTilt(1:firstTilt,:) = sortrows(dosePerTilt(1:firstTilt,:),-2);
      exposure = VECTOR_OP(1)./nPrjs
      for iExposure = 1:nPrjs
        dosePerTilt(iExposure,3) = iExposure .* exposure;
      end
      dosePerTilt
      for iPrj = 1:nPrjs
        CUMeDOSE = dosePerTilt(find(dosePerTilt(:,1) == TLT(iPrj,1)),3);
        geometry.(tomoList{iTomo})(iPrj,11) = CUMeDOSE
      end
    end
  case 'UpdateTilts'
 
    
    % Usually the tilt angle change is pretty small, so re-calculating the 
    % weights may not be strictly necessary -- consider an option to not 
    % clear them.
%     system('rm cache/*.wgt');
    % This is most important to update the changes to the tilt angles in
    % the TLT geometry which affect future tomoCPR.
    if ( conserveDiskSpace )
      system('rm cache/*.wgt*');
      system('rm cache/*.rec*');
    end
    for iTomo = 1:nTomograms
      STACK_PRFX = ...
               masterTM.mapBackGeometry.tomoName.(tomoList{iTomo}).tiltName;

      

      newTLT = sprintf('fixedStacks/ctf/%s_ali%d_ctf.tlt', ...
                                                  STACK_PRFX,mapBackIter+1); 

    
      geometry.(tomoList{iTomo}) = load(newTLT);
      fprintf('Updating TLT %s\n', newTLT);
     

      clear newTLT STACK_PRFX
    end
    
  case 'WriteCsv'
    !mkdir -p csv
    for iTomo = 1:nTomograms
      positionList = geometry.(tomoList{iTomo});
      csvOUT = sprintf('./csv/%s_%s_%s.csv',cycleNumber,tomoList{iTomo},fieldPrefix);
      csvID = fopen(csvOUT, 'w');
      for iSubTomo = 1:size(positionList,1)
        if ( VECTOR_OP(1) == -1 && positionList(iSubTomo,26) ~= -9999 ) || ...
           ( ismember(positionList(iSubTomo,26), VECTOR_OP) )
            fprintf(csvID,'%-06.3f %-06.3f %-06.3f %-04d %-06.3f %-06.3f %-04d %-04d %-06.3f %-06.3f   %-07.3f %-07.3f %-07.3f   %-03.3f %-03.3f %-03.3f    %-03.3f %-03.3f %-03.3f %-03.3f %-03.3f %-03.3f %-03.3f %-03.3f %-03.3f    %-4d\n',...
              positionList(iSubTomo,:));
        end
      end
      fclose(csvID);
    end
  case 'RemoveClasses'
     
    classesToKeep = [];
    classesToKill = [];
    % make a 2d array same size as montage
    clear blankClassMont
    % montage is always square, and the first row is always full, so use this to
    % determine the size.
    sizeMontage = max(locations{end}(2),locations{end}(4))
    blankClassMont(sizeMontage, sizeMontage) = single(0);
    
    % put a 1 at each x,y from the model file identifying classes to kill
    VECTOR_OP = round(VECTOR_OP)
    for iKill = 1:size(VECTOR_OP,1)
      blankClassMont(VECTOR_OP(iKill,1),VECTOR_OP(iKill,2)) = 1;
    end
    
    for iClass = 1:length(locations)
      iCut = blankClassMont(locations{iClass}(1):locations{iClass}(2), ...
                            locations{iClass}(3):locations{iClass}(4));
                          
   
      
      if any(iCut(:))
        classesToKill = [classesToKill, iClass];
      else
        classesToKeep = [classesToKeep, iClass];
      end
    end
    

    fileOUT = fopen(sprintf('%s_ClassMods_%s.txt',cycleNumber, halfSet), 'w');
    fprintf(fileOUT, '%s\n','Classes removed:');
    fprintf(fileOUT, '[%g',classesToKill(1));
    fprintf(fileOUT, ',%g', classesToKill(2:end));
    fprintf(fileOUT, ']\n\n');

    fprintf(fileOUT, '%s\n','Classes retained:');
    fprintf(fileOUT, '[%g',classesToKeep(1));
    fprintf(fileOUT, ',%g', classesToKeep(2:end));
    fprintf(fileOUT, '; 1.*ones(1,%d)]\n\n',numel(classesToKeep));

    
    nTotal = 0;
    nRemoved = 0;
    nRemain = 0;
    for iTomo = 1:nTomograms
      positionList = geometry.(tomoList{iTomo});
    
      
      toRemove = ismember(positionList(:,7),  halfNUM) & ...
                 ismember(positionList(:,26), classesToKill);

      nOrig = ismember(positionList(:,7),  halfNUM) & ...
              (positionList(:,26) ~= -9999);
      positionList(toRemove,26) = -9999;
      nRemain = nRemain + sum(ismember(positionList(:,7),  halfNUM) & ...
              (positionList(:,26) ~= -9999))

      nTotal = nTotal + sum(nOrig);
      nRemoved = nRemoved + sum(toRemove)
      
      geometry.(tomoList{iTomo}) = positionList;
    end

    
   fprintf(fileOUT, '\nremoved:\t%d\nremaining:%d\norig:%d\n',nRemoved,nRemain,nTotal);
    fclose(fileOUT);
  
  case 'ShiftAll'
    % No option to shift eve/odd separately.
    for iTomo = 1:nTomograms
      positionList = geometry.(tomoList{iTomo});
      for iSubTomo = 1:size(positionList,1)
        if ismember(positionList(iSubTomo,7), [1,2]) 
          rotMat = reshape(positionList(iSubTomo,17:25),3,3);
          shifts = (rotMat*shiftXYZ)';
          positionList(iSubTomo,11:13) = positionList(iSubTomo,11:13) + shifts;
        end
      end
      geometry.(tomoList{iTomo}) = positionList;
    end 
  
  case 'ShiftBin'
    for iTomo = 1:nTomograms
      positionList = geometry.(tomoList{iTomo});
      for iSubTomo = 1:size(positionList,1)
        if ismember(positionList(iSubTomo,7) , halfNUM) 
          positionList(iSubTomo,11:13) = positionList(iSubTomo,11:13) + shiftXYZ';
        end
      end
      geometry.(tomoList{iTomo}) = positionList;
    end    
    
  case 'ListTomos'
    system(sprintf('mkdir -p tomoList_%s',cycleNumber));
    tomoListOUT = fopen(sprintf('tomoList_%s/tomoList.txt',cycleNumber),'w');
    for iTomo = 1:nTomograms
      positionList = geometry.(tomoList{iTomo});
        includeList = (positionList(:,26) ~= -9999);
        sTremaining = sum(includeList);
        sTtotal = length(includeList);
        figure('Visible','off'); hist(positionList(includeList,1),floor(sTremaining./5)+1);
        title({sprintf('%s',tomoList{iTomo})},'Interpreter','none'); 
        xlabel('CCC'); 
        file_out = sprintf('tomoList_%s/%s.pdf', cycleNumber, tomoList{iTomo});
        saveas(gcf, file_out,'pdf')   
        fprintf(tomoListOUT,'%s\t%d/%d\n',tomoList{iTomo},sTremaining,sTtotal);
    end 
    
  case 'RemoveTomos'
    
    if (listTomos)
      f=fieldnames(masterTM.mapBackGeometry.tomoName);
      tomoFid = fopen('tomoList.txt','w');
      fprintf(tomoFid,'%s\n',f{:});
      fclose(tomoFid);   
    else
      tomoList = importdata(VECTOR_OP);
      f=fieldnames(masterTM.mapBackGeometry.tomoName);
    
      for iOrig = 1:length(f)
        flgRemove = 1;
        for iToKeep = 1:length(tomoList)         
          if strcmp(f{iOrig},tomoList{iToKeep})
            tomoList{iToKeep};
            flgRemove = 0;
            break;
          end
        end
        if (flgRemove)
          if isfield(masterTM.reconGeometry.(f{iOrig}))
            masterTM.reconGeometry = rmfield(masterTM.reconGeometry,f{iOrig});
          end
          if isfield(masterTM.tiltGeometry.(f{iOrig}))
            masterTM.tiltGeometry = rmfield(masterTM.tiltGeometry,f{iOrig});
          end

          if isfield(masterTM.(cycleNumber).RawAlign.f{iOrig})
            masterTM.(cycleNumber).RawAlign = rmfield(masterTM.(cycleNumber).RawAlign,(f{iOrig}));
          end
          if isfield(masterTM.mapBackGeometry.tomoName.(f{iOrig}))
            tN = masterTM.mapBackGeometry.tomoName.(f{iOrig}).tomoNumber;
            tName = masterTM.mapBackGeometry.tomoName.(f{iOrig}).tiltName;
            if isfield(masterTM.mapBackGeometry,tName)
              masterTM.mapBackGeometry.(tName).nTomos = masterTM.mapBackGeometry.(tName).nTomos - 1;
                keepRows = ismember(1:size(masterTM.mapBackGeometry.(tName).coords,1),tN);
                masterTM.mapBackGeometry.(tName).coords = masterTM.mapBackGeometry.(tName).coords(~keepRows,:);
                masterTM.mapBackGeometry.tomoName = rmfield(masterTM.mapBackGeometry.tomoName,f{iOrig});
              if masterTM.mapBackGeometry.(tName).nTomos == 0 
                masterTM.mapBackGeometry = rmfield(masterTM.mapBackGeometry,tName);
              end
            end
          end
        
        end
      end
      geometry = masterTM.tiltGeometry;  
    end   
    
      
       
    
  case 'RandomizeEulers'
    for iTomo = 1:nTomograms
      rng('shuffle')
      positionList = geometry.(tomoList{iTomo});
      for iSubTomo = 1:size(positionList,1)
        if ismember(positionList(iSubTomo,7) , halfNUM)
          
          rotMat = reshape(positionList(iSubTomo,17:25),3,3);
          
          r1 = randi([-1,1].*VECTOR_OP(1));
          r2 = randi([-1,1].*VECTOR_OP(2));
          r3 = randi([-1,1].*VECTOR_OP(3));
          
          randMat = BH_defineMatrix([r1,r2,r3-r1],'SPIDER','inv');
          
          positionList(iSubTomo,17:25) = reshape(rotMat*randMat,1,9);
        end
      end
      geometry.(tomoList{iTomo}) = positionList;
    end   
    
  case 'RemoveFraction'
    % Get distribution of CCC from given cycle rawAlignment
    % Save histogram, also remove given bottom percentage and report CCC cutoff
    if ~(strcmpi(STAGEofALIGNMENT, 'RawAlignment'))
      error('Can only remove fraction at RawAlignment Stage')
    end
    cccValues = [];
    for iTomo = 1:nTomograms
      positionList = geometry.(tomoList{iTomo});
      cccValues = [cccValues;positionList(( positionList(:,26) ~= -9999 ), 1)];
    end   
    
    nSubTomos = length(cccValues);
    
    sortCCC = sort(cccValues,'descend')
    
    if percentCut < 1
      cutIDX = floor(percentCut .* nSubTomos) + 1
    else
      % Keep the top N values
      cutIDX = percentCut
    end
    
    cccCutOff = sortCCC(cutIDX)
    
    nRemoved = 0;
    for iTomo = 1:nTomograms
      positionList = geometry.(tomoList{iTomo});
      removeList = positionList(:, 1) < cccCutOff & positionList(:,26)~=-9999;
      nRemoved = nRemoved + sum(removeList)
      positionList(removeList,26) = -9999;
      geometry.(tomoList{iTomo}) = positionList;
    end   
    
    figure('visible', 'off'), hist(cccValues);
    title(sprintf('CCC distribution\n%2.2f%% cut\n %f CCC\n%d/ %d removed ',100*percentCut,cccCutOff,nRemoved,nSubTomos));
    xlabel('CCC'); ylabel('nSubTomos');
    
    file_out = sprintf('%s-cccCutoff.pdf', outputPrefix);

    saveas(gcf, file_out,'pdf')
    
  case 'ListPercentiles'
    cccVector = [];
    % Gather included scores
    for iTomo = 1:nTomograms
      positionList = geometry.(tomoList{iTomo});
      cccVector = [cccVector;positionList(( positionList(:,26) ~= -9999 ), 1)];
    end   
    nVols = length(cccVector)
    cccVector = sort(cccVector,1,'descend');
    snrVector = 2.* cccVector ./ (1-cccVector);
    percentiles = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9;0,0,0,0,0,0,0,0,0];
    percentiles(2,:) = cccVector(floor(nVols.*(percentiles(1,:))))
    weightedCurve = zeros(nVols,2);
    weightedCurve(:,1) = 1:nVols;
    for iVol = 1:nVols
      weightedCurve(iVol,2) = mean(snrVector(1:iVol,1).*sqrt(iVol));
    end
    length([1:nVols]')
    length(weightedCurve(:,2))
    r = fit([1:nVols]',weightedCurve(:,2),'linear');
    maxVal = find(r(2:nVols+1)-r(1:nVols) <= 0, 1,'first')
    maxCCC = cccVector(maxVal)
    
    figure('visible', 'off'), plot(weightedCurve(:,2));
    title({'mean SNR via CCC vs nVols',sprintf('%0.3f CCC',maxCCC),sprintf( '%d Vols', maxVal),sprintf('%0.2f',100.*maxVal/nVols)});
    file_out = sprintf('%s_percentiles_%s.pdf', cycleNumber, fieldPrefix);
    saveas(gcf, file_out,'pdf')
    
    fOUT = sprintf('./%s_percentiles_%s.txt',cycleNumber,fieldPrefix);
    fID = fopen(fOUT, 'w');    
    fprintf(fID,'%-4.3f\t%-4.3f\t%-4.3f\t%-4.3f\t%-4.3f\t%-4.3f\t%-4.3f\t%-4.3f\t%-4.3f\t\n%-4.3f\t%-4.3f\t%-4.3f\t%-4.3f\t%-4.3f\t%-4.3f\t%-4.3f\t%-4.3f\t%-4.3f\t\n',percentiles');
    fclose(fID);
  otherwise
    error('OPERATION must be WriteCsv, RemoveClasses, ShiftAll, RemoveFraction not %s', OPERATION)
end 

% Redundant for WriteCsv, otherwise update the new geometry, which was backed up
% at the begining as Pre_OPERATION_...
switch STAGEofALIGNMENT
  
  case 'TiltAlignment'
    masterTM.tiltGeometry = geometry;
  
  case 'RawAlignment'
    masterTM.(cycleNumber).RawAlign = geometry;
  case 'Cluster'

    masterTM.(cycleNumber).(clusterGeom) = geometry;

  otherwise
    error('STAGE_ALIGNMENTS: [Tilt,Class,No,Raw]Alignment, not %s', ...
                                                        STAGEofALIGNMENT);
end

    subTomoMeta = masterTM;
    save(emc.('subTomoMeta'), 'subTomoMeta');
end





