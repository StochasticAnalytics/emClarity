function [] = BH_combineProjects( outputName, inp )

% inp= {'gagT8I_1.mat','gagT8I_2.mat'};

nProj = length(inp);

% Folder to copy the local xforms into - local xforms will need to be
% combined 
system('mkdir -p tmpMB');

for iProj = 1:nProj
  
  if iProj == 1
    masterTM = load(inp{1});
    masterCycle = masterTM.subTomoMeta.currentCycle;
    masterTomoCPR =   masterTM.subTomoMeta.currentTomoCPR;
    
    system(sprintf('cp mapBack%d/*.local tmpMB',masterTomoCPR));
    
    % Keep only the last cycle that has alignment information.
    while ~isfield(masterTM.subTomoMeta.(sprintf('cycle%0.3d',masterCycle)),'RawAlign') && masterCycle > 0
      masterCycle = masterCycle - 1;
    end
    
    fprintf('\nFor project %s the last cycle with alignment info is %d\n',inp{iProj},masterCycle);
  
    for iCycle = masterCycle-1:-1:0
      masterTM.subTomoMeta = rmfield(masterTM.subTomoMeta, ...
                             (sprintf('cycle%0.3d',iCycle)));
      
    end
    continue;
  else
    tmpTM = load(inp{iProj});
    system(sprintf('cp mapBack%d/*.local tmpMB',masterTomoCPR));
  end
  
  if (masterTM.subTomoMeta.maxGoldStandard ~= tmpTM.subTomoMeta.maxGoldStandard)
    fprintf('\n\nWarning different template matching resolutions used. Selecting higher.\n')
    fprintf('\nof %f and %f\n',masterTM.subTomoMeta.maxGoldStandard,tmpTM.subTomoMeta.maxGoldStandard);
    masterTM.subTomoMeta.maxGoldStandard = min(masterTM.subTomoMeta.maxGoldStandard,tmpTM.subTomoMeta.maxGoldStandard);
  end
  
  if (masterTM.subTomoMeta.CUTPADDING ~= tmpTM.subTomoMeta.CUTPADDING)
    error('\nThe cutpadding does not match, what?!\n');
  end

  % Take the best resolution for determing how many slabs to reconstruct.
  masterTM.subTomoMeta.currentResForDefocusError = ...
    min(masterTM.subTomoMeta.currentResForDefocusError, ...
           tmpTM.subTomoMeta.currentResForDefocusError);
  
  % Accumulate the starting numbers
  masterTM.subTomoMeta.nSubTomoInitial = ...
    masterTM.subTomoMeta.nSubTomoInitial + ...
       tmpTM.subTomoMeta.nSubTomoInitial;
     
  % Set the current tomoCPR to the most advanced
  if (masterTM.subTomoMeta.currentTomoCPR < tmpTM.subTomoMeta.currentTomoCPR)
    masterTM.subTomoMeta.currentTomoCPR = tmpTM.subTomoMeta.currentTomoCPR;
  end
  
  % Add the mapBack
  tmpTM = load(inp{iProj});
  fn = fieldnames(tmpTM.subTomoMeta.mapBackGeometry);
  for iField = 1:length(fn)
    if strcmpi(fn{iField},'tomoName')
      fnName = fieldnames(tmpTM.subTomoMeta.mapBackGeometry.tomoName);
      for iName = 1:length(fnName)
        masterTM.subTomoMeta.mapBackGeometry.tomoName.(fnName{iName}) = ...
          tmpTM.subTomoMeta.mapBackGeometry.tomoName.(fnName{iName}) ;
      end
    else
      masterTM.subTomoMeta.mapBackGeometry.(fn{iField}) = ...
        tmpTM.subTomoMeta.mapBackGeometry.(fn{iField});
    end
  end
  
  tomoNames = fieldnames(tmpTM.subTomoMeta.reconGeometry);
  nFields = length(tomoNames);
  for iField = 1:nFields
    
    masterTM.subTomoMeta.reconGeometry.(tomoNames{iField}) = ...
      tmpTM.subTomoMeta.reconGeometry.(tomoNames{iField});
    
    masterTM.subTomoMeta.tiltGeometry.(tomoNames{iField}) = ...
      tmpTM.subTomoMeta.tiltGeometry.(tomoNames{iField});  
    
    
  end
  
    tmpCycle = tmpTM.subTomoMeta.currentCycle;
    
    % Keep only the last cycle that has alignment information.
    while ~isfield(tmpTM.subTomoMeta.(sprintf('cycle%0.3d',tmpCycle)),'RawAlign') && tmpCycle > 0
      tmpCycle = tmpCycle - 1;
    end
    
    fprintf('\nFor project %s the last cycle with alignment info is %d\n',inp{iProj},tmpCycle);
    
    for iField = 1:nFields
      masterTM.subTomoMeta.(sprintf('cycle%0.3d',masterCycle)).RawAlign.(tomoNames{iField}) = ...   
         tmpTM.subTomoMeta.(sprintf('cycle%0.3d',tmpCycle)).RawAlign.(tomoNames{iField});
    end
  
end

% Now reset the cycle
masterTM.subTomoMeta.currentCycle = 0;
 
% Increment the tomoCPR by one so no older files are overwritten

masterTM.subTomoMeta.currentTomoCPR = masterTM.subTomoMeta.currentTomoCPR + 1;

system(sprintf('mkdir -p mapBack%d',masterTM.subTomoMeta.currentTomoCPR));


% Write the tilt geometry, for simplicity right now, just overwrite (one
% per tomo)

tiltList = fieldnames(masterTM.subTomoMeta.tiltGeometry);
for iTilt = 1:length(tiltList)
  
  tiltName = sprintf('fixedStacks/ctf/%s_ali%d_ctf.tlt', ...
                      masterTM.subTomoMeta.mapBackGeometry.tomoName.(tiltList{iTilt}).tiltName, ...
                      masterTM.subTomoMeta.currentTomoCPR + 1);
  thisTilt = fopen(tiltName,'w');
  
  fprintf(thisTilt,['%d\t%08.2f\t%08.2f\t%07.3f\t%07.3f\t%07.3f\t%07.7f\t%07.7f\t',...
             '%07.7f\t%07.7f\t%5e\t%5e\t%5e\t%7e\t%5e\t%5e\t%5e\t%5e\t%5e\t',...
             '%d\t%d\t%d\t%8.2f\n'], masterTM.subTomoMeta.tiltGeometry.(tiltList{iTilt})');
  fclose(thisTilt);
  
  rawTLT = sprintf('mapBack%d/%s_ali%d_ctf.tlt', ...
                      masterTM.subTomoMeta.currentTomoCPR, ...
                      masterTM.subTomoMeta.mapBackGeometry.tomoName.(tiltList{iTilt}).tiltName, ...
                      masterTM.subTomoMeta.currentTomoCPR);
                    
  thisRawTLT = fopen(rawTLT,'w');
  tiltAngles = sortrows(masterTM.subTomoMeta.tiltGeometry.(tiltList{iTilt})(:,[1,4]),1);
  fprintf(thisRawTLT,'%f\n',tiltAngles(:,2)');
  fclose(thisRawTLT);

  xf = sprintf('mapBack%d/%s_ali%d_ctf.tltxf', ...
                      masterTM.subTomoMeta.currentTomoCPR, ...
                      masterTM.subTomoMeta.mapBackGeometry.tomoName.(tiltList{iTilt}).tiltName, ...
                      masterTM.subTomoMeta.currentTomoCPR);
                    
  thisXF = fopen(xf,'w');
  for i = 1:length(tiltAngles)
    fprintf(thisXF,'1.0 0.0 0.0 1.0 0.0 0.0\n');
  end
  fclose(thisXF);
  
end
 
outputName = strsplit(outputName, '.mat');
subTomoMeta = masterTM.subTomoMeta;
save(sprintf('%s.mat',outputName{1}),'subTomoMeta');


end
