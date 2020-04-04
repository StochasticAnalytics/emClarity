function BH_geometry_Constraints(PARAMETER_FILE, CYCLE, distCut, angCut, latticeNumber)


flgCSV = 0;
CYCLE = str2num(CYCLE);

pixelSize = str2double(PARAMETER_FILE);
% % % if (isnan(pixelSize))
  pixelSize = PARAMETER_FILE;
  flgCSV = 1;
  if ~isdir('convmap')
    error('No convmap directory found');
  end
  [failedCSV,csvList] = system('ls convmap/*.csv'); % add explicit convmap check TODO
  if (failedCSV)
    error('Did not find your csv files in convmap folder');
  end
  geom = strsplit(csvList)
  nModFiles = length(geom)-1
  
% % % else
% % %   pBH = BH_parseParameterFile(PARAMETER_FILE);
% % %   load(sprintf('%s.mat', pBH.('subTomoMeta')), 'subTomoMeta');
% % %   pixelSize = pBH.('PIXEL_SIZE')*10^10
% % %   
% % % 
% % %   cycleNumber = sprintf('cycle%0.3u', CYCLE);
% % % 
% % %   if (CYCLE)
% % %     geom  = subTomoMeta.(cycleNumber).RawAlign;
% % %   else
% % %     if isfield(subTomoMeta.('cycle000'),'geometry_orig')
% % %       geom = subTomoMeta.('cycle000').('geometry_orig');
% % %     else
% % %       subTomoMeta.('cycle000').('geometry_orig') = ...
% % %       subTomoMeta.('cycle000').('geometry');
% % %       geom = subTomoMeta.('cycle000').('geometry_orig');
% % %     end
% % %   end
% % %   baseNum = fieldnames(geom);
% % %   nModFiles = length(baseNum);
% % %   binVal = pBH.('Tmp_samplingRate');
% % % end


angCut
distCut
pixelSize

angCut = [0,str2num(angCut)]
distCut = str2num(distCut)./pixelSize
modBin = 7;
nNeighbors = str2num(latticeNumber)
nTotal = 0;





for i = 1:nModFiles
  i

  if flgCSV
    listIN = gpuArray(load(geom{i}));
  else
      listIN = gpuArray(geom.(baseNum{i}));
  end

nVol = size(listIN,1);
keepList = zeros(nVol,1);



  for iPt = 1:nVol
    
 
    if ( listIN(iPt,26) ~= -9999 )
      distVect = sqrt(((listIN(:,11)-listIN(iPt,11)).^2 + ...
                    (listIN(:,12)-listIN(iPt,12)).^2 + ...
                    (listIN(:,13)-listIN(iPt,13)).^2));

      closeVect = find(distVect < distCut);

      % make sure there are at least 4 other volumes within 90 angstrom
      flgRemove = 0;
      if (length(closeVect)>=nNeighbors)
      nAngClose = 0;
      particleAxis = reshape(listIN(iPt,17:25),3,3)*[0;0;1];

       for iAng = 1:length(closeVect)
         iAxis = reshape(listIN(closeVect(iAng),17:25),3,3)*[0;0;1];
         iDot = dot(particleAxis,iAxis);
         if (abs(iDot) > 1 && abs(iDot) < 1.0001)
           iDot = fix(iDot);
         end
         iDot;
         iAngDiff = acosd(iDot);
         if abs(iAngDiff) < abs(angCut(2)) && abs(iAngDiff) > abs(angCut(1))
           nAngClose = nAngClose + 1;
         end
       end
        if nAngClose > nNeighbors-1
          keepList(iPt) = iPt;
        end
      end

    end
  end
  
  

  if ( flgCSV )
    listIN = listIN(find(keepList),:);
    [~,fileName,~] = fileparts(geom{i});

    csvOUT = fopen(sprintf('convmap/%s.txt',fileName),'w');
    h = hist(listIN(:,15))
    for iLine = 1:size(listIN,1)

      fprintf(csvOUT,'%4.4f %4.4f %4.4f\n',listIN(iLine,[11:13])./listIN(iLine,2));
    end
    fclose(csvOUT);
    
    system(sprintf('point2model -number 1 -circle 3 -sphere 3 -scat  -color 80,191,255 convmap/%s.txt convmap/%s.mod',fileName,fileName));
  else
    keepList = keepList(keepList ~= 0);
    csvOUT = fopen(sprintf('convmap/%s_cleaned.txt',baseNum{i}),'w');
    
    % Remove the positions
    if (CYCLE) 
      geom.(baseNum{i})(~keepList,26) = -9999;
    else
      geom.(baseNum{i})= geom.(baseNum{i})(keepList,:);
    end
      fprintf(csvOUT,'%6.6f %6.6f %6.6f\n',geom.(baseNum{i})(:,11:13)'./binVal);
      fclose(csvOUT);
      system(sprintf('point2model -number 1 -circle 3 -sphere 3 -scat  -color 80,191,255 convmap/%s_cleaned.txt convmap/%s.mod-cleaned',baseNum{i},baseNum{i}));
  end
  nTotal = nTotal + sum(keepList~=0);
end
%
  nTotal
if ~flgCSV
  if (CYCLE)
    subTomoMeta.(cycleNumber).RawAlign = geom;
  else
    subTomoMeta.('cycle000').('geometry') = geom;
  end

  save(sprintf('%s.mat', pBH.('subTomoMeta')),'subTomoMeta');
end

end % end of function
