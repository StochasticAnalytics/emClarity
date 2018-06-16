function BH_geometry_Constraints(PARAMETER_FILE, CYCLE, distCut, angCut, latticeNumber)


pBH = BH_parseParameterFile(PARAMETER_FILE);
load(sprintf('%s.mat', pBH.('subTomoMeta')), 'subTomoMeta');


angCut = [0,str2num(angCut)]
distCut = str2num(distCut)
modBin = 7;
nNeighbors = str2num(latticeNumber)
nTotal = 0;
CYCLE = str2num(CYCLE);

cycleNumber = sprintf('cycle%0.3u', CYCLE);

if (CYCLE)
  geom  = subTomoMeta.(cycleNumber).RawAlign;
else
  geom = subTomoMeta.('cycle000').('geometry');
end

baseNum = fieldnames(geom);


for i = 1:length(baseNum)
  i

listIN = geom.(baseNum{i});
size(listIN)
return
nVol = size(listIN,1);
keepList = zeros(nVol,2);


  for iPt = 1:nVol
    iPt
    if ( listIN(iPt,26) == -9999 )
      distVect = sqrt((listIN(:,11)-listIN(iPt,11)).^2 + ...
                    (listIN(:,12)-listIN(iPt,12)).^2 + ...
                    (listIN(:,13)-listIN(iPt,13)).^2);

      closeVect = find(distVect < distCut);

      % make sure there are at least 4 other volumes within 90 angstrom
      flgRemove = 0;
      if (length(closeVect)>=nNeighbors)
      nAngClose = 0;
      particleAxis = reshape(listIN(iPt,17:25),3,3)*[0;0;1];

       for iAng = 1:length(closeVect)
         iAxis = reshape(listIN(closeVect(iAng),17:25),3,3)*[0;0;1];
         iAngDiff = acosd(dot(particleAxis,iAxis));
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

  % Remove the positions
  geom.(baseNum{i})(~keepList,26) = -9999;
  keepList = keepList(keepList~=0);

%modOUT = listIN(keepList,11:13)./modBin;
%modFileOUT = fopen('tmp.txt','w');
%fprintf(modFileOUT,'%f %f %f\n',modOUT');
%fclose(modFileOUT);

nTotal = nTotal + length(keepList)
end
%system(sprintf('point2model -circle 3 -sphere 2 -scat -thick 2 -color 80,191,255 -image ../cache/%s.rec tmp.txt %s.mod',modNameIn,modNameIn));
  

if (CYCLE)
  subTomoMeta.(cycleNumber).RawAlign = geom;
else
  subTomoMeta.('cycle000').('geometry') = geom;
end

save(sprintf('%s.mat', pBH.('subTomoMeta')),'subTomoMeta');

end % end of function
