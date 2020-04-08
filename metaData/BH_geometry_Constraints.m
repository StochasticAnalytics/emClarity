function BH_geometry_Constraints(PARAMETER_FILE, CYCLE, distCut, angCut, latticeNumber, varargin)



pixelSize = PARAMETER_FILE;
  
if nargin < 6
  flgCentroid = false;
else
  flgCentroid = true;
  avgRadius = str2num(varargin{1});
  if avgRadius == 1
    avgRadius = [0,avgRadius./pixelSize];
  else
    avgRadius = [avgRadius(1)./pixelSize,avgRadius(2)./pixelSize];
  end
end

if ~isfolder('convmap')
  error('No convmap directory found');
end
[failedCSV,csvList] = system('ls convmap/*.csv'); % add explicit convmap check TODO
if (failedCSV)
  error('Did not find your csv files in convmap folder');
end

geom = strsplit(csvList);

keepCSV = true(1,length(geom));
for iCSV = 1:length(keepCSV)
  if isempty(geom{iCSV})
    keepCSV(iCSV) = false;
  end
end
geom = geom(keepCSV);
nModFiles = length(geom);
  

angCut = [0,str2num(angCut)];
distCut = str2num(distCut)./pixelSize;
nNeighbors = str2num(latticeNumber);

% FIXME from param file.
nWorkers = 8;
delete(gcp('nocreate'));
parpool(nWorkers);
iterList = cell(nWorkers);
geomList = cell(nWorkers);
nTotal = cell(nWorkers);
for i = 1:nWorkers
  geomList{i} = geom;
  iterList{i} = i:nWorkers:nModFiles;
  nTotal{i} = 0;
end

parfor iPar = 1:nWorkers
  
  geom = geomList{iPar};
  
  for i = iterList{iPar}
    
    listIN = gpuArray(load(geom{i}));
    if (flgCentroid)
      [~,p,~]= fileparts(geom{i});
      b = strsplit(p,'_');
      reconName = strjoin(b(1:end-1),'_')
      reconDims = load(sprintf('recon/%s_recon.txt',reconName));
      reconOrigin = reshape(floor(reconDims(1:3)./2) + 1,1,3);
    end


    nVol = size(listIN,1);
    keepList = zeros(nVol,1);


    if (flgCentroid)
      % First loop over and get rid of any points that are not directed
      % radially outward.

        for iPt = 1:nVol
          particleAxis = reshape(listIN(iPt,17:25),3,3)*[0;0;1];
          particleCoords = listIN(iPt,11:13) - reconOrigin;
          particleNorm = norm(particleCoords);
          if (particleNorm > avgRadius(2) || particleNorm < avgRadius(1) || dot(particleCoords./particleNorm, particleAxis) < 0.75)
            listIN(iPt,26) = -9999 ;
          end

        end
    end

      for iPt = 1:nVol


        if ( listIN(iPt,26) ~= -9999 )
          distVect = sqrt(((listIN(:,11)-listIN(iPt,11)).^2 + ...
                        (listIN(:,12)-listIN(iPt,12)).^2 + ...
                        (listIN(:,13)-listIN(iPt,13)).^2));

          closeVect = find(distVect < distCut);

          if (length(closeVect)>=nNeighbors)
          nAngClose = 0;
          particleAxis = reshape(listIN(iPt,17:25),3,3)*[0;0;1];

           for iAng = 1:length(closeVect)
             iAxis = reshape(listIN(closeVect(iAng),17:25),3,3)*[0;0;1];
             iDot = dot(particleAxis,iAxis);
             if (abs(iDot) > 1 && abs(iDot) < 1.0001)
               iDot = fix(iDot);
             end
=             iAngDiff = acosd(iDot);
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
  
  


      listIN = listIN(find(keepList),:);
      [~,fileName,~] = fileparts(geom{i});

      csvOUT = fopen(sprintf('convmap/%s.txt',fileName),'w');
      h = hist(listIN(:,15))
      for iLine = 1:size(listIN,1)

        fprintf(csvOUT,'%4.4f %4.4f %4.4f\n',listIN(iLine,[11:13])./listIN(iLine,2));
      end
      fclose(csvOUT);

      system(sprintf('point2model -number 1 -circle 3 -sphere 3 -scat  -color 80,191,255 convmap/%s.txt convmap/%s.mod',fileName,fileName));

    nTotal{iPar} = nTotal{iPar}  + sum(keepList~=0);
  end
%

end
all_total = 0;
for iWorker = 1:nWorkers
  all_total = all_total + nTotal{iWorker};
end
  fprintf('Total points found %d\n', all_total);
  
end % end of function
