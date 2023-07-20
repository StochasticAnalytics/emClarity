function BH_geometry_Constraints(PARAMETER_FILE, nWorkers, distCut, angCut, latticeNumber, varargin)



pixelSize = PARAMETER_FILE;
nPeaks = 1;
avgRadius = [0,0];
if nargin < 6
  flgCentroid = false;
else
  flgCentroid = true;
  avgRadius = str2num(varargin{1});
  nPeaks = str2num(varargin{2});
  
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
nWorkers = EMC_str2double(nWorkers);
% delete(gcp('nocreate'));
% EMC_parpool(nWorkers);


iterList = cell(nWorkers);
geomList = cell(nWorkers);
nTotal = cell(nWorkers);
for i = 1:nWorkers
  geomList{i} = geom;
  iterList{i} = i:nWorkers:nModFiles;
  nTotal{i} = 0;
end

for iPar = 1:nWorkers
% for iPar = 1:nWorkers

  geom = geomList{iPar};
  
  for i = iterList{iPar}
    
    listIN = gpuArray(load(geom{i}));
    if (flgCentroid)
      [~,p,~]= fileparts(geom{i});
      b = strsplit(p,'_');
      reconName = strjoin(b(1:end-1),'_');
      reconDims = load(sprintf('recon/%s_recon.txt',reconName));
      reconOrigin = reshape(floor(reconDims(1:3)./2) + 1,1,3);
    end


    nVol = size(listIN,1);
    keepList = zeros(nVol,1);
    peakList = ones(nVol,1);

    if (flgCentroid)
      rad_cutoff_sigma = 3;
      com_iter_max = 35;
      % First loop over and get rid of any points that are not directed
      % radially outward.

        % Tomogram may not be centered, either because the estimate in Z is
        % off, or the particle is near the edge and the box was shifted.
        % First loop over selecting likely candidates, and determine their
        % center of mass, then re-loop using this.
        COM = [0,0,0];
        nCOM = 0;
        % Adjust for outliers to the radius.
        mean_radius = 0;
        mean_radius_sq = 0;
        final_loop = false;
        rad_cutoff = 1e6;
        orig_reconOrigin = reconOrigin;
        for comLoop = 1:com_iter_max
          if comLoop > 1
            COM = [COM ./ nCOM] ./ pixelSize;
            calc_avg = mean_radius ./ nCOM;
            calc_std = sqrt(mean_radius_sq ./ nCOM - (mean_radius ./ nCOM)^2);
            rad_cutoff = calc_avg + rad_cutoff_sigma*calc_std;
            if norm(COM) < 2
              final_loop = true;
            end
%             fprintf('Adjusting the origin from recon %f,%f,%f\nTo %f,%f,%f\n',reconOrigin,COM);
            reconOrigin = reconOrigin + COM;
            COM = [0,0,0];
            nCOM = 0;
            mean_radius = 0;
            mean_radius_sq = 0;

           
          end
          for iPt = 1:nVol
            % It might not be the top peak
            foundOnePeak = false;
            for iPeak = 1:nPeaks
              % Vector that points along the Z-axis of the particle
              particleAxis = gather(reshape(listIN(iPt,[17:25] + 26*(iPeak-1)),3,3)*[0;0;1]);
              % Vector that points from the recon origin to the subtomo origin
              particleCoords = gather(pixelSize*(listIN(iPt,[11:13] + 26*(iPeak-1)) - reconOrigin));
              if (norm(particleCoords) > avgRadius(2) || norm(particleCoords) < avgRadius(1) || norm(particleCoords) > rad_cutoff)
                break;
              end
              if acosd(dot(particleCoords./norm(particleCoords), particleAxis)) < 45 
                if final_loop
                  foundOnePeak = true; 
                  peakList(iPt) = iPeak;
                  keepList(iPt) = iPt;
                end
                COM = COM + particleCoords;
                particleWeight = 1;%listIN(iPt,1 + 26*(iPeak-1));
                nCOM = nCOM + particleWeight;
                mean_radius = mean_radius + particleWeight.*norm(particleCoords);
                mean_radius_sq = mean_radius_sq + particleWeight.*norm(particleCoords).^2;
                break;
              end
            end
             
            if final_loop
              if ~(foundOnePeak)
                listIN(iPt,26:26:26*nPeaks) = -9999 ;
%                 keepList(iPt) = iPt;       
              end
            end
          end
          if final_loop
%             % Until I can add something linke this to spike_constraint in
%             % average just save the offset. Better yet, I'll fit an ellipse
%             % or something.
%             com_offset = fopen(sprintf('recon/%s_comOffset.txt',reconName),'w');
%             fprintf(com_offset,'%3.3f %3.3f %3.3f',reconOrigin - orig_reconOrigin);
%             fclose(com_offset);
            break
          end
        end % loop to deterimne center of mass

    else

      for iPt = 1:nVol

        if flgCentroid 
          thisPeak = peakList(iPt);
        else
          thisPeak = 1;
        end

        for iPeak = thisPeak
          if ( listIN(iPt,26 + 26*(iPeak-1)) ~= -9999 )
            
            c1 = sub2ind(size(listIN),[1:size(listIN,1)]',11 + 26.*(peakList-1));
            c2 = sub2ind(size(listIN),[1:size(listIN,1)]',12 + 26.*(peakList-1));
            c3 = sub2ind(size(listIN),[1:size(listIN,1)]',13 + 26.*(peakList-1));

            distVect = sqrt( (listIN(c1)-listIN(iPt,11 + 26*(iPeak-1))).^2 + ...
                             (listIN(c2)-listIN(iPt,12 + 26*(iPeak-1))).^2 + ...
                             (listIN(c3)-listIN(iPt,13 + 26*(iPeak-1))).^2 );

            closeVect = find(distVect < distCut);
            if (length(closeVect)>=nNeighbors)
            nAngClose = 0;
            particleAxis = reshape(listIN(iPt,[17:25] + 26*(iPeak-1)),3,3)*[0;0;1];

             for iAng = 1:length(closeVect)
               iAxis = reshape(listIN(closeVect(iAng),17:25),3,3)*[0;0;1];
               iAngDiff = dot(particleAxis,iAxis);
               if (abs(iAngDiff) > 1) 
                 iAngDiff = fix(iAngDiff);
               end
               if abs(iAngDiff) < abs(angCut(2)) && abs(iAngDiff) > abs(angCut(1))
                 nAngClose = nAngClose + 1;
               end
             end
              if nAngClose > nNeighbors-1
                keepList(iPt) = iPt;
              else
                keepList(iPt) = 0;
              end
            end

          end
        end
      end
    
    end

      listIN = listIN(find(keepList),:);
      [~,fileName,~] = fileparts(geom{i});

      if isempty(listIN)
        fprintf('NO subtomos retained! for %s\n',fileName);
      else
        csvOUT = fopen(sprintf('convmap/%s.txt',fileName),'w');
        h = hist(listIN(:,15))
        for iLine = 1:size(listIN,1)

          fprintf(csvOUT,'%4.4f %4.4f %4.4f\n',listIN(iLine,[11:13])./listIN(iLine,2));
        end
        fclose(csvOUT);

        system(sprintf('point2model -number 1 -circle 3 -sphere 3 -scat  -color 80,191,255 convmap/%s.txt convmap/%s.mod',fileName,fileName));
        nTotal{iPar} = nTotal{iPar}  + sum(keepList~=0);
      end
  end
%

end
all_total = 0;
for iWorker = 1:nWorkers
  all_total = all_total + nTotal{iWorker};
end
  fprintf('Total points found %d\n', all_total);
  
end % end of function
