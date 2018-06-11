function [ outputRefs ] = BH_multi_combineLowResInfo( inputRefs, inputCounts, pixelSize, resCutOff )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

refIDX = BH_multi_isCell( inputRefs{1} )
nRefs = length(refIDX)

[radialGrid,~,~,~,~,~] = BH_multi_gridCoordinates([512,512,512],...
                            'Cartesian','cpu',{'none'},1,0,1);
radialGrid = radialGrid ./ pixelSize;

outputRefs = cell(2,1);
outputRefs{1} = cell(nRefs,1);
outputRefs{2} = cell(nRefs,1);

for iRef = refIDX'
  % For today assume equal contributions, I think I already save this in the
  % meta data, so add this in soon.
%   oddWeight = sum(nExtracted(iClassPos,1)) ./ sum(nExtracted(iClassPos,1:2))
%   eveWeight = sum(nExtracted(iClassPos,2)) ./ sum(nExtracted(iClassPos,1:2))
  oddWeight = inputCounts{1}(2,iRef) ./ (inputCounts{1}(2,iRef) + inputCounts{2}(2,iRef))
  eveWeight = inputCounts{2}(2,iRef) ./ (inputCounts{1}(2,iRef) + inputCounts{2}(2,iRef))

  
  
  [ combPAD ] = BH_multi_padVal( size(inputRefs{1}{iRef}), 512 );
  % Oversample so the cutoff is more accurate, and use double precision for the
  % same reason.
  oddPAD = fftn(BH_padZeros3d(inputRefs{1}{iRef}, ...
                              combPAD(1,:),combPAD(2,:),'cpu','doubleTaper'));
  evePAD = fftn(BH_padZeros3d(inputRefs{2}{iRef}, ...
                              combPAD(1,:),combPAD(2,:),'cpu','doubleTaper'));

  sharedInfo = (radialGrid < 1/resCutOff) .* (oddWeight.*oddPAD + eveWeight.*evePAD);

  oddPAD = real(ifftn(sharedInfo + (radialGrid >= 1/resCutOff).*oddPAD));
  evePAD = real(ifftn(sharedInfo + (radialGrid >= 1/resCutOff).*evePAD));
  clear sharedInfo oddWeight eveWeight

  outputRefs{1}{iRef} = single(oddPAD(1+combPAD(1,1):end-combPAD(2,1), ...
                                      1+combPAD(1,2):end-combPAD(2,2), ...
                                      1+combPAD(1,3):end-combPAD(2,3)));
  clear oddPAD                               

  outputRefs{2}{iRef} = single(evePAD(1+combPAD(1,1):end-combPAD(2,1), ...
                                      1+combPAD(1,2):end-combPAD(2,2), ...
                                      1+combPAD(1,3):end-combPAD(2,3)));
  clear evePAD
end          
clear radialGrid inputRefs
end

