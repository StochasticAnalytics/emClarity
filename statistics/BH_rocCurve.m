function [] = BH_rocCurve(gt_name, mip_name,symmetry,nDefocus, pixelSize, radius, thr)
%pixelSize=1.35 % Ang
%radius = 180 % Ang
%thr = 0.05 % percent radius considered a hit

if ~exist(gt_name,'file')
  error('Did not find your PDB %s\n',gt_name);
end
if ~exist(mip_name,'file')
  error('Did not find your MIP %s\n',mip_name);
end
switch symmetry
  case 'C1'
    nAngles = 1595420;
  case 'O'
    nAngles = 89411;
  otherwise
    error('did not understand C1 or O from your symmetry %s\n',symmetry)
end


mip = OPEN_IMG('single',mip_name);
% Search peaks first
[d1,d2] = size(mip);
ox = floor(d1/2)+1;
oy = floor(d2/2)+1;


com_radius = 2;
gaussSNR=sqrt(2).*erfcinv(2./(d1*d2*nAngles*nDefocus))
thr_radius = max(1,ceil(thr*radius/pixelSize)-1);

fprintf('\nEstimating a gaussian threshold of %2.2f\n',gaussSNR);
fprintf('\nUsing a true peak check/erase radius of %d/%d pixels\n',thr_radius,com_radius);

%mip_name = 'scaled_mip.mrc';
%gt_name = '../../complex_bgal.pdb'

system(sprintf('awk ''{if($1 == "REMARK" && $2==351) print ($3),($4)}'' %s > pts.txt',gt_name))


gt = load('pts.txt');

nPeaks = size(gt,1);

snr = 0.5:0.1:10
nVals = length(snr);

% SNR TP FN FP
roc = zeros(nVals,5);
peakList = zeros(nPeaks,1);


mip = gpuArray(mip);



% Grid for COM and interpolation of peaks

[bx, by] = BH_multi_gridCoordinates(2.*[com_radius,com_radius]+1,'Cartesian','cpu',{'none'},0,1,0);

for iPeak = 1:nPeaks
  
  % First scan the area for a maximum
  dx = gt(iPeak,1)./pixelSize + ox;
  dy = gt(iPeak,2)./pixelSize + oy;
  
  ix = floor(dx);
  iy = floor(dy);
  
  % TODO out of bounds heck
  cSq = mip(ix-thr_radius:ix+thr_radius, iy - thr_radius:iy+thr_radius);
  [m,c] = max(cSq(:));
  
  [i,j] = ind2sub(2.*[thr_radius,thr_radius]+1,c);
  i = i - thr_radius - 1;
  j = j - thr_radius - 1;
  
  comSq = gather(mip(ix-com_radius + i:ix+com_radius + i, iy - com_radius + j:iy+com_radius+ j));
  
  com = [sum(bx(:).*comSq(:)),sum(by(:).*comSq(:))]./sum(comSq(:));
  
  m2 = interpn(bx,by,comSq,com(1),com(2),'cubic');
  fprintf('found a max of %3.3f, %3.3f ,for peak %d\n',m2,m2/m,iPeak);
  
  peakList(iPeak) = m2;
  
  mip(ix-com_radius + i:ix+com_radius + i, iy - com_radius + j:iy+com_radius+ j) = 0;
  
end

eraseNeighbors=0;

for iVal = 1:nVals
  
  
  if (eraseNeighbors)
    
    kernel = ones([3,3],'single','gpuArray');
    % Remove all values that have more than two neighbors
    m1 = mip >= snr(iVal);
    for iN = 8:-1:2
      m2 = convn((m1),kernel,'same');
      m1 = m1 .*  (m1 < iN);
    end
    
    nSum = sum(m1(:));
    % Now add in the central pixel and take diff
    kernel(5) = 1;
    FP = sum(sum(convn(single(m1),kernel,'same'))) - nSum;
    
    fprintf('snr %2.2f, m1 %d, m2 %d, FP %d \n', snr(iVal),sum(mip(:) > snr(iVal)), sum(m1(:)), FP);
  else
    FP = sum(mip(:) > snr(iVal));
  end
  % True peaks are already zeroed)
  TP = sum(peakList >= snr(iVal));
  FN = nPeaks - TP;
  TN = d1*d2 - TP - FN - FP;
  
  MatthewsCCC = (TP * TN - FP * FN) ./ sqrt((TP+FP).*(TP+FN).*(TN+FP).*(TN+FN));
  
  
  RECALL = TP ./ nPeaks;
  PRECISION = TP ./ (TP + FP);
  F1 = 2. * (RECALL*PRECISION)/(RECALL+PRECISION);
  
  roc(iVal,:) = gather( ...
    [ snr(iVal), ...
    RECALL, ...
    PRECISION, ...
    F1,...
    MatthewsCCC]);
  
end

fout = fopen(sprintf('%s_roc.txt',mip_name),'w');
fprintf(fout,'%3.3f %3.3f %3.3f %3.3f %3.3f\n',roc');
fclose(fout);

onefalsePos = find(roc(:,3) < (nPeaks-1)/nPeaks,1,'last');

% Assuming a 2.5/1.25 degree search ~ 1.6e6 angles

[~,cM] = max(roc(:,4));

a = importdata('histogram.txt',' ',2);
a = a.data;
snrMax = max(find(a(:,3)>0,1,'last'),find(a(:,4)>0,1,'last')) + 1;
snrMax = a(snrMax,1);

figure('visible','off'),

subplot(1,3,[1:2]);
plot(roc(:,1),roc(:,2),'r',roc(:,1),roc(:,3),'b', ...
  snr(cM),roc(cM,5),'c*',...
  ones(length(0:0.1:1),1).*roc(onefalsePos,1),0:0.1:1,'k--');
title({sprintf('1 FP @ %2.2f (%2.2f calc) SNR and %2.2f recall\nmCCC,rec,prec,snr\n [%2.2f %2.2f %2.2f %2.2f]',roc(onefalsePos,[1]),gaussSNR,roc(onefalsePos,[2]),roc(cM,5),snr(cM),roc(cM,[2:3]))});
xlabel('SNR'); ylabel('Recall (red) Precision (blue)');

subplot(1,3,3);
plot(a(:,1),log(a(:,3)),'k--',a(:,1),log(a(:,4)),'b');
title({'Survival histogram'})
xlabel('SNR'); ylabel('log(counts)');
xlim([0,snrMax]); ylim([0.9,inf]);


saveas(gcf,sprintf('%s_roc.pdf',mip_name),'pdf');

end







