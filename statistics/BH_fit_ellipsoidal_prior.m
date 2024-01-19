function [ normal_vect, chi2] = BH_fit_ellipsoidal_prior(particle_coords, ...
  particle_radius_Z, ...
  radial_shrink_factor, ...
  display_fit)
% particle_coords = load('gag.txt').*10.*1.33; % input as Angstrom
%
% particle_radius_Z = 50; % Angstrom
box_size = 1.*256.*[1,1,1];

% Get the intial fit, use this to check for outliers from neighboring
% virions.
[ center, ~, ~, ~, ~] = ellipsoid_fit(particle_coords,'');

% The assumption here is that outliers will leave the center mostly
% unaffected, even if the are heavily skewed to one side, since this will
% just stretch the ellipsoid in that direction. If this assumption is true,
% then there should be at least two clusters of radius from the center
ad = sqrt(sum((particle_coords - center').^2,2));
regularizer = std(ad)/mean(ad)/10;
try
  
  gFit = fitgmdist(ad,2,'RegularizationValue',regularizer,'Start','plus',...
    'Replicates',20,'Options', statset('UseParallel', true, ...
    'MaxIter', 1000));
  save('g.mat','gFit','ad');
  
catch
  fprintf('Failed to fit, increasing regularizer\n');
  gFit = fitgmdist(ad,2,'RegularizationValue',regularizer*10);
end

gPosterior = posterior(gFit,ad);

if (display_fit)
  gFit.ComponentProportion
  save('gmFit.mat','gFit');
  figure, histogram(ad)
  hold on
  histogram(ad(gPosterior(:,1) > 1- gFit.ComponentProportion(1)))
  histogram(ad(gPosterior(:,2) > 1- gFit.ComponentProportion(2)))
end

if (abs(diff(gFit.mu)) < 2*particle_radius_Z)
  fprintf('The best mixture model shows a difference (%3.3f) which is less than particle diameter (%3.3f), assuming no outliers\n',abs(diff(gFit.mu)) , 2*particle_radius_Z);
  ignore_outliers = false(size(particle_coords,1),1);
else
  
  [~,minClass] = min(sum(gPosterior));
  gFit.ComponentProportion(minClass)
  ignore_outliers = gPosterior(:, minClass) > 1- gFit.ComponentProportion(minClass);
  fprintf('ignoring %d outliers\n',sum(ignore_outliers));
end


% Get "robust" fit ignoring outliers
[ center, ~, ~, v, chi2] = ellipsoid_fit(particle_coords(~ignore_outliers,:),'');

% Gather a set of points on this ellipsoid
mind = min( particle_coords );
maxd = max( particle_coords );
inc = max(((maxd - mind) + 1)./box_size).*[1,1,1];

[gX, gY, gZ, vX, vY, vZ] = EMC_coordGrids('cartesian', box_size, 'gpu', {'shift',-1.*single(center' ./ inc)});

gX = (gX ) .* inc(1) ;
gY = (gY ) .* inc(2) ;
gZ = (gZ ) .* inc(3) ;
vX = (vX ) .* inc(1) ;
vY = (vY ) .* inc(2) ;
vZ = (vZ ) .* inc(3) ;

% ellipse as a binary mask around -v(10)
Ellipsoid = v(1) *gX.*gX +   v(2) * gY.*gY + v(3) * gZ.*gZ + ...
  2*v(4) *gX.*gY + 2*v(5)*gX.*gZ + 2*v(6) * gY.*gZ + ...
  2*v(7) *gX    + 2*v(8)*gY    + 2*v(9) * gZ;


Ellipsoid_mask = abs(Ellipsoid +v(10)) < 1e-3;
Ellipsoid = Ellipsoid .* Ellipsoid_mask;


Ellipsoid = 1-abs(Ellipsoid - 1);
pointVal = -1.*sum(Ellipsoid(Ellipsoid_mask));
% clear grids
clear gX gY gZ

% gaussian kernel large enough to find particles, base on multiple of
% particle radius

[bhF, searchKernel] = fourierTransformer(gpuArray(BH_multi_gaussian3d(box_size,radial_shrink_factor .* particle_radius_Z / 3)));
searchKernel = conj( bhF.swapPhase(searchKernel,'fwd') );


if (display_fit)
  figure, imshow3D(gather(Ellipsoid)); hold on
  title(sprintf('fit wich chi2 %3.3e', chi2));
end
tic


normal_vect = particle_coords.*0;
xyz_vect = normal_vect;

for iPt = 1:size(particle_coords,1)
  [~, i] = min(abs(vX - particle_coords(iPt,1)));
  [~, j] = min(abs(vY - particle_coords(iPt,2)));
  [~, k] = min(abs(vZ - particle_coords(iPt,3)));
  
  ellipsoid_copy = Ellipsoid;
  ellipsoid_copy(i,j,k) = pointVal;
  
  ellipsoid_copy = Ellipsoid_mask .* real(bhF.invFFT(bhF.fwdFFT(ellipsoid_copy) .* searchKernel));
  
  
  
  [~, mc] = min(ellipsoid_copy(:));
  [i_min,j_min,k_min] = ind2sub(box_size,mc);
  
  
  
  % Calculate the gradient at the closest point on the ellipse
  % I've ommitted a factor of 2 b/c it is set to unit length anyway
  grad = [ v(1).*vX(i_min) + v(4).*vY(j_min) + v(5).*vZ(k_min) + v(7) , ...
    v(2).*vY(j_min) + v(4).*vX(i_min) + v(6).*vZ(k_min) + v(8) , ...
    v(3).*vZ(k_min) + v(5).*vX(i_min) + v(6).*vY(j_min) + v(9) ];
  
  normal_vect(iPt,:) = -gather(grad ./ norm(grad));
  xyz_vect(iPt,:) = gather([vX(i_min),vY(j_min),vZ(k_min)]);
  if (display_fit)
    scatter3(j,i,k,'rx');
    scatter3(j_min,i_min,k_min,'ro');
    
  end
  
end

if (display_fit)
  iPt = 1;
  figure, quiver3(xyz_vect(iPt,2),xyz_vect(iPt,1),xyz_vect(iPt,3),...
    normal_vect(iPt,2),normal_vect(iPt,1),normal_vect(iPt,3),30,'MaxHeadSize',6,'Marker','o','MarkerSize',2);
  hold on
  for iPt = 2:size(particle_coords,1)
    quiver3(xyz_vect(iPt,2),xyz_vect(iPt,1),xyz_vect(iPt,3), ...
      normal_vect(iPt,2),normal_vect(iPt,1),normal_vect(iPt,3),30,'MaxHeadSize',6,'Marker','o','MarkerSize',2);
  end
  fprintf('Chi2 is %3.3e\n',chi2);
  hold off
end
toc ./ iPt;


% single precision copy of ellipse, with particle position set to
% -1*sum(elilipse) to ensure the lowest value after convolution is the
% closest point on the ellipse

% Only search the relevant octant to speed things up
% ma = EMC_convn(gpuArray(single(t)), dilationKernel);

% mask with ellipse

% get coords
%$[i,j,k] = ind2sub(sigZe(e),cv)
