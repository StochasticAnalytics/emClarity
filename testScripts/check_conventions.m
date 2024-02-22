function [ ] = check_conventions()

% Make some basic data so we can test transformations and indexing options.
eve_size=12;
size_val = {[1,1].*eve_size + 1, [1,1].*eve_size};
img = {zeros(size_val{1}, 'single'), zeros(size_val{2}, 'single')};
filename_init = {'/tmp/odd.mrc', '/tmp/eve.mrc'};
filename_xf = {'/tmp/odd_rot.mrc', '/tmp/eve_rot.mrc'};

% Pick a non 1.0 pixel size to look for any issues with the indexing.
pixel_size = 1.4;

% The center of the pixel at this index is the origin of the image (zero).
% It is intended that on any transformations, an object in the continuous data remain centered at this index.
origin = cell(2,1);
for i = 1:2
    origin{i} = emc_get_origin_index(size_val{i});
    % Place a single point at X=3 and then transform the image to make sure we are where we expect to be
    img{i}(origin{i}(1) + 3, origin{i}(2)) = 1;
    SAVE_IMG(img{i}, filename_init{i}, pixel_size);
end


% First test imod transformations
% Define a rotation that rotates around Z axis by 90 degrees, this will place the point at X=3 to Y=3
% Note that if convention were fwd, or invVector, this would rotate the interpolant vector Y=3 to X=3, resulting in the data at X=3 to be at Y=3.
rot_mat = BH_defineMatrix([0,0,90], 'Bah', 'fwdVector');
% Imod uses a row major matrix, so we need to transpose the matrix
rot_mat = rot_mat';

for i = 1:2
    % Imod defines an origin that is always centered in an image, which means for
    % eve images, it is inbetween pixels, and for odd images, it is at the center of a pixel.
    % The origin may be calculated as (float(size) + 1.f)/2.f
    % emClarity uses the more standard convention of the origin being size/2 (zero based indexing, integer division).
    % So for matlab with 1 based, we calculate floor(size/2) + 1
    % The following will produce the shifts needed to shift the origin to the center of the pixel.
    origin_offset = [0.5,0.5,0.5]' * (1-mod(i, 2));
    test_shifts = rot_mat * origin_offset + origin_offset

    % NOTE: this means that if imod defines a transformation between two images (eg with tiltxcorr or tiltalign)
    % say as in the example of 10,7 -> 7,10, it would be fine if odd size, but for even sized, it would instead
    % determine that the transformation was a rotation followed by this shift, so if we want to 
    % apply this in emClarity, 

    f = fopen('/tmp/rot.xf', 'w');
    % From the newstack man page:
    %    A11 A12 A21 A22 DX DY
    %    where the coordinate (X, Y) is transformed to (X', Y') by:
    %    X' = A11 * (X - Xci) + A12 * (Y - Yci) + DX + Xco
    %    Y' = A21 * (X - Xci) + A22 * (Y - Yci) + DY + Yco
    %    where (Xci, Yci) and (Xco, Yco) are the center coordinates of the input
    %    and output images, respectively.
    fprintf(f, '%f %f %f %f %f %f\n', rot_mat([1,2,4,5]), test_shifts(1:2));
    fclose(f);

    err_msg = system(sprintf('newstack -TransformFile %s %s %s > /dev/null', '/tmp/rot.xf', filename_init{i}, filename_xf{i}));
    if (err_msg ~= 0)
        error('Failed to transform the image using newstack');
    end
end

test_img = getVolume(MRCImage('/tmp/rot.mrc'));
if any(size(test_img) ~= size(eve))
    error('The transformed image has a different size than the original');
end
% test_img = round(test_img);
if test_img(origin_eve, origin_eve + 3) ~= 1
    test_img
    error('The transformed image does not have the point at the expected location');
end
end