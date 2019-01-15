%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Lab 3: The geometry of two views 
% (application: photo-sequencing)
%
% Overleaf: https://www.overleaf.com/3753518638ggsrnmfxctzd

addpath('sift'); % ToDo: change 'sift' to the correct path where you have the sift functions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. Compute the fundamental matrix

% Two camera matrices for testing purposes
P1 = eye(3,4);
c = cosd(15); s = sind(15);
R = [c -s 0; s c 0; 0 0 1];
t = [.3 0.1 0.2]';
P2 = [R t];
n = 8;
X = [rand(3,n); ones(1,n)] + [zeros(2,n); 3 * ones(1,n); zeros(1,n)];
x1_test = P1 * X;
x2_test = P2 * X;

% Estimated fundamental matrix
% ToDo: create the following function that estimates F using the normalised 8 point algorithm
F_es = fundamental_matrix(x1_test, x2_test);

% Real fundamental matrix       % Zimmerman: Page 254
T=[0 -t(3) t(2) ; t(3) 0 -t(1) ; -t(2) t(1) 0 ];
F_gt = T*R;  % ToDo: write the expression of the real fundamental matrix for P1 and P2

% Evaluation: these two matrices should be very similar
F_gt = F_gt / norm(F_gt)
F_es = F_es / norm(F_es)

if sum(sum(abs(F_gt-F_es))) < 0.011
    'F_gt and F_es are similar'
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. Robustly fit fundamental matrix

% Read images
im1rgb = imread('Data/0000_s.png');
im2rgb = imread('Data/0001_s.png');
im1 = sum(double(im1rgb), 3) / 3 / 255;
im2 = sum(double(im2rgb), 3) / 3 / 255;

% show images
figure;
subplot(1,2,1); imshow(im1rgb); axis image; title('Image 1');
subplot(1,2,2); imshow(im2rgb); axis image; title('Image 2');

%% Compute SIFT keypoints

% (make sure that the sift folder provided in lab2 is on the path)

[points_1, desc_1] = sift(im1, 'Threshold', 0.01);
[points_2, desc_2] = sift(im2, 'Threshold', 0.01);

%% Match SIFT keypoints between a and b
close all
matches = siftmatch(desc_1, desc_2);
fig = figure;
plotmatches(im1, im2, points_1(1:2,:), points_2(1:2,:), matches, 'Stacking', 'v');
%saveas(fig,'castle_outliers.png')

% p1 and p2 contain the homogeneous coordinates of the matches
p1 = [points_1(1:2, matches(1,:)); ones(1, length(matches))];
p2 = [points_2(1:2, matches(2,:)); ones(1, length(matches))];

% ToDo: create this function (you can use as a basis 'ransac_homography_adaptive_loop.m')
[F, inliers] = ransac_fundamental_matrix(p1, p2, 2.0, 1000); 

% show inliers
fig = figure;
plotmatches(im1, im2, points_1(1:2,:), points_2(1:2,:), matches(:,inliers), 'Stacking', 'v');
title('Inliers');
%saveas(fig,'castle_inliers.png')

vgg_gui_F(im1rgb, im2rgb, F');


%% Plot some epipolar lines

l2 = F*p1; % epipolar lines in image 2 % ToDo
l1 = F'*p2; % epipolar lines in image 1 % ToDo

% choose three random indices
m1 = inliers(10);
m2 = inliers(20);
m3 = inliers(30);

% image 1 (plot the three points and their corresponding epipolar lines)
fig = figure;
imshow(im1rgb);
hold on;
plot(p1(1, m1), p1(2, m1), '+g');
plot_homog_line(l1(:, m1));

plot(p1(1, m2), p1(2, m2), '+g');
plot_homog_line(l1(:, m2));

plot(p1(1, m3), p1(2, m3), '+g');
plot_homog_line(l1(:, m3));
%saveas(fig,'castle_epipolar_lines_1.png')

% image 2 (plot the three points and their corresponding epipolar lines)
fig = figure;
imshow(im2rgb);
hold on;
plot(p2(1, m1), p2(2, m1), '+g');
plot_homog_line(l2(:, m1));

plot(p2(1, m2), p2(2, m2), '+g');
plot_homog_line(l2(:, m2));

plot(p2(1, m3), p2(2, m3), '+g');
plot_homog_line(l2(:, m3));
%saveas(fig,'castle_epipolar_lines_2.png')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. Photo-sequencing with aerial images

% In this part we will compute a simplified version of the algorithm
% explained in the Photo-sequencing paper. 
% Since we do not have two images
% taken from roughly the same viewpoint at two different time instants we
% will manually pick a dynamic point corresponding to a point in a van 
% (identified by index 'idx_car_I1') and the projection of its 3D trajectory 
% in the reference image. Then we will compute the projection (to the reference image) 
% of three points on this 3D trajectory at three different time instants 
% (corresponding to the time when the three other provided images where taken). 

clear all;

% Read images
im1rgb = imread('Data/frame_00000.tif');
im2rgb = imread('Data/frame_00001.tif');
im3rgb = imread('Data/frame_00002.tif');
im4rgb = imread('Data/frame_00003.tif');

im1 = sum(double(im1rgb), 3) / 3 / 255;
im2 = sum(double(im2rgb), 3) / 3 / 255;
im3 = sum(double(im3rgb), 3) / 3 / 255;
im4 = sum(double(im4rgb), 3) / 3 / 255;

% show images
figure;
subplot(2,2,1); imshow(im1rgb); axis image; title('Image 1');
subplot(2,2,2); imshow(im2rgb); axis image; title('Image 2');
subplot(2,2,3); imshow(im3rgb); axis image; title('Image 3');
subplot(2,2,4); imshow(im4rgb); axis image; title('Image 4');

% Compute SIFT keypoints
THRESHOLD = 0.015;          % Do not change this threshold!
[points_1, desc_1] = sift(im1, 'Threshold', THRESHOLD); 
[points_2, desc_2] = sift(im2, 'Threshold', THRESHOLD);
[points_3, desc_3] = sift(im3, 'Threshold', THRESHOLD);
[points_4, desc_4] = sift(im4, 'Threshold', THRESHOLD);

%% Take image im1 as reference image (image 1) and compute the fundamental 
% matrices needed for computing the trajectory of point idx_car_I1
% (use the SIFT keypoints previously computed)

% Compute sift matches
matches_2 = siftmatch(desc_1, desc_2);
matches_3 = siftmatch(desc_1, desc_3);
matches_4 = siftmatch(desc_1, desc_4);

% F matrix of image 2
p1 = [points_1(1:2, matches_2(1,:)); ones(1, length(matches_2))];
p2 = [points_2(1:2, matches_2(2,:)); ones(1, length(matches_2))];
[F2, inliers_2] = ransac_fundamental_matrix(p1, p2, 2.0, 1000);
'F for images 1 and 2 calculated'

% F matrix of image 3
p1 = [points_1(1:2, matches_3(1,:)); ones(1, length(matches_3))];
p2 = [points_3(1:2, matches_3(2,:)); ones(1, length(matches_3))];
[F3, inliers_3] = ransac_fundamental_matrix(p1, p2, 2.0, 1000);
'F for images 1 and 3 calculated'

% F matrix of image 4
p1 = [points_1(1:2, matches_4(1,:)); ones(1, length(matches_4))];
p2 = [points_4(1:2, matches_4(2,:)); ones(1, length(matches_4))];
[F4, inliers_4] = ransac_fundamental_matrix(p1, p2, 2.0, 1000);
'F for images 1 and 4 calculated'

%% Plot the car trajectory (keypoint idx_car_I1 in image 1)

% ToDo: complete the code

% identify the corresponding point of idx_car_I1 in the other images
idx_car_I1 = 1197;
idx_car_I2 = matches_2(2, (find(matches_2(1,:) == idx_car_I1)));
idx_car_I3 = matches_3(2, (find(matches_3(1,:) == idx_car_I1)));
idx_car_I4 = matches_4(2, (find(matches_4(1,:) == idx_car_I1)));

% coordinates (in image 1) of the keypoint idx_car_I1 (point in a van). 
% point1_1 is the projection of a 3D point in the 3D trajectory of the van
point1_1 = [points_1(1:2,idx_car_I1)' 1]';
% coordinates (in image 1) of another 3D point in the same 3D trajectory of
% the van
point1_2 = [334 697 1]'; % (this is a given data)

% l1 is the projection of the 3D trajectory of keypoint idx_car_I1
% line that joins point1_1 and point1_2
l1 = cross(point1_1, point1_2); % ToDo: compute the line

% plot the line
fig = figure;imshow(im1);
hold on;
t=1:0.1:1000;
plot(t, -(l1(1)*t + l1(3)) / l1(2), 'y');
plot(points_1(1,1197), points_1(2,1197), 'y*');

% ToDo: write the homogeneous coordinates of the corresponding point of idx_car_I1 in image 2
point2 = [points_2(1:2,idx_car_I2)' 1]';
% ToDo: compute the epipolar line of point2 in the reference image
l2 = F2'*point2;
% plot the epipolar line
plot(t, -(l2(1)*t + l2(3)) / l2(2), 'c');
% ToDo: compute the projection of point idx_car_I2 in the reference image 
pi2 = cross(l1,l2);
% plot this point
plot(pi2(1)/pi2(3), pi2(2)/pi2(3), 'c*');

% ToDo: write the homogeneous coordinates of the corresponding point of idx_car_I1 in image 3
point3 = [points_3(1:2,idx_car_I3)' 1]';
% ToDo: compute the epipolar line of point3 in the reference image
l3 = F3'*point3;
% plot the epipolar line
plot(t, -(l3(1)*t + l3(3)) / l3(2), 'r');
% ToDo: compute the projection of point idx_car_I3 in the reference image
pi3 = cross(l1,l3);
plot(pi3(1)/pi3(3), pi3(2)/pi3(3), 'r*');

% ToDo: write the homogeneous coordinates of the corresponding point of idx_car_I1 in image 4
point4 = [points_4(1:2,idx_car_I4)' 1]';
% ToDo: compute the epipolar line of point4 in the reference image
l4 = F4'*point4;
% plot the epipolar line
plot(t, -(l4(1)*t + l4(3)) / l4(2), 'g');
% ToDo: compute the projection of point idx_car_I4 in the reference image
pi4 = cross(l1,l4);
plot(pi4(1)/pi4(3), pi4(2)/pi4(3), 'g*');
hold off
%saveas(fig,'city_solution.png')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4. OPTIONAL: Photo-sequencing with your own images

% Read images
im1rgb = imread('Montserrat/IMG_2222_scaled.jpg');
im2rgb = imread('Montserrat/IMG_2223_scaled.jpg');
im3rgb = imread('Montserrat/IMG_2224_scaled.jpg');
im4rgb = imread('Montserrat/IMG_2225_scaled.jpg');
im5rgb = imread('Montserrat/IMG_2227_scaled.jpg');
im6rgb = imread('Montserrat/IMG_2228_scaled.jpg');
im7rgb = imread('Montserrat/IMG_2229_scaled.jpg');

% im1rgb = imrotate(im1rgb,180);
% im2rgb = imrotate(im2rgb,180);
% im3rgb = imrotate(im3rgb,180);
% im4rgb = imrotate(im4rgb,180);
% im5rgb = imrotate(im5rgb,180);
% im6rgb = imrotate(im6rgb,180);
% im7rgb = imrotate(im7rgb,180);

im1 = sum(double(im1rgb), 3) / 3 / 255;
im2 = sum(double(im2rgb), 3) / 3 / 255;
im3 = sum(double(im3rgb), 3) / 3 / 255;
im4 = sum(double(im4rgb), 3) / 3 / 255;
im5 = sum(double(im5rgb), 3) / 3 / 255;
im6 = sum(double(im6rgb), 3) / 3 / 255;
im7 = sum(double(im7rgb), 3) / 3 / 255;

% Compute SIFT keypoints
THRESHOLD = 0.5;          % Do not change this threshold!
[points_1, desc_1] = sift(im1, 'Threshold', THRESHOLD); 
points_1 = floor(points_1+0.5);
[points_2, desc_2] = sift(im2, 'Threshold', THRESHOLD);
points_2 = floor(points_2+0.5);
[points_3, desc_3] = sift(im3, 'Threshold', THRESHOLD);
points_3 = floor(points_3+0.5);
[points_4, desc_4] = sift(im4, 'Threshold', THRESHOLD);
points_4 = floor(points_4+0.5);
[points_5, desc_5] = sift(im5, 'Threshold', THRESHOLD); 
points_5 = floor(points_5+0.5);
[points_6, desc_6] = sift(im6, 'Threshold', THRESHOLD);
points_6 = floor(points_6+0.5);
[points_7, desc_7] = sift(im7, 'Threshold', THRESHOLD);
points_7 = floor(points_7+0.5);

matches_2 = siftmatch(desc_1, desc_2);
matches_3 = siftmatch(desc_1, desc_3);
matches_4 = siftmatch(desc_1, desc_4);
matches_5 = siftmatch(desc_1, desc_5);
matches_6 = siftmatch(desc_1, desc_6);
matches_7 = siftmatch(desc_1, desc_7);

% plotmatches(im1, im2, points_1(1:2,:), points_2(1:2,:), matches_2, 'Stacking', 'v');

% F matrix of image 2
p1 = [points_1(1:2, matches_2(1,:)); ones(1, length(matches_2))];
p2 = [points_2(1:2, matches_2(2,:)); ones(1, length(matches_2))];
[F2, inliers_2] = ransac_fundamental_matrix(p1, p2, 2.0, 1000);

% F matrix of image 3
p1 = [points_1(1:2, matches_3(1,:)); ones(1, length(matches_3))];
p2 = [points_3(1:2, matches_3(2,:)); ones(1, length(matches_3))];
[F3, inliers_3] = ransac_fundamental_matrix(p1, p2, 2.0, 1000);

% F matrix of image 4
p1 = [points_1(1:2, matches_4(1,:)); ones(1, length(matches_4))];
p2 = [points_4(1:2, matches_4(2,:)); ones(1, length(matches_4))];
[F4, inliers_4] = ransac_fundamental_matrix(p1, p2, 2.0, 1000);

% F matrix of image 2
p1 = [points_1(1:2, matches_5(1,:)); ones(1, length(matches_5))];
p2 = [points_5(1:2, matches_5(2,:)); ones(1, length(matches_5))];
[F5, inliers_5] = ransac_fundamental_matrix(p1, p2, 2.0, 1000);

% F matrix of image 3
p1 = [points_1(1:2, matches_6(1,:)); ones(1, length(matches_6))];
p2 = [points_6(1:2, matches_6(2,:)); ones(1, length(matches_6))];
[F6, inliers_6] = ransac_fundamental_matrix(p1, p2, 2.0, 1000);

% F matrix of image 4
p1 = [points_1(1:2, matches_7(1,:)); ones(1, length(matches_7))];
p2 = [points_7(1:2, matches_7(2,:)); ones(1, length(matches_7))];
[F7, inliers_7] = ransac_fundamental_matrix(p1, p2, 2.0, 1000);

%% 4.1: compute trajectory

% identify the corresponding point of idx_car_I1 in the other images
idx_car_I1 = 3332;
idx_car_I2 = matches_2(2, (find(matches_2(1,:) == idx_car_I1)));
idx_car_I3 = matches_3(2, (find(matches_3(1,:) == idx_car_I1)));
idx_car_I4 = matches_4(2, (find(matches_4(1,:) == idx_car_I1)));
idx_car_I5 = matches_5(2, (find(matches_5(1,:) == idx_car_I1)));
idx_car_I6 = matches_6(2, (find(matches_6(1,:) == idx_car_I1)));
idx_car_I7 = matches_7(2, (find(matches_7(1,:) == idx_car_I1)));

point1_1 = [points_1(1:2,idx_car_I1)' 1]';
% coordinates (in image 1) of another 3D point in the same 3D trajectory of
% the van
point1_2 = [294 266 1]'; % (this is a given data)

% l1 is the projection of the 3D trajectory of keypoint idx_car_I1
% line that joins point1_1 and point1_2
l1 = cross(point1_1, point1_2); % ToDo: compute the line

% plot the line
fig = figure;imshow(im1);
hold on;
t=1:0.1:1000;
plot(t, -(l1(1)*t + l1(3)) / l1(2), 'y');
plot(points_1(1,idx_car_I1), points_1(2,idx_car_I1), 'y*');

% ToDo: write the homogeneous coordinates of the corresponding point of idx_car_I1 in image 2
point2 = [points_2(1:2,idx_car_I2)' 1]';
% ToDo: compute the epipolar line of point2 in the reference image
l2 = F2'*point2;
% plot the epipolar line
plot(t, -(l2(1)*t + l2(3)) / l2(2), 'c');
% ToDo: compute the projection of point idx_car_I2 in the reference image 
pi2 = cross(l1,l2);
% plot this point
plot(pi2(1)/pi2(3), pi2(2)/pi2(3), 'c*');

% ToDo: write the homogeneous coordinates of the corresponding point of idx_car_I1 in image 3
point3 = [points_3(1:2,idx_car_I3)' 1]';
% ToDo: compute the epipolar line of point3 in the reference image
l3 = F3'*point3;
% plot the epipolar line
plot(t, -(l3(1)*t + l3(3)) / l3(2), 'r');
% ToDo: compute the projection of point idx_car_I3 in the reference image
pi3 = cross(l1,l3);
plot(pi3(1)/pi3(3), pi3(2)/pi3(3), 'r*');

% % ToDo: write the homogeneous coordinates of the corresponding point of idx_car_I1 in image 4
% point4 = [points_4(1:2,idx_car_I4)' 1]';
% % ToDo: compute the epipolar line of point4 in the reference image
% l4 = F4'*point4;
% % plot the epipolar line
% plot(t, -(l4(1)*t + l4(3)) / l4(2), 'g');
% % ToDo: compute the projection of point idx_car_I4 in the reference image
% pi4 = cross(l1,l4);
% plot(pi4(1)/pi4(3), pi4(2)/pi4(3), 'g*');

% ToDo: write the homogeneous coordinates of the corresponding point of idx_car_I1 in image 4
point5 = [points_5(1:2,idx_car_I5)' 1]';
% ToDo: compute the epipolar line of point4 in the reference image
l5 = F5'*point5;
% plot the epipolar line
plot(t, -(l5(1)*t + l5(3)) / l5(2), 'g');
% ToDo: compute the projection of point idx_car_I4 in the reference image
pi5 = cross(l1,l5);
plot(pi5(1)/pi5(3), pi5(2)/pi5(3), 'g*');

% ToDo: write the homogeneous coordinates of the corresponding point of idx_car_I1 in image 4
point6 = [points_6(1:2,idx_car_I6)' 1]';
% ToDo: compute the epipolar line of point4 in the reference image
l6 = F6'*point6;
% plot the epipolar line
plot(t, -(l6(1)*t + l6(3)) / l6(2), 'w');
% ToDo: compute the projection of point idx_car_I4 in the reference image
pi6 = cross(l1,l6);
plot(pi6(1)/pi6(3), pi6(2)/pi6(3), 'w*');








% [864,896] position of cable car in image 1
% [217,228] position of cable car in image 1
% [213.8,217.4] keypoint cable car in image 1 4322
% new [184,213]
% new [189,216]

% imshow(im1);
% hold on;
% plot(864,896,'r*', 'LineWidth', 2, 'MarkerSize', 10);

% 4.1 Take a set of images of a moving scene from different viewpoints at 
%     different time instants. At least two images have to be taken from
%     roughly the same location by the same camera.
%
% 4.2 Implement the first part (until line 16) of the Algorithm 1 of the 
%     Photo-sequencing paper with a selection of the detected dynamic
%     features. You may reuse the code generated for the previous question.
%
