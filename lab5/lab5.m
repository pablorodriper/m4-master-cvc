%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Lab 5: Reconstruction from uncalibrated viewas

addpath('../lab2/sift'); % ToDo: change 'sift' to the correct path where you have the sift functions
addpath('Data');
iptsetpref('ImshowBorder','tight')  % Save images without border
addpath(genpath('vanishing_points_v0.9/'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 0. Create synthetic data

[X, P1, P2, w, h] = create_synthetic_data();

% visualize as point cloud
ax1 = subplot(1,2,1); hold on;
plot_camera2(P1,w,h);plot_camera2(P2,w,h);
for i = 1:length(X)
  scatter3(X(1,i), X(2,i), X(3,i), 5^2, [0.5 0.5 0.5], 'filled');
end
axis equal;axis vis3d;
title('Point cloud')

% visualize as lines
ax2 = subplot(1,2,2);
visualize_as_lines(X);
plot_camera2(P1,w,h);plot_camera2(P2,w,h);
title('Ground truth')
hold off
% Move both cameras at the same time in the figure
hlink = linkprop([ax1,ax2],{'CameraPosition','CameraUpVector'}); 
rotate3d on

%% Create homogeneous coordinates
% homogeneous 3D coordinates
Xh=[X; ones(1,length(X))];

% homogeneous 2D coordinates
x1 = P1*Xh;
x1(1,:) = x1(1,:)./x1(3,:);
x1(2,:) = x1(2,:)./x1(3,:);
x1(3,:) = x1(3,:)./x1(3,:);
x2 = P2*Xh;
x2(1,:) = x2(1,:)./x2(3,:);
x2(2,:) = x2(2,:)./x2(3,:);
x2(3,:) = x2(3,:)./x2(3,:);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. Projective reconstruction (synthetic data)

% ToDo: create the function 'factorization_method' that computes a
% projective reconstruction with the factorization method of Sturm and
% Triggs '1996

% if flag == 'ones' init lamda = 1
[Pproj_ones, Xproj_ones] = factorization_method(x1, x2, 'ones');
[Pproj, Xproj] = factorization_method(x1, x2, flag);


%% Check projected points (estimated and data points)
%Ones
for i=1:2
    x_proj_ones{i} = euclid(Pproj_ones(3*i-2:3*i,:)*Xproj_ones);
end
x_d_ones{1} = euclid(P1*Xh);
x_d_ones{2} = euclid(P2*Xh);

% image 1

figure; subplot(2,2,1); hold on
plot(x_d_ones{1}(1,:),x_d_ones{1}(2,:),'r*'); 
plot(x_proj_ones{1}(1,:),x_proj_ones{1}(2,:),'bo');
axis equal
title('Image 1 with ones initialization')
hold off

% image 2
subplot(2,2,2); hold on
plot(x_d_ones{2}(1,:),x_d_ones{2}(2,:),'r*');
plot(x_proj_ones{2}(1,:),x_proj_ones{2}(2,:),'bo');
title('Image 2 with ones initialization')
hold off

%Sturm
for i=1:2
    x_proj{i} = euclid(Pproj(3*i-2:3*i,:)*Xproj);
end
x_d{1} = euclid(P1*Xh);
x_d{2} = euclid(P2*Xh);

% image 1
subplot(2,2,3); hold on
plot(x_d{1}(1,:),x_d{1}(2,:),'r*'); 
plot(x_proj{1}(1,:),x_proj{1}(2,:),'bo');
axis equal
title('Image 1 with Sturm initialization')
hold off

% image 2
subplot(2,2,4); hold on
plot(x_d{2}(1,:),x_d{2}(2,:),'r*');
plot(x_proj{2}(1,:),x_proj{2}(2,:),'bo');
title('Image 2 with Sturm initialization')
hold off

% Visualize projective reconstruction
Xaux(1,:) = Xproj(1,:)./Xproj(4,:);
Xaux(2,:) = Xproj(2,:)./Xproj(4,:);
Xaux(3,:) = Xproj(3,:)./Xproj(4,:);

figure;
visualize_as_lines(Xaux);
title('Projective reconstruction')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. Affine reconstruction (synthetic data)

% ToDo: create the function 'vanishing_point' that computes the vanishing
% point formed by the line that joins points xo1 and xf1 and the line 
% that joins points x02 and xf2
%
% [v1] = vanishing_point(xo1, xf1, xo2, xf2)

% Compute the vanishing points in each image
v1 = vanishing_point(x1(:,21),x1(:,22),x1(:,23),x1(:,24));
v2 = vanishing_point(x1(:,21),x1(:,23),x1(:,22),x1(:,24));
v3 = vanishing_point(x1(:,1),x1(:,2),x1(:,4),x1(:,3));

v1p = vanishing_point(x2(:,21),x2(:,22),x2(:,23),x2(:,24));
v2p = vanishing_point(x2(:,21),x2(:,23),x2(:,22),x2(:,24));
v3p = vanishing_point(x2(:,1),x2(:,2),x2(:,4),x2(:,3));

Pproj1 = Pproj(1:3, :);
Pproj2 = Pproj(4:6, :);

% ToDo: use the vanishing points to compute the matrix Hp that 
%       upgrades the projective reconstruction to an affine reconstruction

Hp = affine_reconstruction(v1,v2,v3,v1p,v2p,v3p,Pproj1,Pproj2,w,h);

% check results
Xa = euclid(Hp*Xproj);
figure;visualize_as_lines(Xa);
title('Affine reconstruction')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. Metric reconstruction (synthetic data)

% ToDo: compute the matrix Ha that 
%       upgrades the projective reconstruction to an affine reconstruction
% Use the following vanishing points given by three pair of orthogonal lines
% and assume that the skew factor is zero and that pixels are square

u = vanishing_point(x1(:,2),x1(:,5),x1(:,3),x1(:,6));
v = vanishing_point(x1(:,1),x1(:,2),x1(:,3),x1(:,4));
z = vanishing_point(x1(:,1),x1(:,4),x1(:,2),x1(:,3));
 
Ha = metric_reconstruction(u,v,z,Pproj,Hp);

% check results
Xm = euclid(Ha*Hp*Xproj);
figure;visualize_as_lines(Xm);
title('Metric reconstruction')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4. Projective reconstruction (real data)

%% read images
Irgb{1} = double(imread('Data/0000_s.png'))/255;
Irgb{2} = double(imread('Data/0001_s.png'))/255;

I{1} = sum(Irgb{1}, 3) / 3; 
I{2} = sum(Irgb{2}, 3) / 3;

Ncam = length(I);

% Compute keypoints and matches.
points = cell(2,1);
descr = cell(2,1);
for i = 1:2
    [points{i}, descr{i}] = sift(I{i}, 'Threshold', 0.01);
    points{i} = points{i}(1:2,:);
end

matches = siftmatch(descr{1}, descr{2});

% Plot matches.
figure;
plotmatches(I{1}, I{2}, points{1}, points{2}, matches, 'Stacking', 'v');

% Fit Fundamental matrix and remove outliers.
x1 = points{1}(:, matches(1, :));
x2 = points{2}(:, matches(2, :));

[F, inliers] = ransac_fundamental_matrix(homog(x1), homog(x2), 2.0);

% Plot inliers.
inlier_matches = matches(:, inliers);
figure;
plotmatches(I{1}, I{2}, points{1}, points{2}, inlier_matches, 'Stacking', 'v');

x1 = points{1}(:, inlier_matches(1, :));
x2 = points{2}(:, inlier_matches(2, :));

% homogeneous 2D coordinates
x1 = homog(x1);
x2 = homog(x2);

% ToDo: compute a projective reconstruction using the factorization method
[Pproj, Xproj] = factorization_method(x1, x2, flag);

% ToDo: show the data points (image correspondences) and the projected
% points (of the reconstructed 3D points) in images 1 and 2. Reuse the code
% in section 'Check projected points' (synthetic experiment).

for i=1:2
    x_proj{i} = euclid(Pproj(3*i-2:3*i, :)*Xproj);
end
x_d{1} = euclid(x1);
x_d{2} = euclid(x2);

err_cam1 = sqrt(sum( x_proj{1}-x_d{1}).^2);
err_cam2 = sqrt(sum(  x_proj{2}-x_d{2}).^2);

figure;
h1 = histogram(err_cam1);
hold on
h2 = histogram(err_cam2);
mean_err1 = mean(err_cam1);
mean_err2 = mean(err_cam2);
total_mean = (mean_err1+mean_err2)/2
line([total_mean total_mean], ylim, 'Color','g');
legend('hist error camera 1', 'hist error camera 2', 'mean error value');
hold off

% image 1
figure;subplot(1,2,2);
imshow(I{1});
hold on
plot(x_d{1}(1,:),x_d{1}(2,:),'r*'); 
plot(x_proj{1}(1,:),x_proj{1}(2,:),'bo');
axis equal
hold off

% image 2
subplot(1,2,1);
imshow(I{2});
hold on
plot(x_d{2}(1,:),x_d{2}(2,:),'r*');
plot(x_proj{2}(1,:),x_proj{2}(2,:),'bo');
hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 5. Affine reconstruction (real data)

% ToDo: compute the matrix Hp that updates the projective reconstruction
% to an affine one
%
% You may use the vanishing points given by function 'detect_vps' that 
% implements the method presented in Lezama et al. CVPR 2014
% (http://dev.ipol.im/~jlezama/vanishing_points/)

% This is an example on how to obtain the vanishing points (VPs) from three
% orthogonal lines in image 1

w = size(Irgb{1},1);
h = size(Irgb{1},2);
img_in =  'Data/0000_s.png'; % input image

% folder_out = '.'; % output folder
% manhattan = 1;
% acceleration = 0;
% focal_ratio = 1;
% params.PRINT = 1;
% params.PLOT = 1;
% [horizon, VPs] = detect_vps(img_in, folder_out, manhattan, acceleration, focal_ratio, params);

% We could not compile the source code, so we load the output of detect_vps
params = load('VPs.mat');

% Compute the vanishing points in each image
v1 = [params.VPs_0(:,1); 1];
v2 = [params.VPs_0(:,2); 1];
v3 = [params.VPs_0(:,3); 1];

v1p = [params.VPs_1(:,1); 1];
v2p = [params.VPs_1(:,2); 1];
v3p = [params.VPs_1(:,3); 1];

Pproj1 = Pproj(1:3, :);
Pproj2 = Pproj(4:6, :);

Hp = affine_reconstruction(v1,v2,v3,v1p,v2p,v3p,Pproj1,Pproj2,w,h);


%% Visualize the result

% x1m are the data points in image 1
x1m = x_d{1};
% Xm are the reconstructed 3D points (projective reconstruction)
Xm = Xproj;

r = interp2(double(Irgb{1}(:,:,1)), x1m(1,:), x1m(2,:));
g = interp2(double(Irgb{1}(:,:,2)), x1m(1,:), x1m(2,:));
b = interp2(double(Irgb{1}(:,:,3)), x1m(1,:), x1m(2,:));
Xe = euclid(Hp*Xm);
figure; hold on;
[w,h] = size(I{1});
for i = 1:length(Xe)
    scatter3(Xe(1,i), Xe(2,i), Xe(3,i), 80, [r(i) g(i) b(i)], 'filled');
end
axis equal;title('Affine reconstruction')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 6. Metric reconstruction (real data)

% ToDo: compute the matrix Ha that updates the affine reconstruction
% to a metric one and visualize the result in 3D as in the previous section

u = homog(params.VPs_0(:,1));
v = homog(params.VPs_0(:,2));
z = homog(params.VPs_0(:,3));

Ha = metric_reconstruction(u,v,z,Pproj,Hp);

% check results
Xa = euclid(Ha*Hp*Xproj);

figure; hold on;
for i = 1:length(Xa)
    scatter3(Xa(1,i), -Xa(2,i), Xa(3,i), 70, [r(i) g(i) b(i)], 'filled');
end
axis equal;title('Metric reconstruction')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 7. OPTIONAL: Projective reconstruction from two views

% ToDo: compute a projective reconstruction from the same two views 
% by computing two possible projection matrices from the fundamental matrix
% and one of the epipoles.
% Then update the reconstruction to affine and metric as before (reuse the code).

% The first projection matrix is always 
P1 = [eye(3,3) zeros(3,1)];

% Calculate e' from e'^T*F = 0
[U, D, V] = svd(F); % The fundamental matrix is previously calculated
e = V(:,3) / V(3,3);
skew_e =[0 -e(3) e(2) ; e(3) 0 -e(1) ; -e(2) e(1) 0 ];

P2 = [skew_e*F e];

if (rank(P2) == 3)
    fprintf('OK: rank of P2 = 3\n');
end

%% (reuse the code)
% toDo: update the reconstruction to affine and metric as before
params = load('VPs.mat');

% Compute the vanishing points in each image
v1 = [params.VPs_0(:,1); 1];
v2 = [params.VPs_0(:,2); 1];
v3 = [params.VPs_0(:,3); 1];

v1p = [params.VPs_1(:,1); 1];
v2p = [params.VPs_1(:,2); 1];
v3p = [params.VPs_1(:,3); 1];

Hp = affine_reconstruction(v1,v2,v3,v1p,v2p,v3p,P1,P2,w,h);

% x1m are the data points in image 1
x1m = x_d{1};
% Xm are the reconstructed 3D points (projective reconstruction)

%%%%%%% Probably the error is in this line, but we are not sure on 
%%%%%%% what Xproj use.
Xm = (x2'*P2)';
%%%%%%%%
%%%%%%%%

r = interp2(double(Irgb{1}(:,:,1)), x1m(1,:), x1m(2,:));
g = interp2(double(Irgb{1}(:,:,2)), x1m(1,:), x1m(2,:));
b = interp2(double(Irgb{1}(:,:,3)), x1m(1,:), x1m(2,:));
Xe = euclid(Hp*Xm);
figure; hold on;
[w,h] = size(I{1});
for i = 1:length(Xe)
    scatter3(Xe(1,i), Xe(2,i), Xe(3,i), 80, [r(i) g(i) b(i)], 'filled');
end
axis equal;title('Affine reconstruction')

% 6. Metric reconstruction (real data)

% ToDo: compute the matrix Ha that updates the affine reconstruction
% to a metric one and visualize the result in 3D as in the previous section

u = homog(params.VPs_0(:,1));
v = homog(params.VPs_0(:,2));
z = homog(params.VPs_0(:,3));

Ha = metric_reconstruction(u,v,z,P1,Hp);

% check results
Xa = euclid(Ha*Hp*Xproj);

figure; hold on;
for i = 1:length(Xa)
    scatter3(Xa(1,i), -Xa(2,i), Xa(3,i), 70, [r(i) g(i) b(i)], 'filled');
end
axis equal;title('Metric reconstruction')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 8. OPTIONAL: Projective reconstruction from more than two views

% ToDo: extend the function that computes the projective reconstruction 
% with the factorization method to the case of three views. You may use 
% the additional image '0002_s.png'
% Then update the reconstruction to affine and metric.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 9. OPTIONAL: Any other improvement you may incorporate 

% Add a 4th view, incorporate new 3D points by triangulation, 
% incorporate new views by resectioning, 
% apply any kind of processing on the point cloud, ...)

