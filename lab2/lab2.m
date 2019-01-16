%% init

addpath('Data')
addpath('sift')
plot_figures = true;    % if false, reduce the number of plotted figures
%iptsetpref('ImshowBorder','tight');
%ax = gca; outerpos = ax.OuterPosition; ti = ax.TightInset;  left = outerpos(1) + ti(1); bottom = outerpos(2) + ti(2); ax_width = outerpos(3) - ti(1) - ti(3); ax_height = outerpos(4) - ti(2) - ti(4); ax.Position = [left bottom ax_width ax_height];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Lab 2: Image mosaics

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Tasks 1, 2, 3 and 4 moved to function construct_mosaics()

%% ToDo: mosaic llanes
close all
imargb = imread('Data/llanes/llanes_a.jpg');
imbrgb = imread('Data/llanes/llanes_b.jpg');
imcrgb = imread('Data/llanes/llanes_c.jpg');
ima = sum(double(imargb), 3) / 3 / 255;
imb = sum(double(imbrgb), 3) / 3 / 255;
imc = sum(double(imcrgb), 3) / 3 / 255;

construct_mosaics(ima, imb, imc, imargb, imbrgb, imcrgb, plot_figures, 'llanes')

%% ToDo: compute the mosaic with castle_int images
close all
imargb = imread('Data/castle_int/0016_s.png');
imbrgb = imread('Data/castle_int/0015_s.png');
imcrgb = imread('Data/castle_int/0014_s.png');
ima = sum(double(imargb), 3) / 3 / 255;
imb = sum(double(imbrgb), 3) / 3 / 255;
imc = sum(double(imcrgb), 3) / 3 / 255;

construct_mosaics(ima, imb, imc, imargb, imbrgb, imcrgb, plot_figures, 'castle')

%% ToDo: compute the mosaic with aerial images set 13
close all
imargb = imread('Data/aerial/site13/frame00000.png');
imbrgb = imread('Data/aerial/site13/frame00002.png');
imcrgb = imread('Data/aerial/site13/frame00003.png');
ima = sum(double(imargb), 3) / 3 / 255;
imb = sum(double(imbrgb), 3) / 3 / 255;
imc = sum(double(imcrgb), 3) / 3 / 255;

construct_mosaics(ima, imb, imc, imargb, imbrgb, imcrgb, plot_figures, 'aerial_1')

%% ToDo: compute the mosaic with aerial images set 22
close all
imargb = double(imread('Data/aerial/site22/frame_00001.tif'));
imbrgb = double(imread('Data/aerial/site22/frame_00018.tif'));
imcrgb = double(imread('Data/aerial/site22/frame_00030.tif'));
ima = imargb;
imb = imbrgb;
imc = imcrgb;

construct_mosaics(ima, imb, imc, imargb, imbrgb, imcrgb, plot_figures, 'aerial_2')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 5. OPTIONAL: Calibration with a planar pattern
% 
% clear all;
% 
% %% Read template and images.
% T     = imread('Data/calib/template.jpg');
% I{1}  = imread('Data/calib/graffiti1.tif');
% I{2}  = imread('Data/calib/graffiti2.tif');
% I{3}  = imread('Data/calib/graffiti3.tif');
% %I{4}  = imread('Data/calib/graffiti4.tif');
% %I{5}  = imread('Data/calib/graffiti5.tif');
% Tg = sum(double(T), 3) / 3 / 255;
% Ig{1} = sum(double(I{1}), 3) / 3 / 255;
% Ig{2} = sum(double(I{2}), 3) / 3 / 255;
% Ig{3} = sum(double(I{3}), 3) / 3 / 255;
% 
% N = length(I);
% 
% %% Compute keypoints.
% fprintf('Computing sift points in template... ');
% [pointsT, descrT] = sift(Tg, 'Threshold', 0.05);
% fprintf(' done\n');
% 
% points = cell(N,1);
% descr = cell(N,1);
% for i = 1:N
%     fprintf('Computing sift points in image %d... ', i);
%     [points{i}, descr{i}] = sift(Ig{i}, 'Threshold', 0.05);
%     fprintf(' done\n');
% end
% 
% %% Match and compute homographies.
% H = cell(N,1);
% for i = 1:N
%     % Match against template descriptors.
%     fprintf('Matching image %d... ', i);
%     matches = siftmatch(descrT, descr{i});
%     fprintf('done\n');
% 
%     % Fit homography and remove outliers.
%     x1 = pointsT(1:2, matches(1, :));
%     x2 = points{i}(1:2, matches(2, :));
%     H{i} = 0;
%     [H{i}, inliers] =  ransac_homography_adaptive_loop(homog(x1), homog(x2), 3, 1000);
% 
%     % Plot inliers.
%     figure;
%     plotmatches(Tg, Ig{i}, pointsT(1:2,:), points{i}(1:2,:), matches(:, inliers));
% 
%     % Play with the homography
%     %vgg_gui_H(T, I{i}, H{i});
% end
% 
% %% Compute the Image of the Absolute Conic
% 
% w = ... % ToDo
%  
% %% Recover the camera calibration.
% 
% K = ... % ToDo
%     
% % ToDo: in the report make some comments related to the obtained internal
% %       camera parameters and also comment their relation to the image size
% 
% %% Compute camera position and orientation.
% R = cell(N,1);
% t = cell(N,1);
% P = cell(N,1);
% figure;hold;
% for i = 1:N
%     % ToDo: compute r1, r2, and t{i}
%     r1 = ...
%     r2 = ...
%     t{i} = ...
%     
%     % Solve the scale ambiguity by forcing r1 and r2 to be unit vectors.
%     s = sqrt(norm(r1) * norm(r2)) * sign(t{i}(3));
%     r1 = r1 / s;
%     r2 = r2 / s;
%     t{i} = t{i} / s;
%     R{i} = [r1, r2, cross(r1,r2)];
%     
%     % Ensure R is a rotation matrix
%     [U S V] = svd(R{i});
%     R{i} = U * eye(3) * V';
%    
%     P{i} = K * [R{i} t{i}];
%     plot_camera(P{i}, 800, 600, 200);
% end
% 
% % ToDo: in the report explain how the optical center is computed in the
% %       provided code
% 
% [ny,nx] = size(T);
% p1 = [0 0 0]';
% p2 = [nx 0 0]';
% p3 = [nx ny 0]';
% p4 = [0 ny 0]';
% % Draw planar pattern
% vgg_scatter_plot([p1 p2 p3 p4 p1], 'g');
% % Paint image texture
% surface('XData',[0 nx; 0 nx],'YData',[0 0; 0 0],'ZData',[0 0; -ny -ny],'CData',T,'FaceColor','texturemap');
% colormap(gray);
% axis equal;
% 
% %% Plot a static camera with moving calibration pattern.
% figure; hold;
% plot_camera(K * eye(3,4), 800, 600, 200);
% % ToDo: complete the call to the following function with the proper
% %       coordinates of the image corners in the new reference system
% for i = 1:N
%     vgg_scatter_plot( [...   ...   ...   ...   ...], 'r');
% end
% 
% %% Augmented reality: Plot some 3D points on every camera.
% [Th, Tw] = size(Tg);
% cube = [0 0 0; 1 0 0; 1 0 0; 1 1 0; 1 1 0; 0 1 0; 0 1 0; 0 0 0; 0 0 1; 1 0 1; 1 0 1; 1 1 1; 1 1 1; 0 1 1; 0 1 1; 0 0 1; 0 0 0; 1 0 0; 1 0 0; 1 0 1; 1 0 1; 0 0 1; 0 0 1; 0 0 0; 0 1 0; 1 1 0; 1 1 0; 1 1 1; 1 1 1; 0 1 1; 0 1 1; 0 1 0; 0 0 0; 0 1 0; 0 1 0; 0 1 1; 0 1 1; 0 0 1; 0 0 1; 0 0 0; 1 0 0; 1 1 0; 1 1 0; 1 1 1; 1 1 1; 1 0 1; 1 0 1; 1 0 0 ]';
% 
% X = (cube - .5) * Tw / 4 + repmat([Tw / 2; Th / 2; -Tw / 8], 1, length(cube));
% 
% for i = 1:N
%     figure; colormap(gray);
%     imagesc(Ig{i});
%     hold on;
%     x = euclid(P{i} * homog(X));
%     vgg_scatter_plot(x, 'g');
% end
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% 6. OPTIONAL: Detect the UPF logo in the two UPF images using the 
% %%              DLT algorithm (folder "logos").
% %%              Interpret and comment the results.
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% 7. OPTIONAL: Replace the logo of the UPF by the master logo
% %%              in one of the previous images using the DLT algorithm.
% 
% 

