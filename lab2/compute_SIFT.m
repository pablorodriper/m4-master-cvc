function compute_SIFT(ima, imb, imc, imargb, imbrgb, imcrgb, plot_figures)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. Compute image correspondences

%% Compute SIFT keypoints
[points_a, desc_a] = sift(ima, 'Threshold', 0.01);
[points_b, desc_b] = sift(imb, 'Threshold', 0.01);
[points_c, desc_c] = sift(imc, 'Threshold', 0.01);

if(plot_figures == 1)
    figure;
    subplot(1,3,1)
    imshow(imargb);%image(imargb)
    hold on;
    plot(points_a(1,:), points_a(2,:),'+y');
    subplot(1,3,2)
    imshow(imbrgb);%image(imbrgb);
    hold on;
    plot(points_b(1,:), points_b(2,:),'+y');
    subplot(1,3,3)
    imshow(imcrgb);%image(imcrgb);
    hold on;
    plot(points_c(1,:), points_c(2,:),'+y');
end

%% Match SIFT keypoints 

% between a and b
matches_ab = siftmatch(desc_a, desc_b);
if(plot_figures == 1)
    figure;
    plotmatches(ima, imb, points_a(1:2,:), points_b(1:2,:), matches_ab, 'Stacking', 'v');
end

% between b and c
matches_bc = siftmatch(desc_b, desc_c);
if(plot_figures == 1)
    figure;
    plotmatches(imb, imc, points_b(1:2,:), points_c(1:2,:), matches_bc, 'Stacking', 'v');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. Compute the homography (DLT algorithm) between image pairs

% Compute homography (normalized DLT) between a and b, play with the homography
th = 3;
xab_a = [points_a(1:2, matches_ab(1,:)); ones(1, length(matches_ab))];
xab_b = [points_b(1:2, matches_ab(2,:)); ones(1, length(matches_ab))];
[Hab, inliers_ab] = ransac_homography_adaptive_loop(xab_a, xab_b, th, 1000); % ToDo: complete this function

if(plot_figures == 1)
    figure;
    plotmatches(ima, imb, points_a(1:2,:), points_b(1:2,:), ...
        matches_ab(:,inliers_ab), 'Stacking', 'v');

    vgg_gui_H(imargb, imbrgb, Hab);
end

% Compute homography (normalized DLT) between b and c, play with the homography
xbc_b = [points_b(1:2, matches_bc(1,:)); ones(1, length(matches_bc))];
xbc_c = [points_c(1:2, matches_bc(2,:)); ones(1, length(matches_bc))];
[Hbc, inliers_bc] = ransac_homography_adaptive_loop(xbc_b, xbc_c, th, 1000); 

if(plot_figures == 1)
    figure;
    plotmatches(imb, imc, points_b(1:2,:), points_c(1:2,:), ...
        matches_bc(:,inliers_bc), 'Stacking', 'v');

    vgg_gui_H(imbrgb, imcrgb, Hbc);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. Build the mosaic

corners = [-400 1200 -100 650];
iwb = apply_H_v2(imbrgb, 'ToDo', corners);                  % ToDo: complete the call to the function
iwa = apply_H_v2(imargb, 'ToDo', corners);                  % ToDo: complete the call to the function
iwc = apply_H_v2(imcrgb, 'ToDo', corners);                  % ToDo: complete the call to the function

figure;
imshow(max(iwc, max(iwb, iwa)));%image(max(iwc, max(iwb, iwa)));
axis off;
title('Mosaic A-B-C');

end
