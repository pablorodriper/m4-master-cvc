function compute_SIFT(ima, imb, imc, imargb, imbrgb, imcrgb, plot_figures)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. Compute image correspondences

%% Compute SIFT keypoints
[points_a, desc_a] = sift(ima, 'Threshold', 0.01);
[points_b, desc_b] = sift(imb, 'Threshold', 0.01);
[points_c, desc_c] = sift(imc, 'Threshold', 0.01);

if(plot_figures == 0)
    figure;
    subplot(1,3,1)
    imshow(imargb);%image(imargb)
    hold on;
    %plot(points_a(1,:), points_a(2,:),'+y');
    title('Image a');
    subplot(1,3,2)
    imshow(imbrgb);%image(imbrgb);
    hold on;
    %plot(points_b(1,:), points_b(2,:),'+y');
    title('Image b');
    subplot(1,3,3)
    imshow(imcrgb);%image(imcrgb);
    hold on;
    %plot(points_c(1,:), points_c(2,:),'+y');
    title('Image c');
    hold off;
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
iwb = apply_H_v2(imbrgb, eye(size(Hab)), corners);                  % ToDo: complete the call to the function
iwa = apply_H_v2(imargb, Hab, corners);                  % ToDo: complete the call to the function
iwc = apply_H_v2(imcrgb, inv(Hbc), corners);                     % ToDo: complete the call to the function

figure;
imshow(max(iwc, max(iwb, iwa)));%image(max(iwc, max(iwb, iwa)));
axis off;
title('Mosaic A-B-C');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4. Refine the homography with the Gold Standard algorithm

% Homography ab

x = points_a(1:2, matches_ab(1,inliers_ab));  %ToDo: set the non-homogeneous point coordinates of the 
xp = points_b(1:2, matches_ab(2,inliers_ab)); %      point correspondences we will refine with the geometric method
Xobs = [ x(:) ; xp(:) ];     % The column vector of observed values (x and x')
P0 = [ Hab(:) ; x(:) ];      % The parameters or independent variables

Y_initial = gs_errfunction( P0, Xobs ); % ToDo: create this function that we need to pass to the lsqnonlin function
% NOTE: gs_errfunction should return E(X) and not the sum-of-squares E=sum(E(X).^2)) that we want to minimize. 
% (E(X) is summed and squared implicitly in the lsqnonlin algorithm.) 
err_initial = sum( sum( Y_initial.^2 ));

options = optimset('Algorithm', 'levenberg-marquardt');
P = lsqnonlin(@(t) gs_errfunction(t, Xobs), P0, [], [], options);

Hab_r = reshape( P(1:9), 3, 3 );
f = gs_errfunction( P, Xobs ); % lsqnonlin does not return f
err_final = sum( sum( f.^2 ));

% we show the geometric error before and after the refinement
fprintf(1, 'Gold standard reproj error initial %f, final %f\n', err_initial, err_final);

% See differences in the keypoint locations

% ToDo: compute the points xhat and xhatp which are the correspondences
% returned by the refinement with the Gold Standard algorithm

xhat = reshape(P0(10:end), [2, length(P0(10:end))/2]);
xhat_hom = [xhat; ones(1, size(xhat, 2))]; %bring it to homogeneous coordinates

xhatp = euclid(Hab_r*xhat_hom);

figure;
imshow(imargb);%image(imargb);
hold on;
plot(x(1,:), x(2,:),'+y');
plot(xhat(1,:), xhat(2,:),'+c');
title('Image A: original (yellow) and refined (blue) correspondences'); 
hold off

figure;
imshow(imbrgb);%image(imbrgb);
hold on;
plot(xp(1,:), xp(2,:),'+y');
plot(xhatp(1,:), xhatp(2,:),'+c');
title('Image B: original (yellow) and refined (blue) correspondences');
hold off

%%  Homography bc

% ToDo: refine the homography bc with the Gold Standard algorithm

x = points_b(1:2, matches_bc(1,inliers_bc));  %ToDo: set the non-homogeneous point coordinates of the 
xp = points_c(1:2, matches_bc(2,inliers_bc)); %      point correspondences we will refine with the geometric method
Xobs = [ x(:) ; xp(:) ];     % The column vector of observed values (x and x')
P0 = [ Hbc(:) ; x(:) ];      % The parameters or independent variables

Y_initial = gs_errfunction( P0, Xobs ); % ToDo: create this function that we need to pass to the lsqnonlin function
% NOTE: gs_errfunction should return E(X) and not the sum-of-squares E=sum(E(X).^2)) that we want to minimize. 
% (E(X) is summed and squared implicitly in the lsqnonlin algorithm.) 
err_initial = sum( sum( Y_initial.^2 ));

options = optimset('Algorithm', 'levenberg-marquardt');
P = lsqnonlin(@(t) gs_errfunction(t, Xobs), P0, [], [], options);

Hbc_r = reshape( P(1:9), 3, 3 );
f = gs_errfunction( P, Xobs ); % lsqnonlin does not return f
err_final = sum( sum( f.^2 ));

% we show the geometric error before and after the refinement
fprintf(1, 'Gold standard reproj error initial %f, final %f\n', err_initial, err_final);




%% See differences in the keypoint locations

% ToDo: compute the points xhat and xhatp which are the correspondences
% returned by the refinement with the Gold Standard algorithm

xhat = reshape(P0(10:end), [2, length(P0(10:end))/2]);
xhat_hom = [xhat; ones(1, size(xhat, 2))]; %bring it to homogeneous coordinates

xhatp = euclid(Hbc_r*xhat_hom);

figure;
imshow(imbrgb);%image(imbrgb);
hold on;
plot(x(1,:), x(2,:),'+y');
plot(xhat(1,:), xhat(2,:),'+c');
title('Image B: original (yellow) and refined (blue) correspondences'); 
hold off;

figure;
imshow(imcrgb);%image(imcrgb);
hold on;
plot(xp(1,:), xp(2,:),'+y');
plot(xhatp(1,:), xhatp(2,:),'+c');
title('Image C: original (yellow) and refined (blue) correspondences'); 
hold off

%% Build mosaic
corners = [-400 1200 -100 650];
iwb = apply_H_v2(imbrgb, eye(size(Hab_r)), corners); % ToDo: complete the call to the function
iwa = apply_H_v2(imargb, Hab_r, corners); % ToDo: complete the call to the function
iwc = apply_H_v2(imcrgb, inv(Hbc_r), corners); % ToDo: complete the call to the function

figure;
imshow(max(iwc, max(iwb, iwa)));%image(max(iwc, max(iwb, iwa)));axis off;
title('Mosaic A-B-C');

end
