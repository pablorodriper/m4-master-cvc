function Y_initial = gs_errfunction( P0, Xobs )

%Find homography from P0
Hab_r = reshape( P0(1:9), 3, 3 );

%Find points x and x'
num_points = size(Xobs,1)/2;
x = reshape(Xobs(1:num_points), [2, num_points/2]);
xp = reshape(Xobs(num_points+1:end), [2, length(Xobs(num_points+1:end))/2]);

%Find xhat and xhatp
xhat = reshape(P0(10:end), [2, length(P0(10:end))/2]);
xhat_hom = [xhat; ones(1, size(xhat, 2))]; %bring it to homogeneous coordinates

xhatp = Hab_r*xhat_hom;

%Compute error between x and xhat, and xp and xhatp

errorx = abs(x - xhat);
errorxp = abs(xp - euclid(xhatp));

Y_initial = [errorx(:); errorxp(:)];


end