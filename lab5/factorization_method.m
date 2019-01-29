function [Pproj, Xproj] = factorization_method(Xh, x1, x2, flag)
%
% Q(p) :unknown homogeneus coordinate vectors of the 3D points
% P(i) unknown 3x4 image projection matrices
% q(i,p): measured homogeneus coordinate vectors of the image points.
%       -Projections of Q_p
% p = 1,...,n labels points n=Npoints
% i = 1,...,m labels images m=Ncam
% lambda(i,p): unknown scale factors
%
% W=lambda_ip*q_ip=P_i*Q_p
%
% If we could recover the depths (lambda)hs, we could apply a SVD
% factorization and thereby recover both 3D structure and camera motion for
% the scene.
%
% Returns:
%       Pproj: 3*Ncam x 4 matrix containing the camera matrices
%       Xproj: 4 x Npoints matrix of homogeneous coordinates of 3D points
% 
% As a convergence criterion you may compute the Euclidean
% distance (d) between data points and projected points in both images 
% and stop when (abs(d - d_old)/d) < 0.1 where d_old is the distance
% in the previous iteration.


% 1. Normalize the image coordinates, by applying similarity transformations
[norm_x1, T1] = normalise2dpts(x1);
[norm_x2, T2] = normalise2dpts(x2);


% 2.Estimate the fundamental matrices and epipoles with the method of [Har95]
F11= fundamental_matrix(x1, x1);
[U, D, V] = svd(F11);
%e proportional to last column of V
e11 = V(:,3) / V(3,3);
criteria = 1;

F21 = fundamental_matrix(x1, x2);
[U, D, V] = svd(F21');
e21 = V(:,3) / V(3,3);

% 3. Determine the scale factors (using equation (3))

Ncam=2;
Npoints=length(x1);
lambda= zeros(Ncam, Npoints);


%Initialization 
if isequal(flag,'ones')
    lambda= ones(Ncam, Npoints);
else
    lambda(1, :) = ones(1,Npoints);
    
    for j = 1:Npoints
        num1 = (x1(:, j)'*F11*cross(e11, x1(:,j)))*lambda(1, j);
        denom1 = norm(cross(e11, x1(:,j))).^2;
        lambda(1,j) = num1/denom1;
    end


     for j = 1:Npoints
        num2 = (x1(:, j)'*F21*cross(e21, x2(:,j)))*lambda(1, j);
        denom2 = norm(cross(e21, x2(:,j))).^2;
        lambda(2,j) = num2/denom2;
     end
     
end

while(criteria < 0.1)
    
    % 4. Build the rescaled measurement matrix W
    W = zeros(3*Ncam, Npoints);
    A = lambda(1,:) .* x1; 
    B= lambda(2,:) .* x2;
    W = cat(1, A, B);

    % 5. Balance W by column-wise and “triplet-of-rows”-wise scalar 
    % mutliplication

    % 6. Compute the SVD of the balanced matrix W
    [U,D,V] = svd(M);
    
    %7. From the SVD, recover projective motion and shape
    
    %8. Adapt projective motion, to account for the normalization
    %transformation of step 1
    
    d_old = d;
    d=0;
    criteria = (abs(d - d_old)/d); 
 
end

    
end

