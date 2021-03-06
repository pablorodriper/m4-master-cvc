function H = homography2d(x1_points, x2_points)
% This function calculates H such as x2 = H*x1

% Normalise each set of points
[x1_norm, T1] = normalise2dpts(x1_points);
[x2_norm, T2] = normalise2dpts(x2_points);

% Assemble the 2nx9 matrices A_i
n_points = length(x1_norm);
A = zeros(2*n_points,9);    % SVD input must be 2-D

vector_zero = [0 0 0];
for n = 1:n_points
    X = x1_norm(:,n)';
    x = x2_norm(1,n); y = x2_norm(2,n); w = x2_norm(3,n);
    A(3*n-1,:) = [vector_zero  -w*X         y*X];
    A(3*n,:) = [w*X          vector_zero  -x*X];
end

% We compute the solution by SVD
[U,D,V] = svd(A);

% H = last column of V
H = reshape(V(:,9),3,3)';

% Denormalization
H = T2\H*T1;

return