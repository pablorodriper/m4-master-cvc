%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Lab 1: Image rectification


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. Applying image transformations

% ToDo: create the function  "apply_H" that gets as input a homography and
% an image and returns the image transformed by the homography.
% The size of the transformed image has to be automatically set so as to 
% contain the whole transformed image.
% At some point you will need to interpolate the image values at some points,
% you may use the Matlab function "interp2" for that.


%% 1.1. Similarities
clear all
I=imread('Data/0005_s.png');

% ToDo: generate a matrix H which produces a similarity transformation

% Parameters
phi = pi/4;
s = 0.8;

% Rotation matrix
R  = [cos(phi) -sin(phi);sin(phi) cos(phi)];
R_mirror = [1 0; 0 -1];  

H = [s.*R zeros(length(R),1); zeros(1,length(R)) 1];
H2 = [s.*R*R_mirror zeros(length(R),1); zeros(1,length(R)) 1];

I_rotate = apply_H(I, H);
I_rotate_mirror = apply_H(I, H2);

% Plot
figure; imshow(I); title('Original image 0005_s.png');
figure; imshow(uint8(I_rotate)); title('45 degrees');
figure; imshow(uint8(I_rotate_mirror)); title('45 degrees + horizontal mirroring');


%% 1.2. Affinities

% ToDo: generate a matrix H which produces an affine transformation

I=imread('Data/0005_s.png'); % we have to be in the proper folder

%% 1.2. Affinities

I=imread('Data/0005_s.png'); % we have to be in the proper folder

% ToDo: generate a matrix H which produces an affine transformation

A = [2 -1; 2 1]

H1 = [A [1 0]'; zeros(1,length(A)) 1];

I_affin = apply_H(I, H1);

% Plot 
%figure; imshow(I); title('Original image 0005_s.png');
figure; imshow(uint8(I_affin)); title('Affine transformation');

% ToDo: decompose the affinity in four transformations: two
% rotations, a scale, and a translation

% Perform the SVD decomposition
% https://www.lucidar.me/en/mathematics/singular-value-decomposition-of-a-2x2-matrix/
[U, SIG, V] = svd2x2(A);

% Using the following reference to perform the transformations
% https://es.mathworks.com/help/images/matrix-representation-of-geometric-transformations.html
% Rotation 1
R1 = [V' zeros(length(V'),1); zeros(1,length(V')) 1];
% Scale
S = [SIG zeros(length(SIG),1); zeros(1,length(SIG)) 1];
% Rotation 2
R2 =[U zeros(length( U),1); zeros(1,length(U)) 1];
% Translation
T = [1 0 0; 0 1 0; 1 0 1];

% ToDo: verify that the product of the four previous transformations
% produces the same matrix H as above
H2 = T'*R2*S*R1; H2 = double(int64(H2));
if H1 == H2
    disp('OK: H1 and H2 are equal')
else
    disp('ERROR: H1 and H2 are not equal')
    return
end

% ToDo: verify that the proper sequence of the four previous
% transformations over the image I produces the same image I2 as before
I4 = apply_H(I, H2);
if int64(I_affin) == int64(I4)
    disp('OK: I3 and I4 are equal')
else
    disp('ERROR: I3 and I4 are not equal')
    return
end
%figure; imshow(I); 
figure; imshow(uint8(I_affin)); title('Affine transformation from SVD');


