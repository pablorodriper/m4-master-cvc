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
I=imread('Data/0005_s.png');

% ToDo: generate a matrix H which produces a similarity transformation

% Parameters
sigma = pi/4;
s = 0.8;

% Rotation matrix
R  = [cos(sigma) -sin(sigma);sin(sigma) cos(sigma)];
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