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
figure; imshow(I); title('Original image 0005\_s.png');
figure; imshow(uint8(I_rotate)); title('45 degrees');
figure; imshow(uint8(I_rotate_mirror)); title('45 degrees + horizontal mirroring');


%% 1.2. Affinities
clear all
I=imread('Data/0005_s.png'); % we have to be in the proper folder

% ToDo: generate a matrix H which produces an affine transformation

A = [2 -1; 2 1];

H1 = [A [1 0]'; zeros(1,length(A)) 1];

I_affin = apply_H(I, H1);

% Plot 
%figure; imshow(I); title('Original image 0005\_s.png');
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

figure; imshow(uint8(I_affin)); title('Affine transformation from SVD');

%% 1.3 Projective transformations (homographies)
clear all
I=imread('Data/0005_s.png'); % we have to be in the proper folder

% ToDo: generate a matrix H which produces a projective transformation

H = [14.707 0.586 1; 2.707 7.242 2; 0.01 0.01 1];

if det(H) ~= 0
    disp('OK: The matrix H is non-singular')
else
    disp('ERROR: The matrix H is not non-singular')
    return
end

I_projective = apply_H(I, H);

% Plot
%figure; imshow(I); title('Original image 0005\_s.png');
figure; imshow(uint8(I_projective)); title('Projective Transformation');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. Affine Rectification
clear all

% Choose the image points
I=imread('Data/0000_s.png');
A = load('Data/0000_s_info_lines.txt');

% Indices of lines
i = 227;
p1 = [A(i,1) A(i,2) 1]';
p2 = [A(i,3) A(i,4) 1]';
i = 367;
p3 = [A(i,1) A(i,2) 1]';
p4 = [A(i,3) A(i,4) 1]';
i = 534;
p5 = [A(i,1) A(i,2) 1]';
p6 = [A(i,3) A(i,4) 1]';
i = 576;
p7 = [A(i,1) A(i,2) 1]';
p8 = [A(i,3) A(i,4) 1]';

% ToDo: compute the lines l1, l2, l3, l4, that pass through the different pairs of points
l1 = cross(p1, p2);
l2 = cross(p3, p4);
l3 = cross(p5, p6);
l4 = cross(p7, p8);


% Show the chosen lines in the image
figure; imshow(I); title('Original image 0000\_s.png');
hold on;
t=1:0.1:1000;
plot(t, -(l1(1)*t + l1(3)) / l1(2), 'y');
plot(t, -(l2(1)*t + l2(3)) / l2(2), 'y');
plot(t, -(l3(1)*t + l3(3)) / l3(2), 'y');
plot(t, -(l4(1)*t + l4(3)) / l4(2), 'y');

% ToDo: compute the homography that affinely rectifies the image

%Vanishing points provided by each pair of parallel lines
v1 = cross(l1, l2);
v2 = cross(l3, l4);

%Vanishing line
l_inf = cross(v1, v2);
l_inf = [l_inf(1)/l_inf(3) l_inf(2)/l_inf(3) 1]';

H = [1 0 0; 0 1 0; l_inf(1) l_inf(2) l_inf(3)];
[I2, min_row2, min_col2] = apply_H(I, H);
figure; imshow(uint8(I2)); title('Affine Rectification')

% ToDo: compute the transformed lines lr1, lr2, lr3, lr4
lr1 = H.'\l1;
lr2 = H.'\l2;
lr3 = H.'\l3;
lr4 = H.'\l4;

% show the transformed lines in the transformed image
figure;imshow(uint8(I2)); title('Affine Rectification with Lines')
hold on;
t=1:0.1:1000;
plot(t+1-min_col2, 1-min_row2-(lr1(1)*t + lr1(3)) / lr1(2), 'y');
plot(t+1-min_col2, 1-min_row2-(lr2(1)*t + lr2(3)) / lr2(2), 'y');
plot(t+1-min_col2, 1-min_row2-(lr3(1)*t + lr3(3)) / lr3(2), 'y');
plot(t+1-min_col2, 1-min_row2-(lr4(1)*t + lr4(3)) / lr4(2), 'y');

% ToDo: to evaluate the results, compute the angle between the different pair 
% of lines before and after the image transformation

% Angles for original lines 
%l1 and l2 are parallel in Euclidean frame
%l3 and l4 are parallel in Euclidean frame
%Compute normal vectors and angles from there

n1 = [l1(1)/l1(3), l1(2)/l1(3)];
n2 = [l2(1)/l2(3), l2(2)/l2(3)];
n3 = [l3(1)/l3(3), l3(2)/l3(3)];
n4 = [l4(1)/l4(3), l4(2)/l4(3)];

cos_1 = (n1*n3')/(sqrt(n1*n1')*sqrt(n3*n3'))
cos_2 = (n2*n4')/(sqrt(n2*n2')*sqrt(n4*n4'))
cos_3 = (n1*n2')/(sqrt(n1*n1')*sqrt(n2*n2'))
cos_4 = (n4*n3')/(sqrt(n4*n4')*sqrt(n3*n3'))

%Angles for transformed lines

nr1 = [lr1(1)/lr1(3), lr1(2)/lr1(3)];
nr2 = [lr2(1)/lr2(3), lr2(2)/lr2(3)];
nr3 = [lr3(1)/lr3(3), lr3(2)/lr3(3)];
nr4 = [lr4(1)/lr4(3), lr4(2)/lr4(3)];

cos_1_r = (nr1*nr3')/(sqrt(nr1*nr1')*sqrt(nr3*nr3'))
cos_2_r = (nr2*nr4')/(sqrt(nr2*nr2')*sqrt(nr4*nr4'))
cos_3_r = (nr1*nr2')/(sqrt(nr1*nr1')*sqrt(nr2*nr2'))
cos_4_r = (nr4*nr3')/(sqrt(nr4*nr4')*sqrt(nr3*nr3'))

% Angles between lines should not change when applying affine rectification


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. Metric Rectification

%% 3.1 Metric rectification after the affine rectification (stratified solution)

% ToDo: Metric rectification (after the affine rectification) using two non-parallel orthogonal line pairs
%       As evaluation method you can display the images (before and after
%       the metric rectification) with the chosen lines printed on it.
%       Compute also the angles between the pair of lines before and after
%       rectification.
% choose the image points

clear all
I=imread('Data/0000_s.png');
A = load('Data/0000_s_info_lines.txt');

% indices of lines
i = 227;
p1 = [A(i,1) A(i,2) 1]';
p2 = [A(i,3) A(i,4) 1]';
i = 367;
p3 = [A(i,1) A(i,2) 1]';
p4 = [A(i,3) A(i,4) 1]';
i = 534;
p5 = [A(i,1) A(i,2) 1]';
p6 = [A(i,3) A(i,4) 1]';
i = 576;
p7 = [A(i,1) A(i,2) 1]';
p8 = [A(i,3) A(i,4) 1]';


% ToDo: compute the lines l1, l2, l3, l4, that pass through the different pairs of points
l1 = cross(p1, p2);
l2 = cross(p3, p4);
l3 = cross(p5, p6);
l4 = cross(p7, p8);

%Vanishing points provided by each pair of parallel lines
v1 = cross(l1, l2);
v2 = cross(l3, l4);

%Vanishing line
l_inf = cross(v1, v2);
l_inf = [l_inf(1)/l_inf(3) l_inf(2)/l_inf(3) 1]';
Ha = [1 0 0; 0 1 0; l_inf(1) l_inf(2) l_inf(3)];

% ToDo: compute the transformed lines lr1, lr2, lr3, lr4
lr1 = Ha'\l1;
lr2 = Ha'\l2;
lr3 = Ha'\l3;
lr4 = Ha'\l4;

%1. l1 l3 be the image of two lines that are orthogonal in the world.
%MX = b
M = [lr1(1)*lr3(1) lr1(1)*lr3(2)+lr1(2)*lr3(1); lr2(1)*lr4(1) lr2(1)*lr4(2)+lr2(2)*lr4(1)];
b = [-lr1(2)*lr3(2) -lr2(2)*lr4(2)];
X = M\b';

% Set matrix S
S = [X(1) X(2); X(2) 1];

%Find K by applying Cholesky
K = chol(S, 'lower');

%Build matrix Ha_s and apply it to image and lines
Hs_a = [K [0; 0]; [0 0] 1];
Ha_s = inv(Hs_a);
[I3, min_row3, min_col3] = apply_H(I, Ha_s);

lrr1 = Ha_s'\l1;
lrr2 = Ha_s'\l2;
lrr3 = Ha_s'\l3;
lrr4 = Ha_s'\l4;

figure; imshow(uint8(I3)); title('Metric Rectification with lines')
hold on;
t=1:0.1:1000;
plot(t+1-min_col3, 1-min_row3-(lrr1(1)*t + lrr1(3)) / lrr1(2), 'y');
plot(t+1-min_col3, 1-min_row3-(lrr2(1)*t + lrr2(3)) / lrr2(2), 'y');
plot(t+1-min_col3, 1-min_row3-(lrr3(1)*t + lrr3(3)) / lrr3(2), 'y');
plot(t+1-min_col3, 1-min_row3-(lrr4(1)*t + lrr4(3)) / lrr4(2), 'y');
hold off

%Compute normal vectors and angles from there

n1 = [l1(1)/l1(3), l1(2)/l1(3)];
n2 = [l2(1)/l2(3), l2(2)/l2(3)];
n3 = [l3(1)/l3(3), l3(2)/l3(3)];
n4 = [l4(1)/l4(3), l4(2)/l4(3)];

cos_1 = (n1*n3')/(sqrt(n1*n1')*sqrt(n3*n3'))
cos_2 = (n2*n4')/(sqrt(n2*n2')*sqrt(n4*n4'))
cos_3 = (n1*n2')/(sqrt(n1*n1')*sqrt(n2*n2'))
cos_4 = (n4*n3')/(sqrt(n4*n4')*sqrt(n3*n3'))

%Angles for transformed lines

nr1 = [lrr1(1)/lrr1(3), lrr1(2)/lrr1(3)];
nr2 = [lrr2(1)/lrr2(3), lrr2(2)/lrr2(3)];
nr3 = [lrr3(1)/lrr3(3), lrr3(2)/lrr3(3)];
nr4 = [lrr4(1)/lrr4(3), lrr4(2)/lrr4(3)];

cos_1_r = (nr1*nr3')/(sqrt(nr1*nr1')*sqrt(nr3*nr3'))
cos_2_r = (nr2*nr4')/(sqrt(nr2*nr2')*sqrt(nr4*nr4'))
cos_3_r = (nr1*nr2')/(sqrt(nr1*nr1')*sqrt(nr2*nr2'))
cos_4_r = (nr4*nr3')/(sqrt(nr4*nr4')*sqrt(nr3*nr3'))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4. Affine and Metric Rectification of the left facade of image 0001
clear all
% ToDo: Write the code that rectifies the left facade of image 0001 with
%       the stratified method (affine + metric). 
%       Crop the initial image so that only the left facade is visible.
%       Show the (properly) transformed lines that use in every step.


% Indices of lines
I=imread('Data/0001_s.png');
A = load('Data/0001_s_info_lines.txt');

I = I(:,1:450,:);
figure; imshow(I); title('Original image 0001\_s resized.png');

% Indices of lines
i = 614;
p1 = [A(i,1) A(i,2) 1]';
p2 = [A(i,3) A(i,4) 1]';
i = 159;
p3 = [A(i,1) A(i,2) 1]';
p4 = [A(i,3) A(i,4) 1]';
i = 645;
p5 = [A(i,1) A(i,2) 1]';
p6 = [A(i,3) A(i,4) 1]';
i = 541;
p7 = [A(i,1) A(i,2) 1]';
p8 = [A(i,3) A(i,4) 1]';

% ToDo: compute the lines l1, l2, l3, l4, that pass through the different pairs of points
l1 = cross(p1, p2);
l2 = cross(p3, p4);
l3 = cross(p5, p6);
l4 = cross(p7, p8);

% Crossed lines
p9 = cross(l1, l3);
p10 = cross(l1, l4);
p11 = cross(l2, l3);
p12 = cross(l2, l4);
l5 = cross(p9, p12);
l6 = cross(p10, p11);

% Show the chosen lines in the image
figure;imshow(I); title('Original image 0001\_s resized with lines.png'); 
hold on;
t=1:0.1:1000;
plot(t, -(l1(1)*t + l1(3)) / l1(2), 'y');
plot(t, -(l2(1)*t + l2(3)) / l2(2), 'y');
plot(t, -(l3(1)*t + l3(3)) / l3(2), 'y');
plot(t, -(l4(1)*t + l4(3)) / l4(2), 'y');
hold off

% ToDo: compute the homography that affinely rectifies the image

% Vanishing points provided by each pair of parallel lines
v1 = cross(l1, l2);
v2 = cross(l3, l4);

% Vanishing line
l_inf = cross(v1, v2);
l_inf = [l_inf(1)/l_inf(3) l_inf(2)/l_inf(3) 1]';

H = [1 0 0; 0 1 0; l_inf(1) l_inf(2) l_inf(3)];
[I2, min_row2, min_col2] = apply_H(I, H);

% ToDo: compute the transformed lines lr1, lr2, lr3, lr4
lr1 = H.'\l1;
lr2 = H.'\l2;
lr3 = H.'\l3;
lr4 = H.'\l4;
lr5 = H.'\l5;
lr6 = H.'\l6;

% Show the transformed lines in the transformed image
figure;imshow(uint8(I2)); title('Affine Rectification with Lines')
hold on;
t=1:0.1:1000;
plot(t+1-min_col2, 1-min_row2-(lr1(1)*t + lr1(3)) / lr1(2), 'y');
plot(t+1-min_col2, 1-min_row2-(lr2(1)*t + lr2(3)) / lr2(2), 'y');
plot(t+1-min_col2, 1-min_row2-(lr3(1)*t + lr3(3)) / lr3(2), 'g');
plot(t+1-min_col2, 1-min_row2-(lr4(1)*t + lr4(3)) / lr4(2), 'g');
plot(t+1-min_col2, 1-min_row2-(lr5(1)*t + lr5(3)) / lr5(2), 'c');
plot(t+1-min_col2, 1-min_row2-(lr6(1)*t + lr6(3)) / lr6(2), 'c');
hold off

%1. l1 l3 be the image of two lines that are orthogonal in the world.
%MX = b
M = [lr1(1)*lr3(1) lr1(1)*lr3(2)+lr1(2)*lr3(1); lr5(1)*lr6(1) lr5(1)*lr6(2)+lr5(2)*lr6(1)];
b = [-lr1(2)*lr3(2) -lr5(2)*lr6(2)];
X = M\b';

% Set matrix S
S = [X(1) X(2); X(2) 1];
K = chol(S, 'lower');

Hs_a = [K [0; 0]; [0 0] 1];
Ha_s = inv(Hs_a);
[I3, min_row3, min_col3] = apply_H(I2, Ha_s);

lrr1 = Hs_a.'*lr1;
lrr2 = Hs_a.'*lr2;
lrr3 = Hs_a.'*lr3;
lrr4 = Hs_a.'*lr4;
lrr5 = Hs_a.'*lr5;
lrr6 = Hs_a.'*lr6;

figure; imshow(uint8(I3)); title('Metric Rectification with Lines')
hold on;
t=1:0.1:1000;
plot(t+1-min_col3, 1-min_row3-(lrr1(1)*t + lrr1(3)) / lrr1(2), 'y');
plot(t+1-min_col3, 1-min_row3-(lrr2(1)*t + lrr2(3)) / lrr2(2), 'y');
plot(t+1-min_col3, 1-min_row3-(lrr3(1)*t + lrr3(3)) / lrr3(2), 'g');
plot(t+1-min_col3, 1-min_row3-(lrr4(1)*t + lrr4(3)) / lrr4(2), 'g');
plot(t+1-min_col3, 1-min_row3-(lrr5(1)*t + lrr5(3)) / lrr5(2), 'c');
plot(t+1-min_col3, 1-min_row3-(lrr6(1)*t + lrr6(3)) / lrr6(2), 'c');

hold off
