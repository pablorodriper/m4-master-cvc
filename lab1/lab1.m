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
I=imread('Data/0005_s.png'); % we have to be in the proper folder
% test data
%I = [1 1 1; 0 2 3;0 0 4]

sigma = pi/4;
s = 1;
R  = [cos(sigma) -sin(sigma);sin(sigma) cos(sigma)];
%R = [-1 0; 0 1];       % mirroring
H = [s.*R zeros(length(R),1); zeros(1,length(R)) ones(1,1)];

I2 = apply_H(I, H);
%I2
figure; imshow(I); figure; imshow(uint8(I2));


%% 1.2. Affinities

% ToDo: generate a matrix H which produces an affine transformation

I=imread('Data/0005_s.png'); % we have to be in the proper folder

A  = [2 2; -1 1];    
H1 = [A [1 0]'; zeros(1,length(A)) ones(1,1)];
I3 = apply_H(I, H1);
%figure; imshow(I); figure; imshow(uint8(I3));

% ToDo: decompose the affinity in four transformations: two
% rotations, a scale, and a translation
%Perform the SVD decomposition
%https://www.lucidar.me/en/mathematics/singular-value-decomposition-of-a-2x2-matrix/
[U, SIG, V] = svd2x2(A);

%Using the following reference to perform the transformations
%https://es.mathworks.com/help/images/matrix-representation-of-geometric-transformations.html
%Rotation 1
R1  = [V' zeros(length(V'),1); zeros(1,length(V')) ones(1,1)];
%Scale
S = [SIG zeros(length(SIG),1); zeros(1,length(SIG)) ones(1,1)];
%Rotation 2
R2 =[ U zeros(length( U),1); zeros(1,length( U)) ones(1,1)];
%Translation
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
if int64(I3) == int64(I4)
    disp('OK: I3 and I4 are equal')
else
    disp('ERROR: I3 and I4 are not equal')
    return
end
%figure; imshow(I); figure; imshow(uint8(I4));


%% 1.3 Projective transformations (homographies)

I=imread('Data/0005_s.png'); % we have to be in the proper folder

% ToDo: generate a matrix H which produces a projective transformation
A  = [0.5 2; -1 1];

H = [A [1 0]'; [0.1 0.1] 1];

if det(H) ~= 0
    disp('OK: The matrix H is non-singular')
else
    disp('ERROR: The matrix H is not non-singular')
    return
end

I2 = apply_H(I, H);
figure; imshow(I); figure; imshow(uint8(I2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. Affine Rectification


% choose the image points
I=imread('Data/0000_s.png');
A = load('Data/0000_s_info_lines.txt');

% indices of lines
i = 424;
p1 = [A(i,1) A(i,2) 1]';
p2 = [A(i,3) A(i,4) 1]';
i = 240;
p3 = [A(i,1) A(i,2) 1]';
p4 = [A(i,3) A(i,4) 1]';
i = 712;
p5 = [A(i,1) A(i,2) 1]';
p6 = [A(i,3) A(i,4) 1]';
i = 565;
p7 = [A(i,1) A(i,2) 1]';
p8 = [A(i,3) A(i,4) 1]';


% ToDo: compute the lines l1, l2, l3, l4, that pass through the different pairs of points

l1 = cross(p1, p2);
l2 = cross(p3, p4);
l3 = cross(p5, p6);
l4 = cross(p7, p8);


% show the chosen lines in the image
figure;imshow(I);
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

H = [1 0 0; 0 1 0; l_inf(1) l_inf(2) l_inf(3)]; %H = H_aff*[]

I2 = apply_H(I, H);
figure; imshow(uint8(I2));

% ToDo: compute the transformed lines lr1, lr2, lr3, lr4

lr1 = H.'\l1;
lr2 = H.'\l2;
lr3 = H.'\l3;
lr4 = H.'\l4;

% show the transformed lines in the transformed image
figure;imshow(uint8(I2));
hold on;
t=1:0.1:1000;
plot(t, -(lr1(1)*t + lr1(3)) / lr1(2), 'y');
plot(t, -(lr2(1)*t + lr2(3)) / lr2(2), 'y');
plot(t, -(lr3(1)*t + lr3(3)) / lr3(2), 'y');
plot(t, -(lr4(1)*t + lr4(3)) / lr4(2), 'y');

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4. Affine and Metric Rectification of the left facade of image 0001

% ToDo: Write the code that rectifies the left facade of image 0001 with
%       the stratified method (affine + metric). 
%       Crop the initial image so that only the left facade is visible.
%       Show the (properly) transformed lines that use in every step.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 5. OPTIONAL: Metric Rectification in a single step
% Use 5 pairs of orthogonal lines (pages 55-57, Hartley-Zisserman book)



