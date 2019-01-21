function X_triang = triangulate(x1, x2, P1, P2, imsize)

nx = imsize(1);
ny = imsize(2);

H = [2/nx 0 -1;
     0  2/ny -1;
     0   0    1];

x1_h = [x1;1];
x2_h = [x2;1];

x1 = H*x1_h;
x2 = H*x2_h;

x1 = [x1(1)/x1(3); x1(2)/x1(3); 1];
x2 = [x2(1)/x2(3); x2(2)/x2(3); 1];

P1 = H*P1;
P2 = H*P2;
     
A = [x1(1)*P1(3,:) - P1(1,:);
     x1(2)*P1(3,:) - P1(2,:);
     x2(1)*P2(3,:) - P2(1,:);
     x2(2)*P2(3,:) - P2(2,:)];
 
 
 [U,S,V] = svd(A);
 
 X_triang = [V(1,4)/V(4,4); V(2,4)/V(4,4); V(3,4)/V(4,4); 1];

end