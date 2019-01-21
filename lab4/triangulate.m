function X_triang = triangulate(x1, x2, P1, P2, imsize)


A = [x1(1)*P1(3,:) - P1(1,:);
     x1(2)*P1(3,:) - P1(2,:);
     x2(1)*P2(3,:) - P2(1,:);
     x2(2)*P2(3,:) - P2(2,:);];
 
 [U,S,V] = svd(A);
 
 X_triang = [V(1,4)/V(4,4); V(2,4)/V(4,4); V(3,4)/V(4,4); 1];

end

