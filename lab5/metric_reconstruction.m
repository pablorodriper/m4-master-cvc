function [Ha] = metric_reconstruction(u,v,z,Pproj,Hp)

A = [u(1)*v(1) u(1)*v(2)+u(2)*v(1) u(1)*v(3)+u(3)*v(1) u(2)*v(2) u(2)*v(3)+u(3)*v(2) u(3)*v(3);
     u(1)*z(1) u(1)*z(2)+u(2)*z(1) u(1)*z(3)+u(3)*z(1) u(2)*z(2) u(2)*z(3)+u(3)*z(2) u(3)*z(3);
     v(1)*z(1) v(1)*z(2)+v(2)*z(1) v(1)*z(3)+v(3)*z(1) v(2)*z(2) v(2)*z(3)+v(3)*z(2) v(3)*z(3);
         0              1                   0              0              0              0;
         1              0                   0              -1             0              0];

[U,D,V] = svd(A);
w = V(:,end);

w = [w(1) w(2) w(3);
     w(2) w(4) w(5);
     w(3) w(5) w(6)];
 
P = Pproj(1:3,:)*inv(Hp);
M = P(:,1:3);
 
A_2 = chol(inv(M'*w*M));

Ha = eye(4,4);
Ha(1:3,1:3) = inv(A_2);

end