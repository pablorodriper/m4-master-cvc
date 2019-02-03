function [Hp] = affine_reconstruction(v1,v2,v3,v1p,v2p,v3p,Pproj1,Pproj2,w,h)

%Corresponding 3D points by triangulation
point_1 = triangulate(euclid(v1), euclid(v1p), Pproj1, Pproj2, [w, h]);
point_2 = triangulate(euclid(v2), euclid(v2p), Pproj1, Pproj2, [w, h]);
point_3 = triangulate(euclid(v3), euclid(v3p), Pproj1, Pproj2, [w, h]);

p = null([point_1'; point_2'; point_3']);

p = p/p(end);

Hp = [1 0 0 0; 
      0 1 0 0; 
      0 0 1 0; 
      p(1) p(2) p(3) p(4)];
end