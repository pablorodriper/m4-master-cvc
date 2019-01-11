function F = fundamental_matrix(x1_test, x2_test)

W = zeros(length(x1_test),length(x1_test)+1);
[norm_x1, T1] = normalise2dpts(x1_test);
[norm_x2, T2] = normalise2dpts(x2_test);

norm_x2_t = norm_x2';
norm_x1_t = norm_x1';

for i=1:length(x1_test)
    W(i,:) = [norm_x1_t(i,1)*norm_x2_t(i,1) norm_x1_t(i,2)*norm_x2_t(i,1) norm_x2_t(i,1) ...
        norm_x1_t(i,1)*norm_x2_t(i,2) norm_x1_t(i,2)*norm_x2_t(i,2) norm_x2_t(i,2) ...
        norm_x1_t(i,1) norm_x1_t(i,2) 1];
end

[U,D,V] = svd(W,0);

% H = last column of V
F_rank3 = reshape(V(:,9),3,3)';

[U2, D2, V2] = svd(F_rank3);
D2(3,3) = 0;

F = U2*D2*V2';

% Denormalization
F = T2'*F*T1;

end
        
      
        
        
