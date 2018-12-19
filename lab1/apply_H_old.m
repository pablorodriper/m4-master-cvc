function I_trans = apply_H(I,H)
    
    %Original image
    [m, n, color] = size(I);
    [X,Y] = meshgrid(1:n,1:m);
    
    %Compute size of destination image
    x_min = n;
    x_max = 0;
    y_min = m;
    y_max = 0;
    
    
    for i = [1,n]
        for j = [1,m]
            edge_pos = H*[i, j, 1]';
            edge_cart = [edge_pos(1)/edge_pos(3), edge_pos(2)/edge_pos(3)];
            if edge_cart(1) < y_min
                y_min = edge_cart(1);
            end
            if edge_cart(1) > y_max
                y_max = edge_cart(1);
            end
            if edge_cart(2) > x_max
                x_max = edge_cart(2);
            end
            if edge_cart(2) < x_min
                x_min = edge_cart(2);
            end
        end
    end
    
    
    %Destination image
    I_trans = zeros(ceil(x_max)-floor(x_min)+1, ceil(y_max)-floor(y_min)+1, color);
    %I_trans = zeros(ceil(y_max)-floor(y_min)+1, ceil(x_max)-floor(x_min)+1, color);
    
    

    %Transform coordinates in destination image to coordinates in original
    %image
    
    for k = 1:color
        pos_x=[];
        pos_y=[];
        for y_trans = floor(y_min):ceil(y_max)
            for x_trans = floor(x_min):ceil(x_max)
                pos = H\[y_trans; x_trans; 1];
                %Coordinates in original image
                pos_cart = [pos(1)/pos(3); pos(2)/pos(3)];
                pos_x =[pos_x pos_cart(1)];
                pos_y =[pos_y pos_cart(2)];

            end
        end

        %Interpolate values to obtain image
       I_trans(:, :, k) = reshape(interp2(X, Y, double(I(:,:,k)), pos_x, pos_y), [ceil(x_max)-floor(x_min)+1, ceil(y_max)-floor(y_min)+1]); 
    end
end