function I2 = apply_H(I, H)

% Original image
[row,column,color] = size(I);

% Destination image
I2_row = 2*(row+column)+1;
I2_column = 2*(row+column)+1;
I2 = zeros(I2_row, I2_column, color);

% Boundaries of destination image
max_row = 1; max_column = 1;
[min_row,min_column] = size(I2);

for i = 1:row
    for j = 1:column
        for k = 1:color
        % Calculate coordinates transformation
            x = H*[i,j,1]';
            dest_row = round((I2_row/2)+x(1))+1;            
            dest_column = round((I2_column/2)+x(2))+1;
            I2(dest_row,dest_column,k) = I(i,j,k);
        end
        
        % save coordinates to cut the final size of the image
        if dest_row < min_row
            min_row = dest_row;
        end
        if dest_row > max_row
            max_row = dest_row;
        end
        if dest_column < min_column
            min_column = dest_column;
        end
        if dest_column > max_column
            max_column = dest_column;
        end
    end
end

I2 = I2(min_row:max_row, min_column:max_column, :);
end