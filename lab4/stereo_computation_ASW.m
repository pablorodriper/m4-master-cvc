function disparity_map = stereo_computation(I_left, I_right, min_disp, max_disp, w_size, cost)
% Function that computes the disparity between a pair of rectified
% images using a local method based on a matching cost 
% between two local windows.
% 
% The input parameters are 5:
% - left image
% - right image
% - minimum disparity
% - maximum disparity
% - window size (e.g. a value of 3 indicates a 3x3 window)
% - matching cost (the user may able to choose between SSD and NCC costs)


%We assume both images have the same size and are in grayscale
[n, m] = size(I_left);
disparity_map = zeros(n, m);

% Create distance matrix
select_cell = [w_size/2 w_size/2];
distance_matrix = zeros(w_size, w_size);
for i = 1:w_size
    for j = 1:w_size
        distance_matrix(i, j) = sqrt((i - select_cell(1))^2 + (j - select_cell(2))^2);
    end
end

%Add padding of half window size around the image, for being able to
%compute disparity on the whole image
if(mod(w_size, 2) == 1)
    pad_size_left = floor(w_size/2);
    pad_size_right = floor(w_size/2);
else
    pad_size_left = w_size/2 -1;
    pad_size_right = w_size/2;
end

I_left_pad = padarray(I_left, [pad_size_left, pad_size_left], 0, 'pre');
I_left_pad = padarray(I_left_pad, [pad_size_right, pad_size_right], 0, 'post');

I_right_pad = padarray(I_right, [pad_size_left, pad_size_left], 0, 'pre');
I_right_pad = padarray(I_right_pad, [pad_size_right, pad_size_right], 0, 'post');

for i = 1+pad_size_left:n+pad_size_left
    for j = 1+pad_size_left:m+pad_size_left
        
        % Compute window values in left image
        I_left_window = I_left_pad(i-pad_size_left:i+pad_size_right, ...
            j-pad_size_left:j+pad_size_right);
        
        %Compute interval for the possible values in the right image, in
        %the same row
        %Left to (i,j)
        min_Ir_pos_left = max(j-max_disp, pad_size_left + 1);
        max_Ir_pos_left = j - min_disp;
        
        left_pos = [min_Ir_pos_left : max_Ir_pos_left];
        
        %Right to (i,j)
        min_Ir_pos_right = j + min_disp;
        max_Ir_pos_right = min(j+max_disp, m + pad_size_left);
        
        right_pos = [min_Ir_pos_right : max_Ir_pos_right];
        
        Ir_pos_interval = [left_pos right_pos];
        
        %Set best score for each cost function
        best_distance = Inf;

        %Go through all the possible values in the right image and compute
        %the distance with the chosen cost for picking minimum one        
        for pos_index = 1:length(Ir_pos_interval)
            j_right = Ir_pos_interval(pos_index);
            I_right_window = I_right_pad(i-pad_size_left:i+pad_size_right, ...
            j_right-pad_size_left:j_right+pad_size_right);
        
            %We assume normal distribution, and therefore all weight
            %equally
            if cost == 'ASW'
                color_difference = abs(I_left_window(round(end/2),round(end/2))-I_left_window(:,:))/12;
                position_difference = distance_matrix/(w_size/2);
                weight_left = exp(-color_difference-position_difference);
                
                color_difference = abs(I_right_window(round(end/2),round(end/2))-I_right_window(:,:))/12;
                weight_right = exp(-color_difference-position_difference);
                
                error = abs(I_left_window-I_right_window);
                error(error(:,:) > max_disp) = max_disp;
                
                distance = sum(sum((weight_left.*weight_right.*error)))/sum(sum((weight_left.*weight_right)));
                if (distance < best_distance)
                    best_distance = distance;
                    j_best = j_right;
                end
            end
        end
        disparity_map(i-pad_size_left, j-pad_size_left) = abs(j-j_best);
    end
end
end

