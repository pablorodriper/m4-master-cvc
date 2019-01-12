function [F, idx_inliers] = ransac_fundamental_matrix(p1, p2, th, max_it)

[Ncoords, Npoints] = size(p1);

% ransac
it = 0;
best_inliers = [];
% probability that at least one random sample set is free of outliers
p = 0.999; 
while it < max_it
    points = randomsample(Npoints, 8);
    F = fundamental_matrix(p1(:,points), p2(:,points));
    inliers = compute_inliers(F, p1, p2, th);
    
    % test if it is the best model so far
    if length(inliers) > length(best_inliers)
        %length(best_inliers)
        best_inliers = inliers;
    end    
    
    % update estimate of max_it (the number of trials) to ensure we pick, 
    % with probability p, an initial data set with no outliers
    fracinliers =  length(inliers)/Npoints;
    pNoOutliers = 1 -  fracinliers^4;
    pNoOutliers = max(eps, pNoOutliers);  % avoid division by -Inf
    pNoOutliers = min(1-eps, pNoOutliers);% avoid division by 0
    max_it = log(1-p)/log(pNoOutliers);
    
    it = it + 1;
end

% compute F from all the inliers
F = fundamental_matrix(p1(:,best_inliers), p2(:,best_inliers));
idx_inliers = best_inliers;


function idx_inliers = compute_inliers(F, p1, p2, th)

    num = zeros(1,length(p1));
	for n = 1:length(p1)
	    num(n) = p2(:,n)'*F*p1(:,n);
	end
    den_1 = F*p1;
    den_2 = F'*p2;
    
    % compute the symmetric geometric error
    d2 = ((num.^2)./(den_1(1,:).^2+den_1(2,:).^2+den_2(1,:).^2+den_2(2,:).^2));
    idx_inliers = find(d2 < th.^2);


function xn = normalise(x)    
    xn = x ./ repmat(x(end,:), size(x,1), 1);

    
function item = randomsample(npts, n)
	a = [1:npts]; 
    item = zeros(1,n);    
    for i = 1:n
	  % Generate random value in the appropriate range 
	  r = ceil((npts-i+1).*rand);
	  item(i) = a(r);       % Select the rth element from the list
	  a(r)    = a(end-i+1); % Overwrite selected element
    end                       % ... and repeat

