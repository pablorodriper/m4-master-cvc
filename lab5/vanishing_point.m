function [v] = vanishing_point(xo1, xf1, xo2, xf2)
%Compute lines that go through the two pair of points
l1 = cross(xo1, xf1);
l2 = cross(xo2, xf2);

%Normalize, set last coordinate to 1
l1 = l1/l1(end);
l2 = l2/l2(end);

%Compute vanishing point
v = cross(l1, l2);
end

