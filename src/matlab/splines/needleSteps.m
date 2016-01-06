% Find image intensity steps along the given needles in the image img.
% The columns of the array p are the points on the spline curve to which
% the needles are attached. Intensity steps are dark to bright going from
% the inside of the region outwards.

function step = needleSteps(needle, img, p)

% Interior and exterior needle endpoints
a = needle.in;
b = needle.out;

% Size of an array that can hold the image intensities along all of the n needles
delta = b - a;
len = sqrt(sum(delta .^ 2, 1));
m = ceil(max(len));
n = size(p, 2);

% Pad the image to make sure we do not hit its boundaries
img = padarray(img, [m m], 'replicate');
a = a + m;

% Points len away from a in the direction of b
c = a + ([1; 1] * ((m - 1) ./ len)) .* delta;

% Interpolation coordinates
xi = interp1([0 m-1], [a(1, :); c(1, :)], 0:(m-1));
yi = interp1([0 m-1], [a(2, :); c(2, :)], 0:(m-1));

% Image values
v = interp2(double(img), xi, yi);

% Differentiate in the vertical direction
sigma = 0.5;
[~, dk] = diffKernels(sigma);
response = conv2(v, dk, 'same');

% Threshold so that a fraction frac of points in response are above it
frac = 0.3;
threshold = prctile(response(:), 100 * (1 - frac));

% Strong local maxima
middle = response(2:(end-1), :);
mx = (response(1:(end-2), :) < middle) & (middle > response(3:end, :)) & ...
     middle >= threshold;
z = false(1, n);
mx = [z; mx; z];

% Image coordinates of local maxima
ind = mx(:);
xs = Inf(m, n);
ys = Inf(m, n);
xs(ind) = xi(ind);
ys(ind) = yi(ind);

% Pick the maximum closest to p for each needle
o = ones(m, 1);
dist2 = (xs - o * p(1, :)) .^ 2 + (ys - o * p(2, :)) .^ 2;
adist2 = (xs - o * a(1, :)) .^ 2 + (ys - o * a(2, :)) .^ 2;
dist2(adist2 > (o * (len .^ 2))) = Inf;
[dist2, row] = min(dist2);
found = dist2 < Inf;
ind = (0:(n-1)) * m + row;
ind = ind(found);

% Find maxima to subpixel resolution by fitting a parabola to the three
% response values adjacent to the pixel-resolution maxima.
% The matrix Q has rows [u_i^2, u_i 1] for u_i = -1, 0 1.
Q = [1 -1 1; 0 0 1; 1 1 1];
R = zeros(3, length(ind));
R(1, :) = response(ind - 1);
R(2, :) = response(ind);
R(3, :) = response(ind + 1);
A = Q \ R;

% Find the values of the parameter u corresponding to maxima of the
% parabolas. Construction of mx ensures that stationary points are strict
% maxima.
u = - A(2, :) ./ A(1, :) / 2;

% Find the image coordinates of the maxima by linear interpolation
step = NaN(2, n);
step(1, found) = interpolate(xi, ind, u);
step(2, found) = interpolate(yi, ind, u);

% Remove image padding offset
step = step - m;

    function xxi = interpolate(xx, ii, uu)
        aa = xx(ii);
        bb = xx(ii + 1);
        xxi = aa + (bb - aa) .* uu;
    end

end
