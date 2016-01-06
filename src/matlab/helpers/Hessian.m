% Compute the upper-triangular part of the (symmetric) Hessian of img with
% smoothing constant sigma. Also return the kernels used in the computation
% if enough output arguments are provided

function [H, skernel, dkernel] = Hessian(img, sigma)

% Assign a default standard deviation if needed
if nargin < 2 || isempty(sigma)
    sigma = 1;
end

nd = ndims(img);

% Convert to doubles if needed
img = double(img);

[g, skernel, dkernel] = grad(img, sigma);

% Allocate storage for the upper triangle of H
H = cell(nd, nd);
for i = 1:nd
    for j = i:nd
        H{i, j} = zeros(size(img));
    end
end

% Six more derivatives
for j = 1:nd
    h = grad(g{j}, sigma, j);
    for i = 1:j
        H{i, j} = h{i};
    end
end