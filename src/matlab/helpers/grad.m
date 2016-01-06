%   function [g, skernel, dkernel] = grad(img, sigma, trunc)
%
% g{i} is the i-th component of the gradient of img;
% sigma is the standard deviation of the Gaussian function used for
% differentiation; if left unspecified, sigma defaults to 1 pixel.
% If enough output arguments are provided, grad also returns the smoothing
% and differentiation kernels skernel and dkernel.
% If trunc < nd, only compute the first trunc components of the gradient.
% This is useful when computing Hessians (see Hessian.m).

function [g, skernel, dkernel] = grad(img, sigma, trunc)

% Assign a default standard deviation if needed
if nargin < 2 || isempty(sigma)
    sigma = 1;
end

transpose = false;
if isvector(img)
    nd = 1;
    img = img';
    transpose = true;
else
    nd = ndims(img);
end

if nargin < 3 || isempty(trunc)
    trunc = nd;
end

if trunc > nd || trunc < 1
    error('trunc must be between 1 and nd')
end

[skernel, dkernel] = diffKernels(sigma);

% Convert to doubles if needed
img = double(img);

% Output g has as many arrays of size size(img) as img has dimensions
g = cell(trunc, 1);
for i = 1:trunc
    g{i} = zeros(size(img));
end

nconv = 0;
g{1} = img;
for i = 1:nd
    % This sequence of filtering operations achieves the minimum number of
    % convolutions. Do not exchange the order of the last two convolutions!
    for j = 1:(i-1)
        if j <= trunc
            g{j} = conv1(g{j}, skernel, i);
            nconv = nconv + 1;
        end
    end
    if i < nd
        if i + 1 <= trunc
            g{i + 1} = conv1(g{i}, skernel, i);
            nconv = nconv + 1;
        end
    end
    if i <= trunc
        g{i} = conv1(g{i}, dkernel, i);
        nconv = nconv + 1;
    end
end
%fprintf(1, '%d 1-D convolutions\n', nconv);

if transpose
    g{1} = g{1}';
end

% Make one-dimensional smoothing and differentiation kernels with parameter
% sigma

    function [skernel, dkernel, d2kernel] = diffKernels(sigma)
        
        % Make the tails of the Gaussian long enough to make truncation
        % unnoticeable
        tail = ceil(3.5 * sigma);
        x = -tail:tail;
        x = x(:);
        
        % A one-dimensional Gaussian and its derivative
        skernel = gauss(x, sigma, 1);
        dkernel = gaussDeriv(x, sigma);
        d2kernel = gaussSecondDeriv(x, sigma);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % With normalized == 0, returns samples of a 1D Gaussian with the given
        % mean and standard deviation. With normalized == 1, adjusts normalization
        % so that the sum of the samples is one.
        function g = gauss(x, sigma, normalized)
            
            if nargin < 4
                normalized = 0;
            end
            
            g = exp(- ((x / sigma) .^ 2) / 2) / sqrt(2 * pi) / sigma;
            
            if normalized
                g = g / sum(g);
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Normalizes so that the inner product of the samples with a unit-slope
        % ramp centered at the origin is minus one
        function d = gaussDeriv(x, sigma)
            
            d = (-x) .* gauss(x, sigma, 0) / sigma^2;
            
            % Normalize
            ramp = -x;
            d = d / sum(ramp .* d);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Normalizes so that the inner product of the samples with a quadratic with
        % second derivative equal to one is equal to one
        function d2 = gaussSecondDeriv(x, sigma)
            
            d2 = ((x/sigma) .^ 2 - 1) .* gauss(x, sigma);
            d2 = d2 / sum(d2 .* (x/2) .^ 2);
        end
    end

    % Compute the convolution of multidimensional array in with the
    % specified one-dimensional kernel along dimension dim.
    % The argument 'left' is useful when kernel is [1, -1]. In that case,
    % left == true (the default) leads to computing the left finite difference,
    % and left == false leads to computing the right finite difference.
    %
    % Should write a faster C version that doesn't require all this data
    % reshuffling

    function out = conv1(in, kernel, dim, left)
        
        if nargin < 4 || isempty(left)
            left = true;
        end
        
        if ~isvector(kernel)
            error('kernel must be a vector')
        end
        
        kernel = kernel(:);
        ndim = ndims(in);
        
        % Bring dimension dim to the front
        order = 1:ndim;
        order([1 dim]) = [dim 1];
        pin = permute(in, order);
        
        % Flatten the input array
        psize = size(pin);
        isize = psize(1);
        restsize = prod(psize(2:end));
        inflat = reshape(pin, [isize restsize]);
        
        % Compute a flat version of the convolution
        outvalid = conv2(inflat, kernel, 'valid');
        
        % Pad with zeros to keep the same size as the input
        % (cannot use 'same' option in conv2 because we don't want spurious
        % rows at the image boundaries)
        outflat = zeros(size(inflat));
        vrows = size(outvalid, 1);
        frows = size(outflat, 1);
        if left
            rs = 1 + ceil((frows - vrows) / 2);
        else
            rs = 1 + floor((frows - vrows) / 2);
        end
        re = rs + vrows - 1;
        outflat(rs:re, :) = outvalid;
        
        % Reshape and permute the result back
        pout = reshape(outflat, psize);
        out = permute(pout, order);
    end

end