% Make one-dimensional smoothing and differentiation kernels with parameter
% sigma. To be used in differentiation functions such as grad, div, etc.

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