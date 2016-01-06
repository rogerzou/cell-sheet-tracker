function [ s ] = optInitSpline_A(gamma, order, k)
%OPTINITSPLINE_A Compute initial spline along a given curve with edge type A
%
% INPUTS
% gamma: (d x n) matrix of n samples in R^d along the curve.
% order: (scalar) order of the spline
% k: (scalar) number of interior control points
%
% OUTPUTS
% s: (struct) initial optimal spline that traces the curve as well as
% possible.
%
% CALLEE functions
%   splineMake
%
% @author Roger Zou
% @date 8/13/15

% convert types
gamma = double(gamma);

% get dimensions
[d, n] = size(gamma);

% test for underdetermined system, resolve by adding interpolation points
if k > (n-2)
%     warning('OPTINITSPLINE: fixing a possibly underdetermined linear system!');
    % interpolate sampling
    gamma1 = gamma(1, :);
    gamma2 = gamma(2, :);
    vq1 = interp1(1:length(gamma1), gamma1, 1:0.25:length(gamma1));
    vq2 = interp1(1:length(gamma2), gamma2, 1:0.25:length(gamma2));
    gamma = [vq1;vq2];
    % update sizes
    [d, n] = size(gamma);
end

% initialize evenly spaced spline control points
initc1 = linspace(0, 1, k+2);
initc2 = zeros(1, length(initc1));
initcontrol = [initc1; initc2];

% Multiplicities of the internal knots
mult = ones(1, k+2-order);
% Spline is open
open = true;
% How many points/needles along the spline
numPts = n;

% create initial spline from arbitrary control points
init_s = splineMake(initcontrol, order, mult, open, numPts, false);

% retrieve polynomial spline matrix P
P = double(init_s.P);

% get spline sampling components
a = gamma(:,1);
gmu = gamma(:,2:end-1);
b = gamma(:,end);

% decompose P
Q1 = P(1+d:end-d, 1:d);
Q2 = P(1+d:end-d, 1+d:end-d);
Q3 = P(1+d:end-d, 1+end-d:end);

% compute delta
delta = gmu(:) - ( Q1*a + Q3*b);

% compute least squares solution for mu
mu = Q2 \ delta;

% obtain optimal spline points
q_hat = [a; mu ; b];
newcontrol = reshape(q_hat, size(initcontrol));

% construct final spline with optimal control points
s = splineMake(newcontrol, order, mult, open, numPts, false);

end
