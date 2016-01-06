function [ sV2, c ] = optIterSpline( sV1, sV20, structA, structC, structD, verbose)
%OPTITERSPLINE Compute the optimal spline between two images
% Given two images I1, I2, where sV1 is the fixed spline corresponding to
% I1 and sV20 is the initial spline corresponding to I2, find the optimal
% sV2 that minimizes the error between the cost of sV1 in I1 and sV2 in I2.
% Optimization performed using Newton's Method.
% This method is used between each iteration in computing the overall cost.
%
% INPUTS
% sV1: (struct) the spline in image I1. This spline is fixed, and used to
% compute the optimal sV2 interior control points.
% sV20: (struct) the initial spline in image I2. The optimality condition on
% the interior control points is a measure of similarity between sV1 in I1
% and sV2 in I2.
% structA: images struct.
%   - I2: (R^2 matrix) image I2, the current image in a sequence.
% structC: parameters struct.
%   - l: (odd scalar) spread of the spline edge in # of pixels.
% structD: functions struct.
%   - estep
% verbose: verbose output
%
% OUTPUTS
% sV: (struct) the "optimal" sV spline in image I2.
% c: (scalar) cost of the optimal spline.
%
% CALLEE functions
%   integrate
%   jacobian
%   D2ecostDs2
%
% @authors Roger Zou, Carlo Tomasi
% @date 8/15/15

% retrieve and validate input variables
l = structC.l;
estep = structD.estep;
d = 2;
if ~mod(l, 2)
    error('OPTITERSPLINE: l must be odd');
end

if verbose
    fprintf('\n############ OPTITERSPLINE START ############\n');
end

% pack initial spline sV20 into a vector
X0 = sV20;

% Compute optimal interior control points with Newton's method 
[X, c] = Newton(@edgeCost, @edgeOpt, @edgeStep, X0, verbose);

% unpack X to retrieve final spline sV2 for image I2
sV2 = X;

if verbose
    fprintf('\n############ OPTITERSPLINE END ############\n');
end


    %% HELPER FUNCTIONS

    % compute edge cost of spline
    function c = edgeCost(X)
        c = ecost(sV1, X, structA, structC, structD);
    end

    % compute analytic gradient and Gramian for Taylor approximation
    function [g, Gr] = edgeOpt(X)
        [g, Gr] = D2ecostDs2(sV1, X, structA, structC, structD);
        g = g(1+d:end-d);
        Gr = Gr(1+d:end-d, 1+d:end-d);
    end

    % newton step - shift interior control points
    function X = edgeStep(X, step)
        X = estep(X, step);
    end

end
