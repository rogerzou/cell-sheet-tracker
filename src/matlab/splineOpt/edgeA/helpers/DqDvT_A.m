function [ DqDaT, DqDbT ] = DqDvT_A( Gr )
%DQDVT_A Jacobian of the interior control point function wrt vertex.
%
% INPUTS
% Gr: the m x m Gramian matrix of the spline edge cost wrt independent
% variables.
% 
% OUTPUTS
% J_eta1: 2k x 2 gradient of the control point function wrt endpoint 1.
% J_eta2: 2k x 2 gradient of the control point function wrt endpoint 2.
%
% @author Roger Zou
% @date 6/12/15

% get input variables
d = 2;
m = size(Gr,1);

% isolate A sub-matrix
A = Gr(1+d:end-d, 1+d:end-d);

% isolate the two B sub-matrix (wrt to either endpoint vertex)
B1 = Gr(1+d:end-d, 1:d);
B2 = Gr(1+d:end-d, end-d+1:end);

% regularize A
lambda = power(eps, 1/3); % regularization constant
A = A + lambda * eye(size(A));

% compute A \ B1 and A \ B2
J_eta1 = - A \ B1;
J_eta2 = - A \ B2;

DqDaT = zeros(m, d);
DqDbT = zeros(m, d);
DqDaT(1:d, :) = eye(d);
DqDbT(1+end-d:end, :) = eye(d);
DqDaT(1+d:end-d, :) = J_eta1;
DqDbT(1+d:end-d, :) = J_eta2;

end

