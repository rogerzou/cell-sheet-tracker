function [ I ] = integrate( p, dr, dc )
%INTEGRATE 2D integration with respect to row and/or column using
%Trapezoidal rule.
%
% INPUTS
% p: (n_1 x n_2 matrix) of values.
% dr: (n_1 - 1 vector) of differences wrt row.
% dc: (n_2 - 1 vector) of differences wrt column.
%
% OUTPUTS
% I: result of integration
%
% @author Roger Zou
% @date 11/18/14

% validate inputs
narginchk(2,3);
compute_row = 0;
compute_col = 0;
if nargin == 2 && ~isempty(dr)
    compute_row = 1;
elseif nargin == 3
    if ~isempty(dr)
        compute_row = 1;
    end
    if ~isempty(dc)
        compute_col = 1;
    end
end

% computes integral using trapezoidal rule
I = p;
if compute_row
    I = trapezoidal(I, dr);
end
if compute_col
    I = trapezoidal(I', dc')';
end

%% HELPER FUNCTIONS
    function [ T ] = trapezoidal( T, dr )
    %TRAPEZOIDAL Implements trapezoidal rule
    % compute difference at each row, along each column
    T1 = [T; zeros(1, size(T,2))];
    T2 = [zeros(1, size(T,2)); T];
    T = T2 + T1;
    T = T(2:end-1, :);
    % trapezoidal rule
    T = bsxfun( @times, dr ./ 2, T );
    T = sum(T, 1);
    end

end
