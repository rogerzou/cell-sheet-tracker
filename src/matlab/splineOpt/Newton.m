function [x, y, X, Y] = Newton(f, df, g, x0, verbose)
%NEWTON Newton's method for optimization. Uses function f to evaluate cost
%value, function g to take Newton step.
% 
% INPUTS
% [c] = f(x): cost function to optimize
% [g, H] = df(x): gradient and Hessian of cost function
% [x, F] = g(x, step): step function, outputs new x and figure frames
% x0: domain of the optimization
% verbose: verbose output
% 
% OUTPUTS
% x: argmin/max of optimization
% y: min/max of optimization
% X: history of x at each iteration
% Y: history of y at each iteration
% 
% @authors Roger Zou and Carlo Tomasi
% @date 7/29/15

% parse input/output arguments
delta = 0.1;    % termination thresholds on newton step
epsilon = delta^2;  % termination thresholds on function value
maxIter = 30;   % max number of iterations

history = ( nargout > nargout('Newton')-2 );
if history
    X = [];
    Y = [];
end

% initial function value
x = x0;
y = f(x);

if verbose
    tStart = tic;
    iterCount = 0;
	fprintf('\n~~~~~~~~~~~~ NEWTON START ~~~~~~~~~~~~\n');
    fprintf('\nInitial function value %g at point\n', y);
end

for ii=1:maxIter

    if verbose
        fprintf('\nIteration %d start\n', ii);
    end
    
    % set y0 to y of previous iteration
    y0 = y;
    
    % compute derivatives
	[grad, Hess] = df(x);
    
	% add base regularization to avoid singular matrices
	Hess = Hess + power(eps, 1/3) * eye(size(Hess));
    
    % compute Newton step
    step = - Hess \ grad;
    
    % display function value, step, eigenvalues, and derivatives
    if verbose
        fprintf('\nNewton variables\n');
        y
        hessrcond = rcond(Hess)
        hessdensity = nnz(Hess)/numel(Hess)
        hesseig = eig(Hess);
        if length(grad) < 20
            Hess
            hesseig
            grad
            step
        else
            hesss = Hess(:);
            mmmnhess = [min(hesss), mean(hesss), max(hesss), norm(hesss)]
            mmmhesseig = [min(hesseig), mean(hesseig), max(hesseig)]
            mmmngrad = [min(grad), mean(grad), max(grad), norm(grad)]
            mmmnstep = [min(step), mean(step), max(step), norm(step)]
        end
    end 
    
    % take temporary Newton step and get new temporary function value
    xtemp = g(x, step);
    ytemp = f(xtemp);
    
    % check if function value increased
    fval_increased = false;
    if ytemp > y0
        fval_increased = true;
    end
    
    % control flow based on if function value increased
    if fval_increased
        % display function value increase info
        if verbose
            warning([sprintf('Function value increased from %g to %g\n', ...
                y0, ytemp), 'Will not take this step.\n']);
        end
        break
    else
        % update x, y
        x = xtemp;
        y = ytemp;
        % update trajectory
        if history
            X = [X x];
            Y = [Y y];
        end
        % display iteration results
        if verbose
            iterCount = iterCount + 1;
            fprintf('\nIteration %d: step', ii);
            if length(grad) < 20
                fprintf(' %10g', step);
            end
            fprintf('\nStep length %g | function value %g \n', norm(step), y);
        end
    end
    
    % check early termination conditions
    normstep = norm(step);
    absytemp = abs(ytemp - y0);
    if normstep < delta     % step size too small BREAK
        if verbose
            fprintf('NOTE: {norm(step)=%g} < {delta=%g}. Early termination.\n', normstep, delta);
        end
        break
    elseif absytemp < epsilon   % function value difference too small BREAK
        if verbose
            fprintf('NOTE: {abs(ytemp-y0)=%g} < {epsilon=%g}. Early termination.\n', absytemp, epsilon);
        end
        break
    end
    
end

% display duration
if verbose
    tEnd = toc(tStart);
    fprintf('\nDONE with Newton after %g iterations, %g seconds.\n', iterCount, tEnd);
    fprintf('\n~~~~~~~~~~~~ NEWTON END ~~~~~~~~~~~~\n');
end
