Copyright (c) 2020 Robotic Exploration Lab

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

function [u,result,Luu,free] = boxQPsolve(Quu,Qu,lower_lim,upper_lim,u0)
% Finds the optimal control within limits to minimize a quadratic cost
% Minimizes 0.5*u'*Quu*u + u'*Qu  s.t. lower_lim <= u <= upper_lim
%
% Inputs:
% ==========================================
% Quu       - control cost Hessian (positive definite)   (m, m)
% Qu        - control cost Jacobian                      (m)
% lower     - control lower limit                        (m)
% upper     - control upper limit                        (m)
% u0        - initial control input for warm-start       (m)
%
% Outputs:
% =====================================
% u         - optimal feed-forward control                (m)
% result    - gives exit criterion (explained below)         
% Luu       - cholesky factor                             (m, m)
% free      - set of free dimensions                      (m)

% Results
% ===========================
%  0: No descent direction found
%  1: Hessian is not positive definite
%  2: Maximum main iterations exceeded       
%  3: Maximum line-search iterations exceeded 
%
%  4: Cost reduction smaller than tolerance     
%  5: Gradient smaller than tolerance    
%  6: All controls are clamped 

m = size(Quu,1);

% Initialize arrays
clamped = false(m,1);      % Indicies of clamped controls
prev_clamped = false(m,1);
free = true(m,1);
delta_u = zeros(m, 1);
grad = zeros(m, 1);
grad_clamped = zeros(m, 1);
u_c = zeros(m, 1);
Luu = zeros(m, m);         % Placeholder to return if Luu not assigned

% Initialize scalars
old_cost = 0;
result = 0;

% Solver Options
max_iters = 100;          % max iterations
min_grad = 1e-8;          % min norm of non-clamped gradient
min_rel_improve = 1e-8;   % min relative improvement
step_dec = 0.6;           % factor for decreasing stepsize
min_step = 1e-20;         % min stepsize for linesearch
armijo_tol = 0.1;   	  % Armijo tolerance (fraction of linear improvement required)

% Initial controls
u = clamp(u0(:), lower_lim, upper_lim);

% Initial cost value
cost = u'*Qu + 0.5*u'*Quu*u;

% Start optimisation
for iter = 1:max_iters
    if result ~=0
        break;
    end
    
    % Check relative cost change for convergence
    if(iter > 1 && (old_cost - cost) < min_rel_improve*abs(old_cost))
        result = 4;
        break;
    end
    old_cost = cost;
    
    % Gradient of cost function
    grad = Qu + Quu*u;
    
    % Find clamped controls
    prev_clamped(:) = clamped;
    clamped(:, 1) = false;
    clamped((u == lower_lim) & (grad > 0)) = true;
    clamped((u == upper_lim) & (grad < 0)) = true;
    free(:) = ~clamped;
    
    % Check if all controls clamped
    if all(clamped)
        result = 6;
        break;
    end
    
    % Cholesky factorize if clamped controls have changed
    if iter == 1
        factorize = true;
    else
        factorize = any(prev_clamped ~= clamped);
    end
    
     % Cholesky (check for non-PD)
    if factorize
        [Luu, indef] = chol_free(Quu(free,free));
        if indef
            result = 1;
            break
        end
    end
    
    % check gradient-norm of free controls
    grad_norm = norm(grad(free));
    if grad_norm < min_grad
        result = 5;
        break;
    end
    
    % get search direction
    grad_clamped = Qu  + Quu*(u.*clamped);
    delta_u(:) = 0;
    delta_u(free) = -chol_solve(Luu, grad_clamped(free)) - u(free); % cholesky solver
    
    % check projected change in cost is a reduction
    expected_change = sum(delta_u.*grad);
    if expected_change >= 0 % (should not happen)
        result = 0;
        break
    end
    
    % Armijo linesearch
    step = 1;
	u_c = clamp(u + step*delta_u, lower_lim, upper_lim);
    cost_c = u_c'*Qu + 0.5*u_c'*Quu*u_c;
    while (cost_c - old_cost)/(step*expected_change) < armijo_tol
        step  = step*step_dec;
		u_c = clamp(u + step*delta_u, lower_lim, upper_lim);
        cost_c = u_c'*Qu + 0.5*u_c'*Quu*u_c;
        if step < min_step
            result = 3;
            break
        end
    end
    
    % Update the control
    u = u_c;
    cost = cost_c;
end

if iter >= max_iters
    result = 2;
end
end

function [clampedVals] = clamp(x, lower, upper)
% Returns array x with all values clamped between lower and upper

clampedVals = max(lower, min(upper, x));

end


function [L, fail] = chol_free(A)
% Wrapper for MATLAB chol for use with auto coder
% Inputs:
%===========
% A - positive semi-definite matrix
[L, fail] = chol(A, 'lower');

end


