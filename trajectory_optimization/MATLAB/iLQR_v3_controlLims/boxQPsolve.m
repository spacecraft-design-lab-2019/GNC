function [u,result,Luu_free,free] = boxQPsolve(Quu,Qu,lower,upper,u0)
% Minimize 0.5*u'*Quu*u + u'*Qu  s.t. lower<=u<=upper
%
%  inputs:
%     Quu          - positive definite matrix   (n * n)
%     Qu           - bias vector                (n)
%     lower        - lower bounds               (n)
%     upper        - upper bounds               (n)
%     u0           - initial state              (n)
%
%  outputs:
%     u            - solution                   (n)
%     result       - result type (roughly, higher is better, see below)
%     Luu          - cholesky factor            (n * n)
%     free         - set of free dimensions     (n)

n = size(Quu,1);

% Initialize arrays
clamped      = false(n,1);
prev_clamped = false(n,1);
free         = true(n,1);
deltaX       = zeros(n, 1);
grad         = zeros(n, 1);
grad_clamped = zeros(n, 1);
uc           = zeros(n, 1);
Luu_free     = zeros(n, n);  % Placeholder to return if Luu_free not assigned

% Initialize scalars
oldvalue     = 0;
result       = 0;

% options
maxIter        = 100;       % maximum number of iterations
minGrad        = 1e-8;      % minimum norm of non-fixed gradient
minRelImprove  = 1e-8;      % minimum relative improvement
stepDec        = 0.6;       % factor for decreasing stepsize
minStep        = 1e-22;     % minimal stepsize for linesearch
Armijo         = 0.1;   	% Armijo parameter (fraction of linear improvement required)

% initial state
u = clamp(u0(:), lower, upper);

% initial objective value
value = u'*Qu + 0.5*u'*Quu*u;

% main loop
for iter = 1:maxIter
    
    if result ~=0
        break;
    end
    
    % check relative improvement
    if( iter>1 && (oldvalue - value) < minRelImprove*abs(oldvalue) )
        result = 3;
        break;
    end
    oldvalue = value;
    
    % get gradient
    grad = Qu + Quu*u;
    
    % find clamped dimensions
    prev_clamped                    = clamped;
    clamped(:, 1)                   = false;
    clamped((u == lower)&(grad>0))  = true;
    clamped((u == upper)&(grad<0))  = true;
    free                            = ~clamped;
    
    % check for all clamped
    if all(clamped)
        result = 5;
        break;
    end
    
    % Cholesky factorize if clamped controls have changed
    if iter == 1
        factorize = true;
    else
        factorize = any(prev_clamped ~= clamped);
    end
    
     % Cholesky (check for non PD)
    if factorize
        [Luu_free, indef] = chol_local(Quu(free, free));
        if indef
            result = -1;
            break
        end
    end
    
    % check gradient norm
    gnorm  = norm(grad(free));
    if gnorm < minGrad
        result = 4;
        break;
    end
    
    % get search direction
    grad_clamped = Qu  + Quu*(u.*clamped);
    deltaX(:, 1) = 0;
    deltaX(free) = -Luu_free\(Luu_free'\grad_clamped(free)) - u(free); % cholesky solver
    
    % check for descent direction
    sdotg = sum(deltaX.*grad);
    if sdotg >= 0 % (should not happen)
        result = 0;
        break
    end
    
    % Armijo linesearch
    step  = 1;
    nstep = 0;
	uc    = clamp(u + step*deltaX, lower, upper);
    vc    = uc'*Qu + 0.5*uc'*Quu*uc;
    while (vc - oldvalue)/(step*sdotg) < Armijo
        step  = step*stepDec;
        nstep = nstep+1;
		uc    = clamp(u + step*deltaX, lower, upper);
        vc    = uc'*Qu + 0.5*uc'*Quu*uc;
        if step < minStep
            result = 2;
            break
        end
    end
    
    % accept candidate
    u = uc;
    value = vc;
end

if iter >= maxIter
    result = 1;
end

% Results
% ===========================
% -1: Hessian is not positive definite
%  0: No descent direction found          (SHOULD NOT OCCUR)
%  1: Maximum main iterations exceeded       
%  2: Maximum line-search iterations exceeded  
%  3: Improvement smaller than tolerance     
%  4: Gradient norm smaller than tolerance    
%  5: All dimensions are clamped 

end

function [clampedVals] = clamp(x, lower, upper)
% Returns array x with all values clamped between lower and upper

clampedVals = max(lower, min(upper, x));

end


function [L, fail] = chol_local(A)
% Wrapper for MATLAB chol for use with auto coder
[L, fail] = chol(A);

end
