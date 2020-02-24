function [u,result,Quufree,free] = boxQPsolve(Quu,Qu,lower,upper,u0)
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
%     Quufree      - subspace cholesky factor   (n_free * n_free)
%     free         - set of free dimensions     (n)


n        = size(Quu,1);
clamped  = false(n,1);
free     = true(n,1);
oldvalue = 0;
result   = 0;
gnorm    = 0;
nfactor  = 0;
Quufree  = zeros(n);


% initial state
u = clamp(u0(:), lower, upper);

% options
maxIter        = 100;       % maximum number of iterations
minGrad        = 1e-8;      % minimum norm of non-fixed gradient
minRelImprove  = 1e-8;      % minimum relative improvement
stepDec        = 0.6;       % factor for decreasing stepsize
minStep        = 1e-22;     % minimal stepsize for linesearch
Armijo         = 0.1;   	% Armijo parameter (fraction of linear improvement required)

% initial objective value
value    = u'*Qu + 0.5*u'*Quu*u;

% main loop
for iter = 1:maxIter
    
    if result ~=0
        break;
    end
    
    % check relative improvement
    if( iter>1 && (oldvalue - value) < minRelImprove*abs(oldvalue) )
        result = 4;
        break;
    end
    oldvalue = value;
    
    % get gradient
    grad = Qu + Quu*u;
    
    % find clamped dimensions
    old_clamped                     = clamped;
    clamped                         = zeros(n,1);
    clamped((u == lower)&(grad>0))  = 1;
    clamped((u == upper)&(grad<0))  = 1;
    free                            = ~clamped;
    
    % check for all clamped
    if all(clamped)
        result = 6;
        break;
    end
    
    % factorize if clamped has changed
    if iter == 1
        factorize = true;
    else
        factorize = any(old_clamped ~= clamped);
    end
    
     % Cholesky (check for non PD)
    if factorize
        [Quufree, indef] = chol_local(Quu(free,free));
        if indef
            result = -1;
            break
        end
        nfactor = nfactor + 1;
    end
    
    % check gradient norm
    gnorm  = norm(grad(free));
    if gnorm < minGrad
        result = 5;
        break;
    end
    
    % get search direction
    grad_clamped = Qu  + Quu*(u.*clamped);
    deltaX = zeros(n,1);
    deltaX(free) = -Quufree\(Quufree'\grad_clamped(free)) - u(free); % cholesky solver
    
    % check for descent direction
    sdotg = sum(deltaX.*grad);
    if sdotg >= 0 % (should not happen)
        break
    end
    
    % armijo linesearch
    step  = 1;
    nstep = 0;
	uc    = clamp(u + step*deltaX, lower, upper);
    vc    = uc'*Qu + 0.5*uc'*Quu*uc;
    while (vc - oldvalue)/(step*sdotg) < Armijo
        step  = step*stepDec;
        nstep = nstep+1;
		uc    = clamp(u + step*deltaX, lower, upper);
        vc    = uc'*Qu + 0.5*uc'*Quu*uc;
        if step<minStep
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


results = { 'Hessian is not positive definite',...          % result = -1
            'No descent direction found',...                % result = 0    SHOULD NOT OCCUR
            'Maximum main iterations exceeded',...          % result = 1
            'Maximum line-search iterations exceeded',...   % result = 2
            'No bounds, returning Newton point',...         % result = 3
            'Improvement smaller than tolerance',...        % result = 4
            'Gradient norm smaller than tolerance',...      % result = 5
            'All dimensions are clamped'};                  % result = 6

fprintf('RESULT: %s.\niterations %d  gradient %-12.6g final value %-12.6g  factorizations %d\n',...
    results{result+2}, iter, gnorm, value, nfactor);
end

function [clampedVals] = clamp(x, lower, upper)
% Returns clamped values

clampedVals = max(lower, min(upper, x));

end


function [L, fail] = chol_local(A)
% Wrapper for MATLAB chol for use with auto coder
[L, fail] = chol(A);

end
