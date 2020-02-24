function [x, u, L, Vx, Vxx, cost, trace, stop] = iLQG(DYNCST, x0, u0, Op)
% iLQG - solve the deterministic finite-horizon optimal control problem.
%
%        minimize sum_i CST(x(:,i),u(:,i)) + CST(x(:,end))
%            u
%        s.t.  x(:,i+1) = DYN(x(:,i),u(:,i))
%
% Inputs
% ======
% DYNCST - A combined dynamics and cost function. It is called in
% three different formats.
%
%  1) step:
%   [xnew,c] = DYNCST(x,u,i) is called during the forward pass. 
%   Here the state x and control u are vectors: size(x)==[n 1],  
%   size(u)==[m 1]. The cost c and time index i are scalars.
%   If Op.parallel==true (the default) then DYNCST(x,u,i) is be 
%   assumed to accept vectorized inputs: size(x,2)==size(u,2)==K
%  
%  2) final:
%   [~,cnew] = DYNCST(x,nan) is called at the end the forward pass to compute
%   the final cost. The nans indicate that no controls are applied.
%  
%  3) derivatives:
%   [~,~,fx,fu,fxx,fxu,fuu,cx,cu,cxx,cxu,cuu] = DYNCST(x,u,I) computes the
%   derivatives along a trajectory. In this case size(x)==[n N+1] where N
%   is the trajectory length. size(u)==[m N+1] with NaNs in the last column
%   to indicate final-cost. The time indexes are I=(1:N).
%   Dimensions match the variable names e.g. size(fxu)==[n n m N+1]
%   note that the last temporal element N+1 is ignored for all tensors
%   except cx and cxx, the final-cost derivatives.
%
% x0 - The initial state from which to solve the control problem. 
% Should be a column vector. If a pre-rolled trajectory is available
% then size(x0)==[n N+1] can be provided and Op.cost set accordingly.
%
% u0 - The initial control sequence. A matrix of size(u0)==[m N]
% where m is the dimension of the control and N is the number of state
% transitions. 
%
%
% Op - optional parameters, see below
%
% Outputs
% =======
% x - the optimal state trajectory found by the algorithm.
%     size(x)==[n N+1]
%
% u - the optimal open-loop control sequence.
%     size(u)==[m N]
%
% L - the optimal closed loop control gains. These gains multiply the
%     deviation of a simulated trajectory from the nominal trajectory x.
%     size(L)==[m n N]
%
% Vx - the gradient of the cost-to-go. size(Vx)==[n N+1]
%
% Vxx - the Hessian of the cost-to-go. size(Vxx)==[n n N+1]
%
% cost - the costs along the trajectory. size(cost)==[1 N+1]
%        the cost-to-go is V = fliplr(cumsum(fliplr(cost)))
%
% lambda - the final value of the regularization parameter
%
% trace - a trace of various convergence-related values. One row for each
%         iteration, the columns of trace are
%         [iter lambda alpha g_norm dcost z sum(cost) dlambda]
%         see below for details.
%
% timing - timing information
%
%
%
% BIBTeX:
%
% @INPROCEEDINGS{
% author={Tassa, Y. and Mansard, N. and Todorov, E.},
% booktitle={Robotics and Automation (ICRA), 2014 IEEE International Conference on},
% title={Control-Limited Differential Dynamic Programming},
% year={2014}, month={May}, doi={10.1109/ICRA.2014.6907001}}

%---------------------- user-adjustable parameters ------------------------
defaults = {'lims',           [],...            control limits
            'parallel',       true,...          use parallel line-search?
            'Alpha',          10.^linspace(0,-3,11),... backtracking coefficients
            'tolFun',         1e-7,...          reduction exit criterion
            'tolGrad',        1e-4,...          gradient exit criterion
            'maxIter',        500,...           maximum iterations            
            'lambda',         1,...             initial value for lambda
            'dlambda',        1,...             initial value for dlambda
            'lambdaFactor',   1.6,...           lambda scaling factor
            'lambdaMax',      1e10,...          lambda maximum value
            'lambdaMin',      1e-6,...          below this value lambda = 0
            'regType',        1,...             regularization type 1: q_uu+lambda*eye(); 2: V_xx+lambda*eye()
            'zMin',           0,...             minimal accepted reduction ratio
            'diffFn',         [],...            user-defined diff for sub-space optimization
            'plot',           1,...             0: no;  k>0: every k iters; k<0: every k iters, with derivs window
            'print',          2,...             0: no;  1: final; 2: iter; 3: iter, detailed
            'plotFn',         @(x)0,...         user-defined graphics callback
            'cost',           [],...            initial cost for pre-rolled trajectory            
            };

% --- initial sizes and controls
n   = size(x0, 1);          % dimension of state vector
m   = size(u0, 1);          % dimension of control vector
N   = size(u0, 2);          % number of state transitions
u   = u0;                   % initial control sequence

% --- proccess options
if nargin < 4
    Op = struct();
end
Op  = setOpts(defaults,Op);

verbosity = Op.print;

switch numel(Op.lims)
    case 0
    case 2*m
        Op.lims = sort(Op.lims,2);
    case 2
        Op.lims = ones(m,1)*sort(Op.lims(:))';
    case m
        Op.lims = Op.lims(:)*[-1 1];
    otherwise
        error('limits are of the wrong size')
end

lambda   = Op.lambda;
dlambda  = Op.dlambda;

% --- initialize trace data structure
trace = struct('iter',nan,'lambda',nan,'dlambda',nan,'cost',nan,...
        'alpha',nan,'grad_norm',nan,'improvement',nan,'reduc_ratio',nan,...
        'time_derivs',nan,'time_forward',nan,'time_backward',nan);
trace = repmat(trace,[min(Op.maxIter,1e6) 1]);
trace(1).iter = 1;
trace(1).lambda = lambda;
trace(1).dlambda = dlambda;

% --- initial trajectory
diverge = true;
for alpha = Op.Alpha
    [x,un,cost]  = forward_pass(x0(:,1),alpha*u,[],[],[],1,DYNCST,Op.lims,[]);
    % simplistic divergence test
    if all(abs(x(:)) < 1e8)
        u = un;
        diverge = false;
        break
    end
end

trace(1).cost = sum(cost(:));

if diverge
    [Vx,Vxx, stop]  = deal(nan);
    L        = zeros(m,n,N);
    cost     = [];
    fprintf('\nEXIT: Initial control sequence caused divergence\n');
    return
end

% constants, timers, counters
flgChange   = 1;
dcost       = 0;
z           = 0;
expected    = 0;


fprintf('\n=========== begin iLQG ===========\n');

for iter = 1:Op.maxIter 
    
    %====== STEP 1: differentiate dynamics and cost along new trajectory
    if flgChange
        [~,~,fx,fu,fxx,fxu,fuu,cx,cu,cxx,cxu,cuu]   = DYNCST(x, [u nan(m,1)], 1:N+1);
        flgChange   = 0;
    end
    
    %====== STEP 2: backward pass, compute optimal control law and cost-to-go
    backPassDone   = 0;
    while ~backPassDone
        [diverge, Vx, Vxx, l, L, dV] = back_pass(cx,cu,cxx,cxu,cuu,fx,fu,fxx,fxu,fuu,lambda,Op.regType,Op.lims,u);
        
        if diverge
            fprintf('Cholesky failed at timestep %d.\n',diverge);
            dlambda   = max(dlambda * Op.lambdaFactor, Op.lambdaFactor);
            lambda    = max(lambda * dlambda, Op.lambdaMin);
            if lambda > Op.lambdaMax
                break;
            end
            continue
        end
        backPassDone      = 1;
    end

    % check for termination due to small gradient
    g_norm         = mean(max(abs(l) ./ (abs(u)+1),[],1));
    if g_norm < Op.tolGrad && lambda < 1e-5
        fprintf('\nSUCCESS: gradient norm < tolGrad\n');
        break;
    end
    
    %====== STEP 3: line-search to find new control sequence, trajectory, cost
    fwdPassDone  = 0;
    if backPassDone
    % serial backtracking line-search
        for alpha = Op.Alpha
            [xnew,unew,costnew]   = forward_pass(x0 ,u+l*alpha, L, x(:,1:N),[],1,DYNCST,Op.lims,Op.diffFn);
            dcost    = sum(cost(:)) - sum(costnew(:));
            expected = -alpha*(dV(1) + alpha*dV(2));
            if expected > 0
                z = dcost/expected;
            else
                z = sign(dcost);
                warning('non-positive expected reduction: should not occur');
            end
            if (z > Op.zMin)
                fwdPassDone = 1;
                break;
            end
        end
        if ~fwdPassDone
            alpha = nan; % signals failure of forward pass
        end
    end
    %====== STEP 4: accept step (or not), draw graphics, print status
 
    if fwdPassDone
        % decrease lambda
        dlambda   = min(dlambda / Op.lambdaFactor, 1/Op.lambdaFactor);
        lambda    = lambda * dlambda * (lambda > Op.lambdaMin);
        
        % accept changes
        u              = unew;
        x              = xnew;
        cost           = costnew;
        flgChange      = 1;
        
        % terminate ?
        if dcost < Op.tolFun
            fprintf('\nSUCCESS: cost change < tolFun\n');
            break;
        end
        
    else % no cost improvement
        % increase lambda
        dlambda  = max(dlambda * Op.lambdaFactor, Op.lambdaFactor);
        lambda   = max(lambda * dlambda, Op.lambdaMin); 
        
        % terminate ?
        if lambda > Op.lambdaMax
            fprintf('\nEXIT: lambda > lambdaMax\n');
            break;
        enddlda
    end
end

if iter == Op.maxIter
    fprintf('\nEXIT: Maximum iterations reached.\n');
end


if ~isempty(iter)
    printf("done")
else
    error('Failure: no iterations completed, something is wrong.')
end


function [xnew,unew,cnew] = forward_pass(x0,u,L,x,du,Alpha,DYNCST,lims,diff)
% parallel forward-pass (rollout)
% internally time is on the 3rd dimension, 
% to facillitate vectorized dynamics calls

n        = size(x0,1);
K        = length(Alpha);
K1       = ones(1,K); % useful for expansion
m        = size(u,1);
N        = size(u,2);

xnew        = zeros(n,K,N);
xnew(:,:,1) = x0(:,ones(1,K));
unew        = zeros(m,K,N);
cnew        = zeros(1,K,N+1);
for i = 1:N
    unew(:,:,i) = u(:,i*K1);
    
    if ~isempty(du)
        unew(:,:,i) = unew(:,:,i) + du(:,i)*Alpha;
    end    
    
    if ~isempty(L)
        dx          = xnew(:,:,i) - x(:,i*K1);
        unew(:,:,i) = unew(:,:,i) + L(:,:,i)*dx;
    end
    
    if ~isempty(lims)
        unew(:,:,i) = min(lims(:,2*K1), max(lims(:,1*K1), unew(:,:,i)));
    end

    [xnew(:,:,i+1), cnew(:,:,i)]  = DYNCST(xnew(:,:,i), unew(:,:,i), i*K1);
end
[~, cnew(:,:,N+1)] = DYNCST(xnew(:,:,N+1),nan(m,K,1),i);
% put the time dimension in the columns
xnew = permute(xnew, [1 3 2]);
unew = permute(unew, [1 3 2]);
cnew = permute(cnew, [1 3 2]);



function [diverge, Vx, Vxx, k, K, dV] = back_pass(cx,cu,cxx,cxu,cuu,fx,fu,fxx,fxu,fuu,lambda,regType,lims,u)
% Perform the Ricatti-Mayne backward pass

N  = size(cx,2);
n  = numel(cx)/N;
m  = numel(cu)/N;

cx    = reshape(cx,  [n N]);
cu    = reshape(cu,  [m N]);
cxx   = reshape(cxx, [n n N]);
cxu   = reshape(cxu, [n m N]);
cuu   = reshape(cuu, [m m N]);

k     = zeros(m,N-1);
K     = zeros(m,n,N-1);
Vx    = zeros(n,N);
Vxx   = zeros(n,n,N);
dV    = [0 0];

Vx(:,N)     = cx(:,N);
Vxx(:,:,N)  = cxx(:,:,N);

diverge  = 0;
for i = N-1:-1:1
    
    Qu  = cu(:,i)      + fu(:,:,i)'*Vx(:,i+1);
    Qx  = cx(:,i)      + fx(:,:,i)'*Vx(:,i+1);
    Qux = cxu(:,:,i)'  + fu(:,:,i)'*Vxx(:,:,i+1)*fx(:,:,i);
    
    Quu = cuu(:,:,i)   + fu(:,:,i)'*Vxx(:,:,i+1)*fu(:,:,i);
    Qxx = cxx(:,:,i)   + fx(:,:,i)'*Vxx(:,:,i+1)*fx(:,:,i);
    
    Vxx_reg = (Vxx(:,:,i+1) + lambda*eye(n)*(regType == 2));  % Regularization type 1
    Qux_reg = cxu(:,:,i)'   + fu(:,:,i)'*Vxx_reg*fx(:,:,i);

    QuuF = cuu(:,:,i)  + fu(:,:,i)'*Vxx_reg*fu(:,:,i) + lambda*eye(m)*(regType == 1);  % (OR!) Regularization type 2
    

    % solve Quadratic Program (for control with limits)
    lower = lims(:,1)-u(:,i);
    upper = lims(:,2)-u(:,i);

    [k_i,result,R,free] = boxQP(QuuF,Qu,lower,upper,k(:,min(i+1,N-1)));
    if result < 1
        diverge  = i;  % Quu non positive definite
        return;
    end

    % Feedback gains in free rows
    K_i    = zeros(m,n);
    if any(free)
        Lfree        = -R\(R'\Qux_reg(free,:));
        K_i(free,:)   = Lfree;
    end
    
    % update cost-to-go approximation
    dV          = dV + [k_i'*Qu  .5*k_i'*Quu*k_i];
    Vx(:,i)     = Qx  + K_i'*Quu*k_i + K_i'*Qu  + Qux'*k_i;
    Vxx(:,:,i)  = Qxx + K_i'*Quu*K_i + K_i'*Qux + Qux'*K_i;
    Vxx(:,:,i)  = .5*(Vxx(:,:,i) + Vxx(:,:,i)');
    
    % save controls/gains
    k(:,i)      = k_i;
    K(:,:,i)    = K_i;
end


