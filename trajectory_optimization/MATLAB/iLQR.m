function [x, u, K, Jhist, trace] = iLQR(DYNAMICS, COST, x0, xg, u0, uLims, Ops)
% Solves finite horizon optimal control problem using the iterative
% linear quadratic redualtor method

% Note that this solver will not work without control limits

% Inputs
% ===========================================
% DYNAMICS - Function handle for the rkstep/ dynamics update
%           (In final implementation remove this and just put dynamcis
%           and cost functions at bottom of this file)
%
% COST - Function handle for cost calculation/ cost derivatives
%
% x0 - The intial trajectory (n, N)
%
% xg - The goal state (n, 1)
%
% u0 - The initial control sequeunce (m, N-1)
%
% uLims - The control limits (m, 2) (lower, upper)
%
% Ops (options):
% -----------------------
% dt, max_iters, exit_tol, grad_tol, zMin, lambda_max, lambda_min, lambda_scaling


% Outputs
% ===========================================
% x - Final nominal trajectory (n, N)
%
% u - Final open-loop controls (m, N-1)
%
% K - Feedback control gains (n, m, N-1)
%
% Jhist - The cost history (convergence)
%
% trace - Error messages



% CONSTANTS
Alphas = 10.^linspace(0, -3, 11);  % line search param
lambda = 1;
dlambda = 1;
N = size(u0, 2) + 1;
Nx = size(x0, 1);
Nu = size(u0, 1);


% Initial Forward rollout
% Returns xtraj, (unew=utraj0), cost
l = zeros(Nu, N-1);
K = zeros(Nu, Nx, N-1);
alpha = 0;
[x,u,fx,fu,cx,cu,cxx,cuu,cost] = forwardRollout(DYNAMICS,COST,x0,xg,u0,l,K,alpha,uLims,Ops.dt);
Jhist(1) = cost;


% Convergence checks and constants
flgChange   = 1;  % True if cost improves (otherwise change lambda)
dcost       = 0;  % Cost change
expectedChange  = 0;  % Expected cost change
z           = 0;  % Ratio of cost change to expected cost change


printf("\n==================Begin iLQR================\n");
flgChange = 1;
for iter = 1:Ops.max_iters
    
    % Backward Pass
    backwardPassDone = 0;
    while ~backwardPassDone
        [l, K, dV, diverge] = backwardPass(fx,fu,cx,cu,cxx,cuu,lambda,uLims,u)
    end


end

if iter == Ops.max_iters
    



end

end

function [xnew,unew,fx,fu,cx,cu,cxx,cuu,cost] = forwardRollout(DYNAMICS,COST,x,xg,u,l,K,alpha,uLims,dt)
% Uses an rk method to roll out a trajectory
% Returns the new trajectory, cost and the derivatives along the trajectory

% If feed-forward and feed-back controls l, K are non-zero
% then the unew returned is the new control sequeunce. Otherwise unew = u

% Sizes
N = size(x,2);
Nx = size(x,1);
Nu = size(u,1);

% Initialize outputs
xnew = zeros(Nx,N);
unew = zeros(Nu,N-1);
fx = zeros(Nx,Nx,N-1);
fu = zeros(Nx,Nu,N-1);
cx = zeros(Nx,1,N);
cu = zeros(Nu,1,N-1);
cxx = zeros(Nx,Nx,N);
cuu = zeros(Nu,Nu,N-1);
cost = 0;

xnew(:,1) = x(:,1);
terminal = 0;
for k = 1:(N-1)
    
    % Update the control during line-search
    % (During inital forward rollout, l and K will be arrays of zeros)
    unew(:,k) = u(:,k) - alpha*l(:,k) - K(:,:,k)*(xnew(:,k) - x(:,k));

    % Ensure control is within limits
    unew(:,k) = min(uLims(:,2), max(uLims(:,1), unew(:,k)));

    % Step the dynamics forward
    [xnew(:,k+1),fx(:,:,k),fu(:,:,k)] = DYNAMICS(xnew(:,k), unew(:,k), dt);
    
    % Calculate the cost
    [c, cx(:,k),cu(:,k), cxx(:,:,k), cuu(:,:,k)] = COST(xnew(:,k), unew(:,k), xg, terminal); 
    cost = cost + c;
end

% Final cost
terminal = 1;
u_temp = zeros(Nu,1);
[c,cx(:,N),~,cxx(:,:,N),~] = COST(xnew(:,N), u_temp, xg, terminal); 
cost = cost + c;

end


function [l, K, dV, diverge] = backwardPass(fx,fu,cx,cu,cxx,cuu,lambda,uLims,u)
% Perfoms the LQR backward pass to find the optimal controls

N = size(u, 2) + 1;
Nx = size(fx,1);
Nu = size(u,1);

% Initialize matrices
l = zeros(Nu,N-1);
K = zeros(Nu,Nx,N-1);
Qx = zeros(Nx,1);
Qu = zeros(Nu,1);
Qxx = zeros(Nx,Nx);
Quu = zeros(Nu,Nu);
Qux = zeros(Nu,Nx);

% Set cost-to-go Jacobain and Hessian equal to final costs
Vx = cx(:, N);
Vxx = cxx(:,:,N);

for k=(N-1):-1:1
    
    % Define cost gradients
    Qx = cx(:,k) + fx(:,:,k)'*Vx;
    Qu = cu(:,k) + fu(:,:,k)'*Vx;
    Qxx = cxx(:,:,k) + fx(:,:,k)'*Vxx*fx(:,:,k);
    Quu = cuu(:,:,k) + fu(:,:,k)'*Vxx*fu(:,:,k);
    Qux = fu(:,:,k)'*Vxx*fx(:,:,k);
    
    % Regularization (for Cholesky positive definiteness)
    Quu = Quu + eye(Nu)*lambda;
    
    % Solve the Quadratic program with control lims
    upper = uLims(:,2) - u(:,k);
    lower = uLims(:,1) - u(:,k);
    [l,result,Luu,freeIdcs] = boxQPsolve(Quu,Qu,lower,upper,u0);
    
    if result < 1
        diverge = k;
        return;
    end
    
    % Solve for feedback gains in non-clamped rows of u
    % (using cholesky factor of Quu)
    if any(freeIdcs)
        K(freeIdcs,:,k) = Luu\Luu'\Qux(freeIdcs,:);
    end
    
    
    
    
end



end



