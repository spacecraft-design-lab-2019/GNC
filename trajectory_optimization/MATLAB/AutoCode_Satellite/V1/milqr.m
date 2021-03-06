function [x,u,K,result] = milqr(x0, xg, u0, u_lims) %#codegen
% Solves finite horizon optimal control problem using a
% multiplicative iterative linear quadratic regulator

% Note that this solver will not work without control limits

% Inputs
% ===========================================
% x0 - The intial trajectory (n, N)
%
% xg - The goal state (n, 1)
%
% u0 - The initial control sequeunce (m, N-1)
%
% u_lims - The control limits (m, 2) (lower, upper)

% Outputs
% ===========================================
% x - Final nominal trajectory (n, N)
%
% u - Final open-loop controls (m, N-1)
%
% K - Feedback control gains (n-1, m, N-1) 
%
% result - Indicates convergence (boolean)


% Options (pass in as array)
dt = 0.03;            % Timestep (Should be highest we can get away with)
max_iters = 500;      % maximum iterations
exit_tol = 1e-7;      % cost reduction exit tolerance
grad_tol = 1e-4;      % gradient exit criterion
lambda_tol = 1e-5;    % lambda criterion for gradient exit
z_min = 0;            % minimum accepted cost reduction ratio
lambda_max = 1e10;    % maximum regularization parameter
lambda_min = 1e-6;    % set lambda = 0 below this value
lambda_scaling = 1.6; % amount to scale dlambda by

% CONSTANTS
Alphas = 10.^linspace(0, -3, 11);  % line search param
lambda = 1;
dlambda = 1;
N = size(u0, 2)+1;
Nx = size(x0, 1);
Nu = size(u0, 1);
Ne = Nx-1;  % error state size (3 param. error representation for attitude)

% Init matrices for update (otherwise MATLAB coder throws an error)
x_n = zeros(Nx,N);
u_n = zeros(Nu,N-1);
fx_n = zeros(Ne,Ne,N-1);
fu_n = zeros(Ne,Nu,N-1);
cx_n = zeros(Ne,N);
cu_n = zeros(Nu,N-1);
cxx_n = zeros(Ne,Ne,N);
cuu_n = zeros(Nu,Nu,N-1);
cost_n = 0;

% Initial Forward rollout
l = zeros(Nu, N-1);
K = zeros(Nu, Ne, N-1);
dV = zeros(1, 2);
alpha = 0;
[x,u,fx,fu,cx,cu,cxx,cuu,cost] = forwardRollout(x0,xg,u0,l,K,alpha,u_lims,dt);

% Convergence check params
expectedChange = 0; % Expected cost change
z = 0;              % Ratio of cost change to expected cost change
result = false;

for iter = 1:max_iters
    
    % Backward Pass
    %=======================================
    backPassDone = false;
    while ~backPassDone
        [l,K,dV,diverge] = backwardPass(fx,fu,cx,cu,cxx,cuu,lambda,u_lims,u);
        
        if diverge
            % Increase regularization parameter (lambda)
            dlambda = max(lambda_scaling * dlambda, lambda_scaling);
            lambda = max(lambda * dlambda, lambda_min);
            if lambda > lambda_max
                break;
            end
            continue;  % Retry with larger lambda
        end
        backPassDone = true;
    end
    
    % Check gradient of control, defined as l/u
    % Terminate if sufficiently small (success)
    g_norm = mean(max(abs(l)./(abs(u)+1),[],1)); % Avg of max grad at each time step
    if g_norm < grad_tol && lambda < lambda_tol
        result = true;
        break;
    end
   
    % Forward Line-Search
    %===========================================
    fwdPassDone = false;
    if backPassDone
        for alpha = Alphas
            [x_n,u_n,fx_n,fu_n,cx_n,cu_n,cxx_n,cuu_n,cost_n] = forwardRollout(x,xg,u,l,K,alpha,u_lims,dt);
            expectedChange = -alpha*(dV(1) + alpha*dV(2));
            if expectedChange > 0
                z = (cost - cost_n)/expectedChange;
            else
                z = sign(cost - cost_n);
            end
            if z > z_min
                fwdPassDone = true;
                break;
            end
        end
    end
    
    % Parameter Updates
    %=============================================
    if fwdPassDone
        % Decrease Lambda
        dlambda = min(dlambda/lambda_scaling, 1/lambda_scaling);
        lambda = lambda * dlambda * (lambda > lambda_min);  % set = 0 if lambda too small
        dcost = cost - cost_n;
        
        % Update trajectory and controls
        x = x_n;
        u = u_n;
        fx = fx_n;
        fu = fu_n;
        cx = cx_n;
        cu = cu_n;
        cxx = cxx_n;
        cuu = cuu_n;
        cost = cost_n;
        
        % Terminate if cost change sufficiently small
        if dcost < exit_tol
            result = true;
            break;
        end
    else
        % No cost reduction (based on z-value)
        % Increase lambda
        dlambda = max(lambda_scaling * dlambda, lambda_scaling);
        lambda = max(lambda * dlambda, lambda_min);
        
        if lambda > lambda_max
            % Lambda too large - solver diverged
            result = false;
            break;
        end
        
    end
end

if iter == max_iters
    % Ddin't converge completely
    result = false;
end
    
end

function [xnew,unew,fx,fu,cx,cu,cxx,cuu,cost] = forwardRollout(x,xg,u,l,K,alpha,u_lims,dt) %#codegen
% Uses an rk method to roll out a trajectory
% Returns the new trajectory, cost and the derivatives along the trajectory

% If feed-forward and feed-back controls l, K are non-zero
% then the unew returned is the new control sequeunce. Otherwise unew = u

% Sizes
N = size(x,2);
Nx = size(x,1);
Nu = size(u,1);
Ne = Nx-1;

% Initialize outputs (Don't overwrite x or u)
xnew = zeros(Nx,N);
unew = zeros(Nu,N-1);
fx = zeros(Ne,Ne,N-1);
fu = zeros(Ne,Nu,N-1);
cx = zeros(Ne,N);
cu = zeros(Nu,N-1);
cxx = zeros(Ne,Ne,N);
cuu = zeros(Nu,Nu,N-1);
cost = 0;

xnew(:,1) = x(:,1);
terminal = 0;
dx = zeros(6,1);
for k = 1:(N-1)
    
    % Find the state error vector dx
    dx(4:6) = xnew(5:7,k) - x(5:7,k);
    dx(1:3) = quat_error(xnew(1:4,k), x(1:4,k));
    
    % Find the new control and ensure it is within the limits
    unew(:,k) = u(:,k) - alpha*l(:,k) - K(:,:,k)*dx;
    unew(:,k) = min(u_lims(:,2), max(u_lims(:,1), unew(:,k)));

    % Step the dynamics forward
    [xnew(:,k+1),fx(:,:,k),fu(:,:,k)] = satellite_step(xnew(:,k), unew(:,k), dt);
    
    % Calculate the cost
    [c, cx(:,k),cu(:,k), cxx(:,:,k), cuu(:,:,k)] = satellite_cost(xnew(:,k), xg, unew(:,k), terminal); 
    cost = cost + c;
end

% Final cost
terminal = 1;
u_temp = zeros(Nu,1);
[c,cx(:,N),~,cxx(:,:,N),~] = satellite_cost(xnew(:,N), xg, u_temp, terminal); 
cost = cost + c;

function [dq] = quat_error(qnew, q_nom)
% Calculate error between qnew and q_nom
% Defined as conj(q_nom)*qnew

q_inv = [1;-1;-1;-1].*q_nom;  % conjugate
q_error = L_mult(q_inv)*qnew;
q_error = q_error/sqrt(q_error'*q_error);  % re-normalize
dq = q_error(2:4) / q_error(1);  % inverse Cayley Map
end

end


function [l, K, dV, diverge] = backwardPass(fx,fu,cx,cu,cxx,cuu,lambda,u_lims,u) %#codegen
% Perfoms the LQR backward pass to find the optimal controls
% Solves a quadratic program (QP) at each timestep for the optimal
% controls given the control limits

N = size(u, 2) + 1;
Ne = size(fx,1);
Nu = size(u,1);

% Initialize matrices (for C)
l = zeros(Nu,N-1);
K = zeros(Nu,Ne,N-1);
Qx = zeros(Ne,1);
Qu = zeros(Nu,1);
Qxx = zeros(Ne,Ne);
Quu = zeros(Nu,Nu);
Qux = zeros(Nu,Ne);

Kk = zeros(Nu,Ne);
result = 0;

% Change in cost
dV = [0 0];

% Set cost-to-go Jacobain and Hessian equal to final costs
Vx = cx(:, N);
Vxx = cxx(:,:,N);

diverge = false;
for k=(N-1):-1:1
    
    % Define cost gradients
    Qx = cx(:,k) + fx(:,:,k)'*Vx;
    Qu = cu(:,k) + fu(:,:,k)'*Vx;
    Qxx = cxx(:,:,k) + fx(:,:,k)'*Vxx*fx(:,:,k);
    Quu = cuu(:,:,k) + fu(:,:,k)'*Vxx*fu(:,:,k);
    Qux = fu(:,:,k)'*Vxx*fx(:,:,k);
    
    % Regularization (for Cholesky positive definiteness)
    QuuF = Quu + eye(Nu)*lambda;
    
    % Solve the Quadratic program with control limits
    upper = u_lims(:,2) - u(:,k);
    lower = u_lims(:,1) - u(:,k);
    l_idx = min(N-1, k+1);
    [lk,result,Luu,free] = boxQPsolve(QuuF,Qu,lower,upper,-1*l(:,l_idx));

    if result < 1
        diverge = true;
        return;
    end
    
    % Solve for feedback gains in non-clamped rows of u
    % (using cholesky factor of Quu)
    Kk(:,:) = 0;
    if any(free)
        Kk(free, :) = -chol_solve(Luu, Qux(free,:));
    end
    
    % Update Cost to Go
    dV  = dV + [lk'*Qu  (1/2)*lk'*Quu*lk];
    Vx  = Qx  + Kk'*Quu*lk + Kk'*Qu  + Qux'*lk;
    Vxx = Qxx + Kk'*Quu*Kk + Kk'*Qux + Qux'*Kk;
    Vxx = (1/2)*(Vxx + Vxx');  % Make sure Hessian is symmetric
    
    % Update Control Vectors
    l(:, k) = -lk;
    K(:,:,k) = -Kk;
    
end

end



