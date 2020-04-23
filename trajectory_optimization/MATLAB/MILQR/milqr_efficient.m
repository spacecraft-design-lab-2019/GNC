 

function [x,u,K,result] = milqr_efficient(x0, xg, u0, u_lims,dt,B_ECI)
% Solves finite horizon optimal control problem using a
% multiplicative iterative linear quadratic regulator

% Note that this solver will not work without control limits

% Inputs
% ===========================================
% x0 - The intial state (n, 1)
%
% xg - The goal state (n, 1)
%
% u0 - The initial control sequeunce (m, N-1)
%
% u_lims - The control limits (m, 2) (lower, upper)
%
% B_ECI - A time sequence of Earth magnetic field vectors in ECI (3, N)


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
% dt = 0.03;            % Timestep (Should be highest we can get away with)
max_iters = 300;      % maximum iterations
cost_tol = 1e-7;      % cost reduction exit tolerance
contr_tol = 1e-4;     % feedforward control change exit criterion
lambda_tol = 1e-5;    % max regularizaion param allowed for exit
c_ratio_min = 0;      % minimum accepted cost reduction ratio
lambda_max = 1e7;    % maximum regularization parameter
lambda_min = 1e-6;    % set lambda = 0 below this value
lambda_scale = 1.6;   % amount to scale dlambda by

% Init optimisation params
alphas = 10.^linspace(0,-3,11);  % line search param vector
lambda = 1;
d_lambda = 1;
N = size(u0, 2)+1;
Nx = size(x0, 1);
Nu = size(u0, 1);
Ne = Nx-1;  % error state size (3 param. error representation for attitude)

% variable initialization
cost = 0;
cost_n = 0;

% initialize x and u
x = zeros(Nx,N);
u = zeros(Nu,N-1);

% initialize the "next" variables
x_n = x0(:,ones(N,1));
u_n = u0;
fx_n = zeros(Ne,Ne,N-1);
fu_n = zeros(Ne,Nu,N-1);
cx_n = zeros(Ne,N);
cu_n = zeros(Nu,N-1);
cxx_n = zeros(Ne,Ne,N);
cuu_n = zeros(Nu,Nu,N-1);

% Initialize Forward rollout
l = zeros(Nu, N-1);
K = zeros(Nu, Ne, N-1);
dV = zeros(1, 2);
alpha = 0;

% Throw out everything after x and u, then in backward pass evalute on the
% spot
[x,u,cost] = forwardRollout(x0(:,ones(N,1)),xg,u0,l,K,alpha,u_lims,dt,B_ECI);

% Convergence check params
expected_change = 0;      % Expected cost change
c_ratio = 0;              % Ratio of cost change to expected cost change
result = false;

fprintf("\n=====Running MILQR Optimisation====\n");
for iter = 1:max_iters
    fprintf("\n---New Iteration---");
    fprintf("\n lambda: ");
    fprintf(string(lambda));
    if exist('dcost')
        fprintf("\n cost change: ");
        fprintf(string(dcost));
    end
    fprintf("\n");

    
    
    % Backward Pass
    %=======================================
    backPassCheck = false;
    while ~backPassCheck
        % don't pass in gradient info, calculate in place
        [l,K,dV,diverged] = backwardPass(lambda,u_lims,u,x,xg,dt,B_ECI);
        if diverged
            % fprintf("---Warning: Cholesky factorizaton failed---\n");
            
            % Increase regularization parameter (lambda)
            lambda = updateLambda(lambda,1,d_lambda,lambda_scale,lambda_min);
            if lambda > lambda_max
                break;
            end
            continue;  % Retry with larger lambda
        end
        backPassCheck = true;
    end
    
    % Check relative change of feedforward control
    % Terminate if sufficiently small (success)
    c_norm = mean(max(abs(l)./(abs(u)+1),[],1)); % Avg over time of max change
    if  lambda < lambda_tol && c_norm < contr_tol
        % fprintf("\n---Success: Control change decreased below tolerance---\n");
        result = true;
        break;
    end
   
    % Forward Line-Search
    %===========================================
    lineSearchCheck = false;
    if backPassCheck
        for alpha = alphas
            [x_n,u_n,cost_n] = ...
                forwardRollout(x,xg,u,l,K,alpha,u_lims,dt,B_ECI);
            expected_change = alpha*dV(1) + (alpha^2)*dV(2);
            if expected_change < 0
                c_ratio = (cost_n - cost)/expected_change;
            else
                % Non-positive expected cost reduction
                % actual cost change must be negative to accept the step
                c_ratio = -sign(cost_n - cost);
            end
            if c_ratio > c_ratio_min
                lineSearchCheck = true;
                break;
            end
        end
    end
    
    % Parameter Updates
    %=============================================
    if lineSearchCheck
        % Decrease Lambda
        lambda = updateLambda(lambda,-1,d_lambda,lambda_scale,lambda_min);
        dcost = cost - cost_n;
        
        % Update the trajectory and controls
        x = x_n;
        u = u_n;
        fx = fx_n;
        fu = fu_n;
        cx = cx_n;
        cu = cu_n;
        cxx = cxx_n;
        cuu = cuu_n;
        cost = cost_n;
        
        % Change in cost small enough to terminate?
        if dcost < cost_tol
            result = true;
            % fprintf('\n---Success: cost change < tolerance---\n');
            break;
        end
        
    else
        % No cost reduction (based on cost change ratio)
        % Increase lambda
        lambda = updateLambda(lambda,1,d_lambda,lambda_scale,lambda_min);
        if lambda > lambda_max
            result = false;
            % fprintf("\n---Diverged: new lambda > lambda_max---\n");
            break;
        end
        
    end
end

if iter == max_iters
    % Ddin't converge completely
    result = false;
%     % fprintf("\n---Warning: Max iterations exceeded---\n");
end

function [lambda] = updateLambda(lambda, direction, delta, l_scale, l_min)
    % Increases or decreases the regularization parameter according
    % to a non-linear scaling regime.
    
    if (direction == 1)  % increase lambda
        delta = max(l_scale * delta, l_scale);
        lambda = max(lambda * delta, l_min);
        
    elseif (direction == -1) % decrease lambda
        delta = min(delta/l_scale, 1/l_scale);
        lambda = lambda * delta * (lambda > l_min);  % set = 0 if lambda too small
    end
end
    
end

% get rid of everything after unew except cost
function [xnew,unew,cost] = ...
    forwardRollout(x,xg,u,l,K,alpha,u_lims,dt,B_ECI)
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
cost = 0;

xnew(:,1) = x(:,1);
terminal = 0;
dx = zeros(6,1);
for k = 1:(N-1)
    % Find the state error vector dx
    dx(4:6) = xnew(5:7,k) - x(5:7,k);
    dx(1:3) = quat_error(xnew(1:4,k),x(1:4,k));
    
    % Find the new control and ensure it is within the limits
    unew(:,k) = u(:,k) - alpha*l(:,k) - K(:,:,k)*dx;
    unew(:,k) = min(u_lims(:,2), max(u_lims(:,1),unew(:,k)));

    % Step the dynamics forward
    % recreate dynamics function to separate into dynamics step and 
    % derivative function to get derivatives
    [xnew(:,k+1)] = satellite_step_efficient(xnew(:,k), unew(:,k), dt, B_ECI(:,k));
    
    % Calculate the cost
    % make second function that's satellite cost derivatives
    [c] = satellite_cost_efficient(xnew(:,k),xg,unew(:,k),terminal); 
    cost = cost + c;
end

% Final cost
terminal = 1;
u_temp = zeros(Nu,1);
[c] = satellite_cost_efficient(xnew(:,N), xg, u_temp, terminal); 
cost = cost + c;

function [dq] = quat_error(qk, q_nom)
% Calculate error between qk and q_nom
% Defined as conj(q_nom)*qnew
% Returns error as Rodrigues parameters (3x1)

q_inv = [1;-1;-1;-1].*q_nom;               % conjugate
q_error = L_mult(q_inv)*qk;
q_error = q_error/sqrt(q_error'*q_error);  % re-normalize
dq = q_error(2:4) / q_error(1);            % inverse Cayley Map
end
end


function [l,K,dV,diverged] = backwardPass(lambda,...
    u_lims,u,x,xg,dt,B_ECI)
% Perfoms the LQR backward pass to find the optimal controls
% Solves a quadratic program (QP) at each timestep for the optimal
% controls given the control limits

N = size(u, 2) + 1;
Ne = size(x,1)-1;
Nu = size(u,1);

% insert section where I evaluate the cost and dynamics derivatives
% function [cx,cu,cxx,cuu,cxu] = cost_derivatives(x,u)
% function [fx,fu] = state_derivatives(x,u)

% Initialize matrices (for C code, not needed in MATLAB)
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

% Set cost-to-go Jacobian and Hessian equal to final costs
terminal = 1;
[cx, ~, cxx, ~] = satellite_cost_derivatives(x(:,N), xg, zeros(Nu,1), terminal);
Vx = cx;
Vxx = cxx;

diverged = false;
for k=(N-1):-1:1
    
    % calculate dynamics derivatives
    [fx, fu] = satellite_derivatives(x(:,k),u(:,k),dt,B_ECI(:,k));
    
    % Define cost gradients and hessians
    % convert cost and dynamics derivatives into functions
    [cx, cu, cxx, cuu] = satellite_cost_derivatives(x(:,k), xg, u(:,k), terminal);
    Qx = cx + fx'*Vx;
    Qu = cu + fu'*Vx;
    Qxx = cxx + fx'*Vxx*fx;
    Quu = cuu + fu'*Vxx*fu;
    Qux = fu'*Vxx*fx;
    
    % Regularization (for Cholesky positive definiteness)
    QuuR = Quu + eye(Nu)*lambda;
    
    % Solve the Quadratic program with control limits
    upper_lim = u_lims(:,2) - u(:,k);
    lower_lim = u_lims(:,1) - u(:,k);
    l_idx = min(N-1, k+1);
    [lk,result,Luu,free] = boxQPsolve(QuuR,Qu,lower_lim,upper_lim,-1*l(:,l_idx));

    if result < 2
        diverged = true;
        % fprintf('\nDiverged with lambda = %f\n',lambda);
        return;
    end
    
    % Solve for feedback gains in non-clamped rows of u
    % (using cholesky factor of Quu)
    Kk(:,:) = 0;
    if any(free)
        Kk(free, :) = -chol_solve(Luu, Qux(free,:));
    end
    
    % Update Cost to Go Jacobian and Hessian
    Vx  = Qx  + Kk'*Quu*lk + Kk'*Qu  + Qux'*lk;
    Vxx = Qxx + Kk'*Quu*Kk + Kk'*Qux + Qux'*Kk;
    Vxx = (1/2)*(Vxx + Vxx');  % Ensure Hessian is symmetric
    
    % Record control cost change to check convergence
    dV  = dV + [lk'*Qu  (1/2)*lk'*Quu*lk];  
    
    % Update Control Vectors
    l(:, k) = -lk;
    K(:,:,k) = -Kk;
    
end

end



