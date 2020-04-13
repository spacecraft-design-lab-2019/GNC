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

function [x,u,K,result] = milqr(x0, xg, u0, u_lims,B_ECI)
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
dt = 0.03;            % Timestep (Should be highest we can get away with)
max_iters = 300;      % maximum iterations
cost_tol = 1e-7;      % cost reduction exit tolerance
contr_tol = 1e-4;     % feedforward control change exit criterion
lambda_tol = 1e-5;    % max regularizaion param allowed for exit
c_ratio_min = 0;      % minimum accepted cost reduction ratio
lambda_max = 1e10;    % maximum regularization parameter
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

% Initial Forward rollout
l = zeros(Nu, N-1);
K = zeros(Nu, Ne, N-1);
dV = zeros(1, 2);
alpha = 0;
[x,u,fx,fu,cx,cu,cxx,cuu,cost] = forwardRollout(x0,xg,u0,l,K,alpha,u_lims,dt,B_ECI);

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
        [l,K,dV,diverged] = backwardPass(fx,fu,cx,cu,cxx,cuu,lambda,u_lims,u);
        if diverged
            fprintf("---Warning: Cholesky factorizaton failed---\n");
            
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
        fprintf("\n---Success: Control change decreased below tolerance---\n");
        result = true;
        break;
    end
   
    % Forward Line-Search
    %===========================================
    lineSearchCheck = false;
    if backPassCheck
        for alpha = alphas
            [x_n,u_n,fx_n,fu_n,cx_n,cu_n,cxx_n,cuu_n,cost_n] = ...
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
            fprintf('\n---Success: cost change < tolerance---\n');
            break;
        end
        
    else
        % No cost reduction (based on cost change ratio)
        % Increase lambda
        lambda = updateLambda(lambda,1,d_lambda,lambda_scale,lambda_min);
        if lambda > lambda_max
            result = false;
            fprintf("\n---Diverged: new lambda > lambda_max---\n");
            break;
        end
        
    end
end

if iter == max_iters
    % Ddin't converge completely
    result = false;
    fprintf("\n---Warning: Max iterations exceeded---\n");
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


function [xnew,unew,fx,fu,cx,cu,cxx,cuu,cost] = ...
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
    dx(1:3) = quat_error(xnew(1:4,k),x(1:4,k));
    
    % Find the new control and ensure it is within the limits
    unew(:,k) = u(:,k) - alpha*l(:,k) - K(:,:,k)*dx;
    unew(:,k) = min(u_lims(:,2), max(u_lims(:,1),unew(:,k)));

    % Step the dynamics forward
    [xnew(:,k+1),fx(:,:,k),fu(:,:,k)] = satellite_step(xnew(:,k), unew(:,k), dt, B_ECI(:,k));
    
    % Calculate the cost
    [c, cx(:,k),cu(:,k),cxx(:,:,k),cuu(:,:,k)] = satellite_cost(xnew(:,k),xg,unew(:,k),terminal); 
    cost = cost + c;
end

% Final cost
terminal = 1;
u_temp = zeros(Nu,1);
[c,cx(:,N),~,cxx(:,:,N),~] = satellite_cost(xnew(:,N), xg, u_temp, terminal); 
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


function [l,K,dV,diverged] = backwardPass(fx,fu,cx,cu,cxx,cuu,lambda,u_lims,u)
% Perfoms the LQR backward pass to find the optimal controls
% Solves a quadratic program (QP) at each timestep for the optimal
% controls given the control limits

N = size(u, 2) + 1;
Ne = size(fx,1);
Nu = size(u,1);

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
Vx = cx(:, N);
Vxx = cxx(:,:,N);

diverged = false;
for k=(N-1):-1:1
    
    % Define cost gradients and hessians
    Qx = cx(:,k) + fx(:,:,k)'*Vx;
    Qu = cu(:,k) + fu(:,:,k)'*Vx;
    Qxx = cxx(:,:,k) + fx(:,:,k)'*Vxx*fx(:,:,k);
    Quu = cuu(:,:,k) + fu(:,:,k)'*Vxx*fu(:,:,k);
    Qux = fu(:,:,k)'*Vxx*fx(:,:,k);
    
    % Regularization (for Cholesky positive definiteness)
    QuuR = Quu + eye(Nu)*lambda;
    
    % Solve the Quadratic program with control limits
    upper_lim = u_lims(:,2) - u(:,k);
    lower_lim = u_lims(:,1) - u(:,k);
    l_idx = min(N-1, k+1);
    [lk,result,Luu,free] = boxQPsolve(QuuR,Qu,lower_lim,upper_lim,-1*l(:,l_idx));

    if result < 2
        diverged = true;
        fprintf('\nDiverged with lambda = %f\n',lambda);
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



