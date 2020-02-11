function [xtraj, utraj, K, Jhist] = iLQRsatellite(x0, xg, utraj0, Q, R, Qf, Qqf, dt, tol, max_iters)
%iLQR Trajectory Optimization for Satellite

Nx = length(x0);
Nu = size(utraj0,1);
N = size(utraj0,2)+1;

xtraj = zeros(Nx, N);
xtraj(:,1) = x0;

A = zeros(Nx, Nx, N-1);
B = zeros(Nx, Nu, N-1);

utraj = utraj0;

%First forward simulate with utraj0 to get initial matrices
quat_cost_sign = zeros(Nx);
J = 0;
for k = 1:(N-1)
    % Find the quaternion-error cost (record the sign for use in backward pass)
    [quat_cost, sign] = calc_quat_cost(xg, xtraj(:,k));
    quat_cost_sign(k) = sign;
    
    % Cumulative Cost
    J = J + quat_cost + (1/2)*(xtraj(:,k)-xg)'*Q*(xtraj(:,k)-xg) + (1/2)*utraj(:,k)'*R*utraj(:,k);
    [xtraj(:,k+1), A(:,:,k), B(:,:,k)] = rkstep(xtraj(:,k),utraj(:,k), dt);
end
% Terminal cost
[final_quat_cost, sign] = calc_quat_cost(xtraj(:,N), xg);
quat_cost_sign(N) = sign;
J = J + Qqf*final_quat_cost + (1/2)*(xtraj(:,N)-xg)'*Qf*(xtraj(:,N)-xg);
Jhist(1) = J;


% Set up backwards pass matrices
S = zeros(Nx, Nx);
s = zeros(Nx,1);
K = zeros(Nu,Nx,N-1);
l = (tol+1)*ones(Nu,N-1);
Snew = zeros(Nx,Nx);  % temp matrices
snew = zeros(Nx,1);
unew = zeros(Nu,N-1); % Line search temp matrices
xnew = zeros(Nx,N); 

% Backward Pass
iter = 0;
while max(abs(l)) > tol

    iter = iter + 1;
    % Add failure condition if max_iters exceeded
    
    %Set up backwards LQR pass
    S = Qf;
    s = Qf*(xtraj(:,N)-xg);  %!!!!! Fix this for quaternion !!!!!
    for k = (N-1):-1:1
        
        %Calculate cost gradients for this time step
        % !!!!!!!!!!Need to update this for quaternion !!!!!!!!!!
        q = Qw*(xtraj(:,k)-xg);
        r = R*utraj(:,k);
        
        %Calculate l and K
        l(:,k) = (R + B(:,:,k)'*S*B(:,:,k))\(r + B(:,:,k)'*s);
        K(:,:,k) = (R + B(:,:,k)'*S*B(:,:,k))\(B(:,:,k)'*S*A(:,:,k));
        
        %Calculate new S and s        
        Snew = Q + K(:,:,k)'*R*K(:,:,k) + (A(:,:,k)-B(:,:,k)*K(:,:,k))'*S*(A(:,:,k)-B(:,:,k)*K(:,:,k));
        snew = q - K(:,:,k)'*r + K(:,:,k)'*R*l(:,k) + (A(:,:,k)-B(:,:,k)*K(:,:,k))'*(s- S*B(:,:,k)*l(:,k));
        S = Snew;
        s = snew;

    end

    %Now do forward pass line search with new l and K
    xnew(:,1) = xtraj(:,1);
    alpha = 1;
    Jnew = J+1;
    while Jnew > J
        Jnew = 0;
        quat_cost_sign = zeros(Nx);
        for k = 1:N-1
            unew(:,k) = utraj(:,k) - alpha*l(:,k) - K(:,:,k)*(xnew(:,k)-xtraj(:,k));
            [xnew(:,k+1), A(:,:,k), B(:,:,k)] = rkstep(xnew(:,k),unew(:,k), dt);
            
            % Find the quaternion-error cost (record the sign for use in backward pass)
            [quat_cost, sign] = calc_quat_cost(xg, xnew(:,k));
            quat_cost_sign(k) = sign;
            
            % Cumulative Cost
            Jnew = Jnew + quat_cost + (1/2)*(xnew(:,k)-xg)'*Q*(xnew(:,k)-xg) + (1/2)*utraj(:,k)'*R*utraj(:,k);
        end
        % Terminal cost
        [quat_cost_final, sign] = calc_quat_cost(xnew(:,N), xg);
        quat_cost_sign(N) = sign;
        Jnew = Jnew + Qqf*quat_cost_final + (1/2)*(xnew(:,N)-xg)'*Qf*(xnew(:,N)-xg);
        alpha = (1/2)*alpha;    
    end
    xtraj = xnew;
    utraj = unew;
    J = Jnew;
    Jhist(iter+1) = J;
    
    disp([max(abs(l)) 2*alpha])
end

end

function [x, A, B] = rkstep(x0,u0, dt)
    %Explicit midpoint step from x_k to x_{k+1}
    
    Nx = length(x0);
    Nu = length(u0);

    [xdot1, dxdot1] = satellite_dynamics(0,x0,u0);
    [xdot2, dxdot2] = satellite_dynamics(0,x0+.5*dt*xdot1,u0);
    x1 = x0 + dt*xdot2;
    
    % Normalize the quaternion
    q1 = x1(1:4);
    x = [q1/sqrt(q1'*q1); x1(5:7)];
    
    dqnorm = eye(4)/sqrt(q1'*q1) - (q1*q1')/((q1'*q1)^(3/2));
    Nk = [dqnorm, zeros(4, 3);
        zeros(3, 4), eye(3)];
    
    A1 = dxdot1(:,(1:Nx));
    B1 = dxdot1(:,Nx+(1:Nu));
    
    A2 = dxdot2(:,(1:Nx));
    B2 = dxdot2(:,Nx+(1:Nu));
    
    A = Nk*(eye(Nx) + dt*A2 + 0.5*dt*dt*A2*A1);
    B = Nk*(dt*B2 + 0.5*dt*dt*A2*B1);
end


function [quat_cost, sign] = calc_quat_cost(x, xg)
    % Find the linear quaternion-error cost (qcost = min(1+qg'q, 1-qg'q))
    % Also record the sign for use in the backward pass
    
    quat_cost = xg(1:4)'*x(1:4);
    if (1.0 + quat_cost) < (1.0 - quat_cost)
        quat_cost = (1.0 + quat_cost);
        sign = 1;
    else
        quat_cost = (1.0 - quat_cost);
        sign = -1;
    end
    
end








