function [xtraj, utraj, K, Jhist] = iLQRsatellite(x0, xg, utraj0, Qw, R, Qwf, Qqf, dt, tol)
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
    
    % Find the linear quaternion-error cost
    % Also record the sign for use in the backward pass
    quat_cost = xg(1:4)'*xtraj(1:4,k);
    if (1.0 + quat_cost) < (1.0 - quat_cost)
        quat_cost = (1.0 + quat_cost);
        quat_cost_sign(k) = 1;
    else
        quat_cost = (1.0 - quat_cost);
        quat_cost_sign(k) = -1;
    end
    
    J = J + quat_cost + (1/2)*(xtraj(5:7,k)-xg(5:7))'*Qw*(xtraj(5:7,k)-xg(5:7)) + (1/2)*utraj(:,k)'*R*utraj(:,k);
    [xtraj(:,k+1), A(:,:,k), B(:,:,k)] = rkstep(xtraj(:,k),utraj(:,k), dt);
end
qN = xtraj(1:4,N);
wN = xtraj(5:7,N);
J = J + Qqf*min((1+qg'*qN), (1-qg'*qN)) + (1/2)*(wN-wg)'*Qwf*(wN-wg);
Jhist(1) = J;


% Set up backwards pass matrices
S = zeros(Nx, Nx);
s = zeros(Nx);
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
    
    %Set up backwards LQR pass
    S = Qf;
    s = Qf*(xtraj(:,N)-xg);  %!!!!! Fix this
    for k = (N-1):-1:1
        
        %Calculate cost gradients for this time step
        % !!!!!!!!!!Need to update this for quaternion
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
        for k = 1:N-1
            unew(:,k) = utraj(:,k) - alpha*l(:,k) - K(:,:,k)*(xnew(:,k)-xtraj(:,k));
            [xnew(:,k+1), A(:,:,k), B(:,:,k)] = rkstep(xnew(:,k),unew(:,k), dt);
            
            % Split state vector into quaternion and angular velocity
            qk = xtraj(1:4,k);
            wk = xtraj(5:7,k);
    
            % Find the linear quaternion-error cost
            % Also record the sign for use in the backward pass
            quat_cost = qg'*qk;
            if (1.0 + quat_cost) < (1.0 - quat_cost)
                quat_cost = (1.0 + quat_cost);
                quat_cost_sign(k) = 1;
            else
                quat_cost = (1.0 - quat_cost);
                quat_cost_sign(k) = -1;
            end
            
            Jnew = Jnew + quat_cost + (1/2)*(wk-wg)'*Qw*(wk-wg) + (1/2)*unew(:,k)'*R*unew(:,k);
        end
        qN = xnew(1:4,N);
        wN = xnew(5:7,N);
        Jnew = Jnew + Qqf*min((1+qg'*qN), (1-qg'*qN)) + (1/2)*(wN-wg)'*Qwf*(wN-wg);
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
