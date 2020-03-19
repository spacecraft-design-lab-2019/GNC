function [x, fx, fu] = satellite_step(x0,u0,dt)
    % Steps the dynamics forward using a 2nd order rk-method
    % Returns: new state, discrete time Jacobians
    
    Nx = length(x0);
    Nu = length(u0);
    
    % Step dynamics
    %Explicit midpoint step from x_k to x_{k+1}
    [xdot1, dxdot1] = satellite_dynamics(x0,u0);
    [xdot2, dxdot2] = satellite_dynamics(x0+.5*dt*xdot1,u0);
    x1 = x0 + dt*xdot2;
    
    % Normalize the quaternion
    q1 = x1(1:4);
    x = [q1/sqrt(q1'*q1); x1(5:7)];
    
    % Continuous time Jacobians
    A1 = dxdot1(:,(1:Nx));
    B1 = dxdot1(:,Nx+(1:Nu));
    A2 = dxdot2(:,(1:Nx));
    B2 = dxdot2(:,Nx+(1:Nu));
    
    % Discrete time Jacobians
    % First form the state error Jacobians E(x_k) & E(x_k+1)
    E0 = [-x0(2:4)', zeros(1,3);
          x0(1)*eye(3) + skew_mat(x0(2:4)), zeros(3,3);
          zeros(3,3), eye(3)];
    E1 = [-x1(2:4)', zeros(1,3);
          x1(1)*eye(3) + skew_mat(x1(2:4)), zeros(3,3);
          zeros(3,3), eye(3)];
     
    fx = E1'*(eye(Nx) + dt*A2 + 0.5*dt*dt*A2*A1)*E0;
    fu = E1'*(dt*B2 + 0.5*dt*dt*A2*B1);
    
end


function [xdot, dxdot] = satellite_dynamics(x,u)
% Calculates the continuous time state derivative and Jacobians

J = 0.01*eye(3); % kgm^2
Jinv = inv(J);

% Angular velocity
w = x(5:7);

% Quaternion components
s = x(1);
v = x(2:4);

% Non-linear dynamics
qdot = 0.5*[-v'; s*eye(3) + skew_mat(v)]*w;
wdot = Jinv*(u - skew_mat(w)*J*w);
xdot = [qdot; wdot];

% Jacobians 
A = 0.5*[0, -w', -v';
         w, -skew_mat(w), s*eye(3)+skew_mat(v);
         zeros(3,4), -2*Jinv*(skew_mat(w)*J - skew_mat(J*w))];

B = [zeros(4,3);
     Jinv];   %!!!!! This changes slightly with magnetorquers !!!!

dxdot = [A, B];

end
