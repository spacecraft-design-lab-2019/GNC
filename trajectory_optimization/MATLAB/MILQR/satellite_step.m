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

function [x, fx, fu] = satellite_step(x0,u0,dt,B_ECI)
% Steps the dynamics forward using a 2nd order rk-method
% Returns: new state, discrete time Jacobians

Nx = length(x0);
Nu = length(u0);

% Step dynamics
% Explicit midpoint step from x_k to x_{k+1}
[xdot1, dxdot1] = satellite_dynamics(x0,u0,B_ECI);
[xdot2, dxdot2] = satellite_dynamics(x0+.5*dt*xdot1,u0,B_ECI);
x1 = x0 + dt*xdot2;

% Re-normalize the quaternion
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


function [xdot, dxdot] = satellite_dynamics(x,u,B_ECI)
% Calculates the continuous time state derivative and Jacobians

% Inputs
%=======================
% x - the current state
% u - the current control (magentic moment vector)
% B - Earth magnetic field vector in ECI co-ordinates (3x1)

% TODO: ideally, we'd like to pass in a B_B (body frame magnetic field
% vector)

J = 0.01*eye(3); % kgm^2
Jinv = inv(J);

w = x(5:7);  % Angular velocity
q = x(1:4);  % Quaternion

% Magnetic field section 
% given B_ECI (ECI magnetic field at the time step)
B_B = qrot(x(1:4),B_ECI);

% Non-linear dynamics
qdot = 0.5*[-q(2:4)'; q(1)*eye(3) + skew_mat(q(2:4))]*w;

% wdot = Jinv*(u - skew_mat(w)*J*w);
wdot = Jinv*(-cross(B_B,u) - skew_mat(w)*J*w); % with magnetorquers
xdot = [qdot; wdot];

% Jacobians 
A = 0.5*[0, -w', -q(2:4)';
         w, -skew_mat(w), q(1)*eye(3)+skew_mat(q(2:4));
         zeros(3,4), -2*Jinv*(skew_mat(w)*J - skew_mat(J*w))];

% B = [zeros(4,3);
%      Jinv];
B = [zeros(4,3);  % with magnetorquers:
     -Jinv*skew_mat(B_B)];
 
dxdot = [A, B];

end

function output = qrot(q,vec)
% this function multiplie s a vector by a quaternion rotation

q_cross = [q(1);-q(2:4)];

output = qmult(q_cross,qmult([0;vec],q));
output = output(2:4);

end

function q_out = qmult(q1,q2)
% this function multiplies quaternions

q_out = [q1(1)*q2(1) - dot(q1(2:4),q2(2:4));...
    q1(1)*q2(2:4) + q2(1)*q1(2:4) + cross(q1(2:4),q2(2:4))];

end
