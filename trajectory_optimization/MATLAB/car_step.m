function [x,fx,fu] = car_step(x0,u0,dt)
% RK step for car
% returns x, fx, fu

% === states and controls:
% x = [x y th v]' = [x; y; car_angle; front_wheel_velocity]
% u = [delta a]'     = [front_wheel_angle; acceleration]

% constants
l = 2.0;       % wheelbase

% states
x = x0(1);
y = x0(2);
th = x0(3);
v = x0(4);

% controls
delta  = u0(1); % front wheel angle
a  = u0(2);     % front wheel acceleration

xdot = zeros(4, 1);
xdot(1) = v*cos(delta)*cos(th);  % xdot
xdot(2) = v*cos(delta)*sin(th);  % ydot
xdot(3) = v*sin(delta)/l;        % thetadot
xdot(4) = a;                     % vdot

% Euler step
x = x0 + xdot * dt;

% Jacobians
A = zeros(4,4);
A(1,3) = -v*cos(delta)*sin(th);
A(1,4) = cos(delta)*cos(th);
A(2,3) = v*cos(delta)*cos(th);
A(2,4) = cos(delta)*sin(th);
A(3,4) = sin(delta)/l;

B = zeros(4,2);
B(1,1) = -v*sin(delta)*cos(th);
B(2,1) = -v*sin(delta)*sin(th);
B(3,1) = v*cos(delta)/l;
B(4,2) = 1;

fx = eye(4) + dt*A;
fu = dt*B;

end
