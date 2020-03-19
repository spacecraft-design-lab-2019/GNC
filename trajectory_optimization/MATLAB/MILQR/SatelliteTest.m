% Testing satllite attitude control

clear
close all
clc

% Parameters
% Sim Params
N = 500;  % num steps
Nx = 7;
Nu = 3;

% Initial State trajectory
theta = pi/2;  % [rad]
% q0 = [cos(theta/2); 0; 0; sin(theta/2)];  % 90 degree rotation about z-axis
r = [1;-2;3]/norm([1,-2, 3]);
q0 = [cos(theta/2); r*sin(theta/2)];
w0 = [0; 0; 0];  % [rad/s]
x0 = [q0; w0];
x0 = x0(:,ones(N,1));

% Goal state
theta_g = pi;
qg = [cos(theta_g/2); 0; 0; sin(theta_g/2)];
wg = [0; 0; 0];
xg = [qg; wg];

% Initial Control trajectory
u0 = zeros(Nu, N-1);
u_lims = [-.1 .1;           % torque limits
          -.1 .1;
          -.1 .1];

% Run iLQR
[x,u,K,result] = milqr(x0, xg, u0, u_lims);

% Plot results
figure(1)

% Quaternion
subplot(3,1,1)
plot(x(1,:))
hold on
plot(x(2,:))
plot(x(3,:))
plot(x(4,:))
legend('q1','q2','q3','q4');

% Angular velocities
subplot(3,1,2)
plot(x(5,:))
hold on
plot(x(6,:))
plot(x(7,:))
legend('omega1','omega2','omega3');

% Control torques
subplot(3,1,3)
plot(u(1,:))
hold on
plot(u(2,:))
plot(u(3,:))
legend('u1','u2','u3');


