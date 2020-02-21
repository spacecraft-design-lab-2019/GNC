% Nick Goodson
% 5th Feb 2020
%
% Tests the quaternion iLQR function for controlling the
% attitude of a satellite

clear;
clc;
close all;

% Sim Params
N = 5000;  % num steps
dt = 0.01;
tol = 1e-3;
max_iters = 1000;
Nx = 7;
Nu = 3;

% Initial State
theta = pi/2;  % [rad]
q0 = [cos(theta/2); 0; 0; sin(theta/2)]  % 90 degree rotation about z-axis
w0 = [0; 0; 0];  % [rad/s]
x0 = [q0; w0];

% Goal state
theta_g = 3*pi/2;
qg = [cos(theta_g/2); 0; 0; sin(theta_g/2)]
wg = [0; 0; 0];
xg = [qg; wg];

% Initial Control trajectory
utraj0 = zeros(Nu, N-1);

% Cost matrices
% Cumulative
Q = zeros(Nx, Nx);
Qw = 0.1*eye(3);
Q(5:7, 5:7) = Qw;
R = 0.5*eye(3);

% Terminal
Qf= zeros(Nx, Nx);
Qwf = 1*eye(3);
Qf(5:7, 5:7) = Qwf;
Qqf = 200;  % Final cost for attitude error

[xtraj, utraj, K, Jhist, success] = iLQRsatellite(x0, xg, utraj0, Q, R, Qf, Qqf, dt, tol, max_iters);

% print final state
xtraj(:, end)

% Plot the quaternion time evolution
figure(1);
sgtitle("Quaternion")
subplot(4,1,1) 
plot(xtraj(1,:));
subplot(4,1,2)
plot(xtraj(2, :));
subplot(4,1,3)
plot(xtraj(3,:));
subplot(4,1,4)
plot(xtraj(4,:));

%plot the controls
figure(2)
subplot(3,1,1)
sgtitle("Controls")
plot(utraj(1,:));
subplot(3,1,2)
plot(utraj(2,:));
subplot(3,1,3)
plot(utraj(3,:));

%plot omega
figure(3)
sgtitle("Omega")
subplot(3,1,1)
plot(xtraj(5, :))
subplot(3,1,2)
plot(xtraj(6, :))
subplot(3,1,3)
plot(xtraj(7,:))

%plot the cost(
figure(4)
sgtitle("Cost")
plot(Jhist)






