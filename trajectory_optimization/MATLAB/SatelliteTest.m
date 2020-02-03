clear;
clc;
close all;

% Sim steps
N = 250;

% Initial State
theta = pi/2;  % [rad]
q0 = [cos(theta/2); 0; 0; sin(theta/2)];  % 90 degree rotation about z-axis
w0 = [0; 0; 0];  % [rad/s]
x0 = [q0; w0];

% Goal state
theta_g = 3*pi/2;
xg = [cos(theta_g/2); 0; 0; sin(theta_g/2)];

utraj0 = zeros(3, N-1);

% Cost matrices
Qw = 0.05*eye(3);
R = 0.3*eye(3);
Qf = 10*eye(3);

[xtraj, utraj, K, jhist] = iLQRsatellite(@satellite_dynamics, x0, xg, utraj0, Q, R, Qf, .01, 1e-3);






