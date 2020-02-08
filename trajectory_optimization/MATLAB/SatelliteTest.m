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
qg = [cos(theta_g/2); 0; 0; sin(theta_g/2)];
wg = [0; 0; 0];
xg = [qg; wg];

utraj0 = zeros(3, N-1);

% Cost matrices
Qw = 0.05*eye(3);
R = 0.3*eye(3);
Qwf = 10*eye(3);
Qqf = 30;

[xtraj, utraj, K, jhist] = iLQRsatellite(x0, xg, utraj0, Qw, R, Qwf, Qqf, .01, 1e-3);






