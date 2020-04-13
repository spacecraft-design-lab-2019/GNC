% Testing satllite attitude control with MILQR

clear
close all
clc
addpath('utils')

% Satellite Model Params
Nx = 7;
Nu = 3;
J = 0.01*eye(3); % Inertia kgm^2
Js = [J, inv(J)];

% Initial State trajectory
theta = pi/2;  % [rad] 
r0 = [0;0;1];
q0 = [cos(theta/2); r0*sin(theta/2)];  % 90 degree rotation about z-axis
w0 = [0; 0; 0];  % [rad/s]
x0 = [q0; w0];
x0 = x0(:,ones(N,1));

% Goal state
theta_g = pi;
rg = [0;0;1];
qg = [cos(theta_g/2); rg*sin(theta_g/2)];
wg = [0; 0; 0];
xg = [qg; wg];

% Initial Control trajectory
u0 = zeros(Nu, N-1);
u_lims = [-100 100;           % magnetic moment limits
          -100 100;
          -100 100];

% Solver options
N = 5000;             % num steps
dt = 0.03;            % Timestep (Should be highest we can get away with)
max_iters = 300;      % maximum iterations
cost_tol = 1e-7;      % cost reduction exit tolerance
contr_tol = 1e-4;     % feedforward control change exit criterion
lambda_tol = 1e-5;    % max regularizaion param allowed for exit
c_ratio_min = 0;      % minimum accepted cost reduction ratio
lambda_max = 1e10;    % maximum regularization parameter
lambda_min = 1e-6;    % set lambda = 0 below this value
lambda_scale = 1.6;   % amount to scale dlambda by

options = [N, dt, max_iters, cost_tol, contr_tol, lambda_tol, ...
           c_ratio_min, lambda_max, lambda_min, lambda_scale];

% magnetic field (ECI)
% to test, let's just get some sinusoids out here
B_ECI = [40E-6*sin(.04*(1:N)); 60E-6*sin(.04*(1:N)); 60E-6*cos(.04*(1:N))];

% Run MILQR
[x,u,K,result] = milqr(x0, xg, u0, u_lims, B_ECI, Js, options);

% Plot results
%==========================================================================
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


