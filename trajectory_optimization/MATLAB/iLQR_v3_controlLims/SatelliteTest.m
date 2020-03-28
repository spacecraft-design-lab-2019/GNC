% Testing satllite attitude control

clear
close all
clc

% Parameters
% Sim Params
N = 500;  % num steps
Nx = 7;
Nu = 3;

% Initial State
theta = pi/2;  % [rad]
% q0 = [cos(theta/2); 0; 0; sin(theta/2)];  % 90 degree rotation about z-axis
r = [1;-2;3]/norm([1,-2, 3]);
q0 = [cos(theta/2); r*sin(theta/2)];
w0 = [0; 0; 0];  % [rad/s]
x0 = zeros(Nx, N);
x0(:,1) = [q0; w0];

% Goal state
theta_g = pi;
qg = [cos(theta_g/2); 0; 0; sin(theta_g/2)];
wg = [0; 0; 0];
xg = [qg; wg];

% Initial Control trajectory
u0 = zeros(Nu, N-1);
u_lims = [-.5 .5;           % torque limits
          -.5 .5;
          -.5 .5];

% Options
Ops.dt = 0.03;            % Timestep
Ops.max_iters = 500;      % maximum iterations
Ops.exit_tol = 1e-7;      % cost reduction exit tolerance
Ops.grad_tol = 1e-4;      % gradient exit criterion
Ops.lambda_tol = 1e-5;    % lambda criterion for gradient exit
Ops.z_min = 0;            % minimum accepted cost reduction ratio
Ops.lambda_max = 1e10;    % maximum regularization parameter
Ops.lambda_min = 1e-6;    % set lambda = 0 below this value
Ops.lambda_scaling = 1.6; % amount to scale dlambda by


% Run iLQR
[x,u,K,Jhist,result] = iLQR(@satellite_step, @satellite_cost, x0, xg, u0, u_lims, Ops);

% Plot results
figure(1)

% Quaternion
subplot(2,1,1)
plot(x(1,:))
hold on
plot(x(2,:))
plot(x(3,:))
plot(x(4,:))
legend('q1','q2','q3','q4');

% Control torques
subplot(2,1,2)
plot(u(1,:))
hold on
plot(u(2,:))
plot(u(3,:))
legend('u1','u2','u3');


