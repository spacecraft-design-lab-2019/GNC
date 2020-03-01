% Testing a bicyle model car with 1st order dynamics

clear
close all
clc

% Paramters
T = 501;                    % time horizon
x0 = zeros(4, T);
x0(:,1) = [1;1;pi*3/2;0];   % initial state
% xg = zeros(4,1);            % goal state
xg = [0;0;0;0];
u0 = zeros(2,T);            % initial controls
u_lims = [-.5 .5;           % wheel angle limits (radians)
          -2  2];           % acceleration limits (m/s^2)

% Options
Ops.dt = 0.03;            % Timestep (Should match MCU Hz)
Ops.max_iters = 500;      % maximum iterations
Ops.exit_tol = 1e-7;      % cost reduction exit tolerance
Ops.grad_tol = 1e-4;      % gradient exit criterion
Ops.lambda_tol = 1e-5;    % lambda criterion for gradient exit
Ops.z_min = 0;            % minimum accepted cost reduction ratio
Ops.lambda_max = 1e10;    % maximum regularization parameter
Ops.lambda_min = 1e-6;    % set lambda = 0 below this value
Ops.lambda_scaling = 1.6; % amount to scale dlambda by


% Run iLQR
[x,u,K,Jhist,result] = iLQR(@car_step, @car_cost2, x0, xg, u0, u_lims, Ops);

% Plot results
% Car path
figure(1)
plot(x(1,:),x(2,:));
grid on
hold on
s = 0.5;
quiver(x0(1,1),x0(2,1), cos(x0(3,1)), sin(x0(3,1)), 'g');
quiver(xg(1),xg(2), s*cos(xg(3)), s*sin(xg(3)), 'r')
quiver(x(1,end), x(2,end), cos(x(3,end)), sin(x(3,end)), 'b')


% State and control trajectories
figure(2)
subplot(2,1,1)
hold on
plot(x(1, :));
plot(x(2, :));
plot(x(3, :));
plot(x(4, :));
grid on
legend("x", "y", "theta", "v");
% Controls
subplot(2,1,2)
plot(u(1,:));
hold on
plot(u(2,:));
legend("Steer angle", "Accel");

% Cost history
figure(3)
plot(Jhist);
grid on


