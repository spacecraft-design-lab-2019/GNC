clear
close all
clc

% Paramters
T = 501;                    % time horizon
x0 = zeros(4, T);
x0(:,1) = [1;1;pi*3/2;0];   % initial state
xg = zeros(4,1);            % goal state
u0 = zeros(2,T-1);            % initial controls
u_lims = [-.5 .5;           % wheel angle limits (radians)
          -2  2];           % acceleration limits (m/s^2)

% Run iLQR
[x,u,K,result] = ilqrCar(x0, xg, u0, u_lims);

result

figure(1)
plot(x(1,:), x(2,:))
