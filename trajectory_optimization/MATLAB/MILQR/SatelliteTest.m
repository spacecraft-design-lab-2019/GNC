% Testing satllite attitude control

clear
close all
clc

% Sim Params
N = 5000;  % num steps
Nx = 7;
Nu = 3;

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
u_lims = [-10 10;           % magnetic moment limits
          -10 10;
          -10 10];
      
% magnetic field (ECI)
% to test, let's just get some sinusoids out here
% B_ECI = [40E-6*sin(.004*(1:N)); 60E-6*sin(.004*(1:N)); 60E-6*cos(.004*(1:N))];
B_ECI = rand(3, N)*0.1;

% Run MILQR
[x,u,K,result] = milqr(x0, xg, u0, u_lims,B_ECI);

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


