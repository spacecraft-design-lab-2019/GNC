% Copyright (c) 2020 Robotic Exploration Lab
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.

% Testing satllite attitude control

clear
close all
clc

addpath('utils')
addpath('sgp4') % not needed but just in case ;)

% Sim Params
N = 1E4;  % num steps
Nx = 7;
Nu = 3;
dt = 1;
t = 0:dt:N*dt;

% Initialize Parameters
Earth = InitializeEarth();

% Initial State trajectory
theta = pi/2;  % [rad] 
axis = [0;0;1];
r0 = [5.7198;-0.7556;4.1212]*1E3; % km
v0 = [-4.4390;-1.0905;5.945]; % km/s
q0 = [cos(theta/2); axis*sin(theta/2)];  % 90 degree rotation about z-axis
w0 = [0; 0; 0];  % [rad/s]
x0 = [q0; w0];
MJD0 = 58947.77014; % modified julian date on April 8, 2020
t_MJD = t/60/60/24+MJD0;

% Goal state
theta_g = pi;
rg = [0;0;1];
qg = [cos(theta_g/2); rg*sin(theta_g/2)];
wg = [0; 0; 0];
xg = [qg; wg];

% magnetic field (ECI)

% IGRF magnetic field
[B_eci] = get_magnetic_field_series([r0;v0;x0],t_MJD);

% B_eci = [40E-6*sin(.04*(1:N)); 60E-6*sin(.04*(1:N)); 60E-6*cos(.04*(1:N))];

plot(B_eci');

% Initial Control trajectory
u0 = zeros(Nu, N-1);
u_lims = [-100 100;           % magnetic moment limits
          -100 100;
          -100 100];
      
% Run MILQR
x0 = x0(:,ones(N,1)); % extend x0 to entire trajectory to fill
[x,u,K,result] = milqr(x0, xg, u0, u_lims, B_eci);

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


