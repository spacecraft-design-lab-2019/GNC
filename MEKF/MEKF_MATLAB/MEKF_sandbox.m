% Prototype script for Matlab MEKF
% To be autocoded into C
% Paul DeTrempe
clear; clc;

%% Test rotation matrix
% q0 = .5;
% q1 = .5;
% q2 = .5;
% q3 = .5;
% 
% % Test Jayden's stuff vs. regular quaternion rotation
% q = [q0; q1; q2; q3];
% R1 = quat2DCM(q);
% R2 = quat2DCM2(q);
% % ^^ These don't agree....

%% Load data
load('noise.mat')
load('simstates.mat')
load('predictions.mat')
load('sensors.mat')
load('bias.mat')
dt = .1;

% Tune Q and R matrices
% Q = Q;

% gyro noise
% R(4:6, 4:6) = 1e3*R(4:6, 4:6);


qtrue = simstates(:,4:7)';
N = size(qtrue, 2);
btrue = bias';
whist = simstates(:,11:13)';

%% Load Zach's data
% clear;
% load mekf_truth
% load mekf_inputs
% Q = W;
% R = V;
% N = size(qtrue, 2);
times = 0:dt:(N-1)*dt;
% 
% % Make quaternions scalar 1st
% q1 = qtrue(:,1)
% qtrue = circshift(qtrue,1,1);
% q1_shifted = qtrue(:,1)

%% Run the MEKF
% Initial conditions
x_k = [qtrue(:,1); .04*ones(3,1)];
P_k = blkdiag((.1*pi/180)^2*eye(3), (1*pi/180)^2*eye(3)); %10 deg. and 10 deg/sec 1-sigma uncertainty
w_k = whist(:,1);

% Preallocate
b_hist = zeros(size(btrue));
q_hist = zeros(size(qtrue));
P_hist = zeros(size(P_k,1), N);

for i = 1:N
    % measure some stuff
    r_B_body = sensors(i,1:3)';
    r_sun_body = sensors(i,7:9)';
    w_k = sensors(i,4:6)';
    r_sun_inert = predictions(i,4:6)';
    r_B_inert = predictions(i,1:3)';
    
    % Update our beliefs
    [x_k1, P_k1] = MEKFstep(x_k, P_k, w_k, r_sun_body, r_B_body,...
                                 r_sun_inert, r_B_inert, Q, R, dt);
                             
    % record history
    q_hist(:,i) = x_k1(1:4);
    b_hist(:,i) = x_k1(5:7);
    P_hist(:,i) = diag(P_k1);
    
    % Iterate
    x_k = x_k1;
    P_k = P_k1;
end

%% Plots
figure;
title('Quaternions')
subplot(4,1,1)
plot(times, qtrue(1,:))
hold on
plot(times, q_hist(1,:));
legend('q_1, true', 'q_1, estimate')
subplot(4,1,2)
plot(times, qtrue(2,:))
hold on
plot(times, q_hist(2,:));
legend('q_2, true', 'q_2, estimate')
subplot(4,1,3)
plot(times, qtrue(3,:))
hold on
plot(times, q_hist(3,:));
legend('q_3, true', 'q_3, estimate')
subplot(4,1,4)
plot(times, qtrue(4,:))
hold on
plot(times, q_hist(4,:));
legend('q_4, true', 'q_4, estimate')

figure;
plot(times, btrue(1,:), 'b--')
hold on
plot(times, b_hist(1,:), 'b-')
plot(times, btrue(2,:), 'r--')
plot(times, b_hist(2,:), 'r-')
plot(times, btrue(3,:), 'k--')
plot(times, b_hist(3,:), 'k-')
title('Estimated Gyro Bias vs. True Bias')



%% MEKF functions

function [x_k1, P_k1] = MEKFstep(x_k, P_k, w_k, r_sun_body, r_B_body,...
                                 r_sun_inert, r_B_inert, Q, R, dt)
    [x_pred, P_pred] = predict(x_k, P_k, w_k, dt, Q);
    [z, S, C] = innovation(x_pred, P_pred, r_sun_body, r_B_body, r_sun_inert, r_B_inert, R);
    L = getKalmanGain(P_pred, C, S);
    [dx, x_k1, P_k1] = update(x_pred, P_pred, z, L, C, R);
end

function [x_pred, P_pred] = predict(x_k, P_k, w_k, dt, Q)
    %------------Predict State-----------------
    % Unpack state
    q_k = x_k(1:4);
    b_k = x_k(5:7);

    % "control input"
    theta = norm(w_k-b_k)*dt;
    r = (w_k-b_k)/norm(w_k-b_k);
    s_k = [cos(theta/2); r*sin(theta/2)];   % error quaternion
    
    % propagate quaternion
    q_k1 = quatmult(q_k, s_k);
    
    % update bias
    b_k1 = b_k;
    
    % pack up state
    x_pred = [q_k1; b_k1];
    
    %--------------------Predict Covariance--------------------
%     A_k1 = getAmatrix(s_k, dt);
    A_k = getAmatrix2(s_k, dt);    % try 279C version
    P_pred = A_k*P_k*A_k' + Q;    
end

function [z, S, C] = innovation(x_k, P_k, r_sun_body, r_B_body, r_sun_inert, r_B_inert, R)
    %------------------- Innovation step for state -------------------
    % Get rotation matrix based on quaternion giving rotation from body to
    % inertial reference frame
    q = x_k(1:4);
    R_body2inert =  qtodcm(q);
    
    % Transpose that to get the matrix we actually need for our prefit
    % residuals
    R_N2B = R_body2inert';
    
    % prefit residuals of sun sensor and magnetometer measurements
    z = [r_sun_body; r_B_body] - [R_N2B,    zeros(3);...
                                  zeros(3), R_N2B]* [r_sun_inert; r_B_inert];
    
    % --------------------- Innovation/pre-fit covariance -------------
    C = [skew_mat(R_N2B*r_sun_inert), zeros(3);
         skew_mat(R_N2B*r_B_inert)  , zeros(3)];
    
    S = C*P_k*C' + R;
end

function L = getKalmanGain(P, C, S)
    L = P*C'/S;%inv(S);
end

function [dx, x_k1, P_k1] = update(x_k, P_k, z_k, L, C, R)
    % unpack state
    q_k = x_k(1:4);
    b_k = x_k(5:7);
    
    % post-fit residuals
    dx = L*z_k;
    phi = dx(1:3);
    db = dx(4:6);
    
    dq = [sqrt(1-phi'*phi); phi/2];
     
    % -------------------- State Update -----------------------
    q_k1 = quatmult(q_k, dq);
    b_k1 = b_k + db;
    x_k1 = [q_k1; b_k1];
    
    % -------------------- Covariance Update ------------------
    P_k1 = (eye(6) - L*C) * P_k * (eye(6)-L*C)' + L*R*L';
end


%% Utility functions


% function A = getAmatrix(s_k, dt)
%     L = getL(s_k);
%     R = getR(s_k);
%     V = getV();
%     
%     A = zeros(6,6);
%     A(1:3, 1:3) = V*L'*R*V';
%     A(1:3, 4:6) = -.5*dt*eye(3);
%     A(4:6, 4:6) = eye(3);
% end

function A = getAmatrix2(s_k, dt)
    % version from 279C notes
    A = zeros(6,6);
    A(1:3, 1:3) = qtodcm(s_k)'; % V*L'*R*V';
    A(1:3, 4:6) = -dt*eye(3);
    A(4:6, 4:6) = eye(3);
end

function V = getV()
    % V matrix used to extract vector part of quaternion
    V = zeros(3,4);
    V(:,2:4) = eye(3);
end

function q_out = quatmult(q1, q2)
    % equivalent to quaternion multiplication q1*q2 for scalar first
    % quaternions
    L = getL(q1);
%     R = getR(q1);
    q_out = L*q2;
end

% function rotationMatrix = quat2DCM(q)
%     % returns the rotation matrix from body to inertial frame if a body to
%     % inertial quaternion is given
%     L = getL(q);
%     R = getR(q);
%     
%     temp = L*R';
%     rotationMatrix = temp(2:4, 2:4);    
% end

% Fairly certain Jayden's quat2DCM is backwards
function rotationMatrix = quat2DCM2(q)
    % Try stuff from Jayden's code too
    Q = skew_mat(q(2:4));
    rotationMatrix = ( q(1)^2 - q(2)^2 - q(3)^2 - q(4)^2)*eye(3)...
                     -2*q(1)*Q + 2*q(2:4)*q(2:4)';    
end

function L = getL(q)
    % for left quaternion multiplication q1*q2 = L(q1)q2
    L = zeros(4,4);
    L(1,:) = [q(1), -q(2), -q(3), -q(4)];
    L(:,1) = [q(1);  q(2);  q(3);  q(4)];
    L(2:4, 2:4) = q(1)*eye(3) + skew_mat(q(2:4));
end

function R = getR(q)
    % for right quaternion multiplication q1*q2 = R(q2)*q1
    R = zeros(4,4);
    R(1,:) = [q(1), -q(2), -q(3), -q(4)];
    R(:,1) = [q(1);  q(2);  q(3);  q(4)];
    R(2:4, 2:4) = q(1)*eye(3) - skew_mat(q(2:4));
end

function R = qtodcm(q)

    R = eye(3) + 2*hat(q(2:4))*(hat(q(2:4)) + eye(3)*q(1));
    
end

function [x_skew] = skew_mat(x)
% Returns skew symmetric - cross porduct matrix of a vector

x_skew = [0 -x(3) x(2); x(3) 0 -x(1); -x(2) x(1) 0];

end