% Prototype script for Matlab MEKF
% To be autocoded into C
% Paul DeTrempe
clear; close all; clc;

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
load('mekf_inputs.mat')
load('mekf_truth.mat')
Q = W;  % Process Noise
R = V;  % Measurement Noise
N = size(qtrue, 2);
times = 0:dt:(N-1)*dt;
% load('output_mekf.mat')

%% Run the MEKF
% Initial conditions
x_k = [qtrue(:,1); btrue(:,1)];
P_k = Q;
w_k = whist(:,1);
r_sun_inert = rN1;
r_B_inert = rN2;

% Preallocate
b_hist = zeros(size(btrue));
q_hist = zeros(size(qtrue));
P_hist = zeros(size(P_k,1), N);

for i = 1:N
    % measure some stuff
    w_k = whist(:,i);
    r_sun_body = rB1hist(:,i);
    r_B_body = rB2hist(:,i);
    
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
plot(times, qtrue(1,:))
hold on
plot(times, q_hist(1,:));



%% MEKF functions

function [x_k1, P_k1] = MEKFstep(x_k, P_k, w_k, r_sun_body, r_B_body,...
                                 r_sun_inert, r_B_inert, Q, R, dt)
    [x_pred, P_pred] = predict(x_k, P_k, w_k, dt, Q);
    [z, S, C] = innovation(x_pred, P_pred, r_sun_body, r_B_body, r_sun_inert, r_B_inert, R);
    L = getKalmanGain(P_pred, C, S);
    [dx, x_k1, P_k1] = update(x_k, P_k, z, L, C, R);
end

function [x_pred, P_pred] = predict(x_k, P_k, w_k, dt, Q)
    %------------Predict State-----------------
    % Unpack state
    q_k = x_k(1:4);
    b_k = x_k(5:7);

    % "control input"
    u_k = w_k + b_k;
    theta = norm(u_k-b_k)*dt;
    r = (u_k-b_k)/norm(u_k-b_k);
    s_k = [cos(theta/2); r*sin(theta/2)];   % error quaternion
    
    % propagate quaternion
    q_k1 = quatmult(q_k, s_k);
    
    % update bias
    b_k1 = b_k;
    
    % pack up state
    x_pred = [q_k1; b_k1];
    
    %--------------------Predict Covariance--------------------
    A_k = getAmatrix(s_k, dt);
    P_pred = A_k*P_k*A_k' + Q;    
end

function [z, S, C] = innovation(x_k, P_k, r_sun_body, r_B_body, r_sun_inert, r_B_inert, R)
    %------------------- Innovation step for state -------------------
    % Get rotation matrix based on quaternion giving rotation from body to
    % inertial reference frame
    q = x_k(1:4);
    R_body2inert = quat2DCM(q);
    
    % Transpose that to get the matrix we actually need for our prefit
    % residuals
    R_N2B = R_body2inert';
    
    % prefit residuals of sun sensor and magnetometer measurements
    z = [r_sun_body; r_B_body] - [R_N2B,    zeros(3);...
                                  zeros(3), R_N2B]* [r_sun_inert; r_B_inert];
    
    % --------------------- Innovation/pre-fit covariance -------------
    C = [2*skew_mat(r_sun_body), zeros(3);
         2*skew_mat(r_B_body)  , zeros(3)];
    
    S = C*P_k*C' + R;
end

function L = getKalmanGain(P, C, S)
    L = P*C'*inv(S);
end

function [dx, x_k1, P_k1] = update(x_k, P_k, z_k, L, C, R)
    % unpack state
    q_k = x_k(1:4);
    b_k = x_k(5:7);
    
    % post-fit residuals
    dx = L*z_k;
    phi = dx(1:3);
    db = dx(4:6);
    
    dq = [sqrt(1-phi'*phi); phi];
    
    % -------------------- State Update -----------------------
    q_k1 = quatmult(q_k, dq);
    b_k1 = b_k + db;
    x_k1 = [q_k1; b_k1];
    
    % -------------------- Covariance Update ------------------
    P_k1 = (eye(6) - L*C) * P_k * (eye(6)-L*C)' + L*R*L';
end


%% Utility functions


function A = getAmatrix(s_k, dt)
    L = getL(s_k);
    R = getR(s_k);
    V = getV();
    
    A = zeros(6,6);
    A(1:3, 1:3) = V*L'*R*V';
    A(1:3, 4:6) = .5*dt*eye(3);
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
    R = getR(q1);
    q_out = L*R'*q2;
end

function rotationMatrix = quat2DCM(q)
    % returns the rotation matrix from body to inertial frame if a body to
    % inertial quaternion is given
    L = getL(q);
    R = getR(q);
    
    temp = L*R';
    rotationMatrix = temp(2:4, 2:4);    
end

% Fairly certain Jayden's quat2DCM is backwards
% function rotationMatrix = quat2DCM2(q)
%     % Try stuff from Jayden's code too
%     Q = skew_mat(q(2:4));
%     rotationMatrix = ( q(1)^2 - q(2)^2 - q(3)^2 - q(4)^2)*eye(3)...
%                      -2*q(1)*Q + 2*q(2:4)*q(2:4)';    
% end

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

function [x_skew] = skew_mat(x)
% Returns skew symmetric - cross porduct matrix of a vector

x_skew = [0 -x(3) x(2); x(3) 0 -x(1); -x(2) x(1) 0];

end