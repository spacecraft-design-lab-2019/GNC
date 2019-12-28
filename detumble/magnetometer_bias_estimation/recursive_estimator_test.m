% Script for testing new iterative least squares method

% Script for simulating bias estimation
% Paul DeTrempe
clear;close all;clc;

%% Magnetometer characteristics
bias_mag = 40000;       % [nanoTesla]
cov = 1000;             % covariance [nanoTesla^2], this number was made up!!!!!!!!!
scaleF = .02;           % scale factor [??]
caSense = .02;          % cross-axis sensitivity [??]


%% Loop this a bunch of times and compare errors:
bias_true = 45000 * randn(3,1);        % true bias, [nanoTesla]
% true_vec(i,:) = bias_true';

[T,scaling_matrix,misalignment_matrix] = get_T_matrix(scaleF,caSense);   % Scaling+Misalignment matrix [-]



%% Spoof a bunch of measurements and record sequential bias estimate
N = 5000;   % Number of measurements
B_mat = zeros(N,3);
B_mag_mean = 45000;
B_mag_std_dev = 5000;

estimates = zeros(N,3);
P = zeros(9,9);
q = zeros(9,1);
for i = 1:N
    % TODO: figure out actual distribution of magnitudes in polar orbit
    B_true_mag = B_mag_mean + normrnd(0,B_mag_std_dev);
    B_vec = randn(3,1);
    B_true_vec = B_true_mag * B_vec/norm(B_vec);
    B_mat(i,:) = ( measure(B_true_vec,T,bias_true,cov) )';
    
    % add to P every time
    x = B_mat(i,1);
    y = B_mat(i,2);
    z = B_mat(i,3);
    a = [ x .* x + y .* y - 2 * z .* z, ...
        x .* x + z .* z - 2 * y .* y, ...
        2 * x .* y, ...
        2 * x .* z, ...
        2 * y .* z, ...
        2 * x, ...
        2 * y, ...
        2 * z, ...
        1 + 0 * x ]';
    y = x^2 + y^2 + z^2;
    
    P = P + a*a';
    q = q + y*a;
    
    % initialize P_inv after some time
    start_idx = 50;
    if i>=start_idx
        P_inv = inv(P);
        [P_inv_updated, q_updated, bias_updated] = UpdateEllipse(P_inv, q, B_mat(i,:)');
        estimates(i,:) = bias_updated';
    end
end



%% Try with new bias estimator
v = ellipsoid_fit2( B_mat);

%% Plot for comparison
figure;
plot(abs(estimates(start_idx:end,1)-bias_true(1)));
hold on
plot(abs(estimates(start_idx:end,2)-bias_true(2)));
plot(abs(estimates(start_idx:end,3)-bias_true(3)));


%% Functions
function [T,scaling_matrix, misalignment_matrix] = get_T_matrix(scaleF,caSense)
scaling_matrix = eye(3) + diag(normrnd(0,scaleF,[3,1]));
misalignment_matrix = normrnd(0,caSense,[3,3]);
% make skew symmetric
misalignment_matrix(2:3,1) = -misalignment_matrix(1,2:3)';
misalignment_matrix(3,2) = -misalignment_matrix(2,3);
misalignment_matrix = misalignment_matrix - diag(diag(misalignment_matrix)); % zero out items on diagonal
T = scaling_matrix + misalignment_matrix;
end

function B_measured = measure(B_true,T,bias,covariance)
multiplicative_noise = T;
additive_noise = normrnd(0,covariance,[3,1]);
B_measured = multiplicative_noise*B_true + additive_noise + bias;
end


