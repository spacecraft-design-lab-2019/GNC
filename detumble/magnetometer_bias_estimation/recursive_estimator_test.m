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
numTrials = 100;


for count = 1:numTrials
bias_true = 45000 * randn(3,1);        % true bias, [nanoTesla]
% true_vec(i,:) = bias_true';

[T,scaling_matrix,misalignment_matrix] = get_T_matrix(scaleF,caSense);   % Scaling+Misalignment matrix [-]



    %% Spoof a bunch of measurements
    N = 100;   % Number of measurements
    B_mat = zeros(N,3);
    B_mag_mean = 45000;
    B_mag_std_dev = 5000;
    for i = 1:N
        % TODO: figure out actual distribution of magnitudes in polar orbit
        B_true_mag = B_mag_mean + normrnd(0,B_mag_std_dev);
        B_vec = randn(3,1);
        B_true_vec = B_true_mag * B_vec/norm(B_vec);
        B_mat(i,:) = ( measure(B_true_vec,T,bias_true,cov) )';
    end

    %% Solve using old bias estimator
    [ center, radii, evecs, v, chi2 ] = ellipsoid_fit( B_mat);
    center;
        
    %% Try with new bias estimator
    center_2 = ellipsoid_fit2( B_mat);
    
    %% Calculate error
    
    abs_err_old = (center-bias_true)';
    rel_err_old(count) = norm(abs_err_old)/norm(bias_true);
    
    abs_err_new = (center_2-bias_true)';
    rel_err_new(count) = norm(abs_err_new)/norm(bias_true);
    
end

%% Plot to compare
figure;
plot(rel_err_old)
hold on
plot(rel_err_new)
legend('old','new')



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


