% Script for simulating bias estimation
% Paul DeTrempe
clear;close all;clc;

%% Magnetometer characteristics
bias_mag = 40000;       % [nanoTesla]
cov = 1000;             % covariance [nanoTesla^2], this number was made up!!!!!!!!!
scaleF = .02;           % scale factor [??]
caSense = .02;          % cross-axis sensitivity [??]
bias_true = 45000 * randn(3,1)        % true bias, [nanoTesla]

[T,scaling_matrix,misalignment_matrix] = get_T_matrix(scaleF,caSense);   % Scaling+Misalignment matrix [-]

%% Spoof a bunch of measurements
N = 100;   % Number of measurements
B_mat = zeros(N,3);
B_mag_mean = 45000;
B_mag_std_dev = 5000;
for i = 1:N
    % TODO: figure out actual distribution of magnitudes in polar orbit
    B_true_mag = B_mag_mean; % + 5000*randn();
    B_vec = randn(3,1);
    B_true_vec = B_true_mag * B_vec/norm(B_vec);
    B_mat(i,:) = ( measure(B_true_vec,T,bias_true,cov) )';
end

%% Perform fitting using CVX
cvx_begin SDP

% Primal problem
variable A(3,3) symmetric
variable b(3)

maximize det_rootn(A)

subject to
for i = 1:N
    norm(A*(B_mat(i,:))'-b)<=1;
end

A>0

cvx_end

A
b

% find matrix square root for mapping from unit ball
[V,D] = eig(A);

A_inv = V*(diag(diag(D).^-1))*V';




%% Perform fitting using linear least squares
[ center, radii, evecs, v, chi2 ] = ellipsoid_fit(B_mat, '' );

%% Estimate T matrix from each ellipsoid
[U_sdp,S_sdp,V_sdp] = svd(A_inv)

%% Output Error Metrics
bias_estimated_sdp = A_inv*b
abs_bias_error_sdp = (bias_estimated_sdp-bias_true)'
rel_bias_error_sdp = norm(abs_bias_error_sdp)/norm(bias_true)

bias_estimated_lls = center         % nanoTesla
abs_bias_error_lls = (bias_estimated_lls-bias_true)'
rel_bias_error_lls = norm(abs_bias_error_lls)/norm(bias_true)

%% Manipulate SDP data for plotting

% Unpack as a matrix centered at
[X_unsc,Y_unsc,Z_unsc] = ellipsoid(0,0,0,1,1,1);

X_sc = zeros(size(X_unsc));
Y_sc = zeros(size(Y_unsc));
Z_sc = zeros(size(Z_unsc));

p = size(X_sc,1);   % number of points for plotting

for i = 1:p
    for j = 1:p
        point_unsc = [X_unsc(i,j);Y_unsc(i,j);Z_unsc(i,j)];
        point_sc = A_inv*(point_unsc+b);
        X_sc(i,j) = point_sc(1);
        Y_sc(i,j) = point_sc(2);
        Z_sc(i,j) = point_sc(3);
    end
end

%% Plot
figure;
% raw data
plot3(B_mat(:,1),B_mat(:,2),B_mat(:,3),'bo')
axis equal
grid on
hold on
xlabel('X, [nT]')
ylabel('Y, [nT]')
zlabel('Z, [nT]')
% least squares fit
[X,Y,Z] = ellipsoid(center(1),center(2),center(3),radii(1),radii(2),radii(3));
surf(X,Y,Z,'FaceAlpha',0.5,'EdgeColor','none','FaceColor','b')
% SDP fit
surf(X_sc,Y_sc,Z_sc,'FaceAlpha',0.4,'EdgeColor','none','FaceColor','c')
% Biases
plot3(0,0,0,'ko','Linewidth',5)
plot3([0,bias_true(1)],[0,bias_true(2)],[0,bias_true(3)],'r','Linewidth',2)
plot3([0,bias_estimated_lls(1)],[0,bias_estimated_lls(2)],[0,bias_estimated_lls(3)],'g--','Linewidth',2)
plot3([0,bias_estimated_sdp(1)],[0,bias_estimated_sdp(2)],[0,bias_estimated_sdp(3)],'k--','Linewidth',2)
legend('Measurements','Least-squares Fit','Minimum Bounding Ellipsoid','Origin','True bias','Estimated bias (lls)','Estimate bias (sdp)')



%% Functions
function [T,scaling_matrix, misalignment_matrix] = get_T_matrix(scaleF,caSense)
    scaling_matrix = eye(3) + scaleF*diag(randn(3,1));
    misalignment_matrix = mvnrnd(zeros(1,3),caSense*eye(3),1);
    misalignment_matrix = misalignment_matrix - diag(diag(misalignment_matrix)); % zero out items on diagonal
    T = scaling_matrix + misalignment_matrix;
end

function B_measured = measure(B_true,T,bias,covariance)
    multiplicative_noise = T;
    additive_noise = covariance * randn(3,1);
    B_measured = multiplicative_noise*B_true + additive_noise + bias;
end

