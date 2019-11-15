% Script for simulating bias estimation
% Paul DeTrempe
clear;close all;clc;

%% Magnetometer characteristics
bias_mag = 40000;       % [nanoTesla]
cov = 1000;             % covariance [nanoTesla^2], this number was made up!!!!!!!!!
scaleF = .02;           % scale factor [??]
caSense = .02;          % cross-axis sensitivity [??]
bias_true = 45000 * randn(3,1)        % true bias, [nanoTesla]

T = get_T_matrix(scaleF,caSense);   % Scaling+Misalignment matrix [-]

%% Spoof a bunch of measurements
N = 1000;   % Number of measurements
B_mat = zeros(N,3);
for i = 1:N
    B_true_mag = 45000 + 20000*randn();
    B_vec = randn(3,1);
    B_true_vec = B_true_mag * B_vec/norm(B_vec);
    B_mat(i,:) = ( measure(B_true_vec,T,bias_true,cov) )';
end

%% Perform fitting using CVX
% cvx_begin SDP
% 
% % Primal problem
% variable A(3,3) symmetric nonnegative
% variable b(3)
% 
% minimize log(det_inv(A))
% 
% subject to
% for i = 1:N
%     norm(A*(B_mat(i,:))'-b)<=1;
% end
% 
% A>0
% 
% cvx_end
% 
% A
% b

%% Perform fitting using linear least squares
[ center, radii, evecs, v, chi2 ] = ellipsoid_fit(B_mat, '' );
bias_estimated = center         % nanoTesla
abs_bias_error = (bias_estimated-bias_true)'
rel_bias_error = norm(abs_bias_error)/norm(bias_true)

%% Plot
figure;
% raw data
plot3(B_mat(:,1),B_mat(:,2),B_mat(:,3),'bo')
axis equal
grid on
hold on
xlabel('X, [nT]')
ylabel('X, [nT]')
zlabel('X, [nT]')
% least squares fit
[X,Y,Z] = ellipsoid(center(1),center(2),center(3),radii(1),radii(2),radii(3));
surf(X,Y,Z,'FaceAlpha',0.5,'EdgeColor','none','FaceColor','c')
% Biases
plot3(0,0,0,'ko','Linewidth',5)
plot3([0,bias_true(1)],[0,bias_true(2)],[0,bias_true(3)],'r--','Linewidth',3)
plot3([0,bias_estimated(1)],[0,bias_estimated(2)],[0,bias_estimated(3)],'k--','Linewidth',3)
legend('Measurements','Least-squares Fit','Origin','True bias','Estimated bias')



%% Functions
function T = get_T_matrix(scaleF,caSense)
    scaling_matrix = eye(3) + diag(wgn(3,1,scaleF,'linear'));
    misalignment_matrix = mvnrnd(zeros(1,3),caSense*eye(3),1);
    misalignment_matrix = misalignment_matrix - diag(diag(misalignment_matrix)); % zero out items on diagonal
    T = scaling_matrix + misalignment_matrix;
end

function B_measured = measure(B_true,T,bias,covariance)
    multiplicative_noise = T;
    additive_noise = covariance * randn(3,1);
    B_measured = multiplicative_noise*B_true + additive_noise + bias;
end

