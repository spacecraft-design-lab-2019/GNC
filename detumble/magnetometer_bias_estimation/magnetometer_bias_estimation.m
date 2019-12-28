% Script for simulating bias estimation
% Paul DeTrempe
clear;close all;clc;

%% Magnetometer characteristics
bias_mag = 40000;       % [nanoTesla]
cov = 1000;             % covariance [nanoTesla^2], this number was made up!!!!!!!!!
scaleF = .02;          % scale factor [??]
caSense = .02;        % cross-axis sensitivity [??]
bias_true = 45000 * randn(3,1);        % true bias, [nanoTesla]

[T,scaling_matrix,misalignment_matrix] = get_T_matrix(scaleF,caSense);   % Scaling+Misalignment matrix [-]

%% Spoof a bunch of measurements
N = 1000;   % Number of measurements
B_mat = zeros(N,3);
B_mag_mean = 45000
B_mag_std_dev = 5000
for i = 1:N
    % TODO: figure out actual distribution of magnitudes in polar orbit
    B_true_mag = B_mag_mean + normrnd(0,B_mag_std_dev);
    B_vec = randn(3,1);
    B_true_vec = B_true_mag * B_vec/norm(B_vec);
    B_mat(i,:) = ( measure(B_true_vec,T,bias_true,cov) )';
end

%% Perform fitting using CVX
% cvx_begin SDP quiet
% 
% % Primal problem
% variable A(3,3) symmetric
% variable b(3)
% 
% maximize det_rootn(A)
% 
% subject to
% for i = 1:N
%     norm(A*(B_mat(i,:))'-b) <= 1;
% end
% 
% A>0;
% 
% cvx_end
% 
% A
% b
% 
% % find matrix inverse for mapping from unit ball
% [V,D] = eig(A);
% 
% A_inv = V*(diag(diag(D).^-1))*V';

%% Solve using old bias estimator
[ center, radii, evecs, v, ~ ] = ellipsoid_fit( B_mat);

%% Try with new bias estimator
% u = ellipsoid_fit2( B_mat);

%% Output Error Metrics
% bias_estimated_sdp = A_inv*b;
% abs_bias_error_sdp = (bias_estimated_sdp-bias_true)';
% rel_bias_error_sdp = norm(abs_bias_error_sdp)/norm(bias_true)

bias_estimated_lls = center;       % nanoTesla
abs_bias_error_lls = (bias_estimated_lls-bias_true)';
rel_bias_error_lls = norm(abs_bias_error_lls)/norm(bias_true);

% Unpack as a matrix centered at
[X_unsc,Y_unsc,Z_unsc] = ellipsoid(0,0,0,1,1,1);

%% Manipulate SDP data for plotting


% 
% X_sc = zeros(size(X_unsc));
% Y_sc = zeros(size(Y_unsc));
% Z_sc = zeros(size(Z_unsc));
% 
% p = size(X_sc,1);   % number of points for plotting
% 
% for i = 1:p
%     for j = 1:p
%         point_unsc = [X_unsc(i,j);Y_unsc(i,j);Z_unsc(i,j)];
%         point_sc = A_inv*(point_unsc+b);
%         X_sc(i,j) = point_sc(1);
%         Y_sc(i,j) = point_sc(2);
%         Z_sc(i,j) = point_sc(3);
%     end
% end

%% Manipulate LLS data for plotting
[X,Y,Z] = ellipsoid(0,0,0,radii(1),radii(2),radii(3));

q = size(X,1);
for i = 1:q
    for j = 1:q
        point_unrotated = [X(i,j);Y(i,j);Z(i,j)];
        point_rotated = evecs*point_unrotated;
        X_rot(i,j) = point_rotated(1)+center(1);
        Y_rot(i,j) = point_rotated(2)+center(2);
        Z_rot(i,j) = point_rotated(3)+center(3);
    end
end

%% Manipulate true T_matrix data for plotting
for i = 1:q
    for j = 1:q
        point_unsc = [X_unsc(i,j);Y_unsc(i,j);Z_unsc(i,j)];
        point_sc = B_mag_mean*T*(point_unsc)+bias_true;
        X_T(i,j) = point_sc(1);
        Y_T(i,j) = point_sc(2);
        Z_T(i,j) = point_sc(3);
    end
end

%% Compare true T matrix to estimate
T
evecs
radii

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
surf(X_rot,Y_rot,Z_rot,'FaceAlpha',0.5,'EdgeColor','none','FaceColor','b')
% SDP fit
% surf(X_sc,Y_sc,Z_sc,'FaceAlpha',0.4,'EdgeColor','none','FaceColor','c')
% T matrix scaled up by mean magnitude
surf(X_T,Y_T,Z_T,'FaceAlpha',.2,'EdgeColor','none','FaceColor','r')
% Biases
plot3(0,0,0,'ko','Linewidth',5)
plot3([0,bias_true(1)],[0,bias_true(2)],[0,bias_true(3)],'r','Linewidth',2)
plot3([0,bias_estimated_lls(1)],[0,bias_estimated_lls(2)],[0,bias_estimated_lls(3)],'g--','Linewidth',2)
% plot3([0,bias_estimated_sdp(1)],[0,bias_estimated_sdp(2)],[0,bias_estimated_sdp(3)],'k--','Linewidth',2)
legend('Measurements','Least-squares Fit','Minimum Bounding Ellipsoid','Scaled T matrix','Origin','True bias','Estimated bias (LLS)','Estimated bias (SDP)')

% surf(X_T,Y_T,Z_T,'FaceAlpha',.2,'EdgeColor','none','FaceColor','y')


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


