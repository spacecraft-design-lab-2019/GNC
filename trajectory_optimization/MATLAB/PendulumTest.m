close all;
clear all;

Qf = 100*eye(2);
Q = 0.01*eye(2);

R = .3;

x0 = [0 0]';
xg = [pi 0]';
u0 = zeros(1, 399);

% Using MATLAB
[xhist, uhist, K] = iLQRv1(x0, xg, u0, Q, R, Qf, .01, 1e-3);

figure(1);
subplot(3,1,1)
plot(xhist(1,:));
ylabel('q');

subplot(3,1,2)
plot(xhist(2,:));
ylabel('qdot');

subplot(3,1,3);
plot(uhist);
ylabel('u');
% 
% figure(2);
% semilogy(J);
% ylabel('Cost');
% xlabel('Iteration');
% 
% 
% % Compare to CPP iLQR
% % Load file output from cpp iLQR
% cpp_data = readmatrix("iLQR_pendulum_data.csv");
% xtraj = cpp_data(:, 1:2);
% utraj = cpp_data(:, 3);
% Jhist = cpp_data(:, 4);
% 
% % Plot cpp iLQR results
% figure(3);
% sgtitle('cpp iLQR results');
% subplot(3,1,1)
% plot(xtraj(:, 1));
% ylabel('q');
% 
% subplot(3,1,2)
% plot(xtraj(:, 2));
% ylabel('qdot');
% 
% subplot(3,1,3)
% plot(utraj);
% ylabel('u');
% 
% figure(4);
% title('cpp iLQR cost');
% semilogy(Jhist);
% ylabel('Cost');
% 
% figure(5);
% subplot(1,2,1)
% sgtitle('MATLAB iLQR nominal trajectory vs CPP')
% plot(xhist(1, :));
% ylabel('Angle, \phi [rad]');
% xlabel('step k')
% 
% subplot(1,2,2)
% plot(xtraj(:, 1));
% xlabel('step k')