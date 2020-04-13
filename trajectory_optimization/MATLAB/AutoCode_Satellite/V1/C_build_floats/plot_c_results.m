% Script for plotting the results from auto-coded ilqr

clear;
close all;
clc;


% Read the data from the csv
c_data = readmatrix("c_milqr_run1.csv");
q = c_data(:,2:5);  % quaternion
om = c_data(:,6:8);  % omega


% Plot results
figure(1)

% Quaternion
subplot(3,1,1)
plot(q(:,1))
hold on
plot(q(:,2))
plot(q(:,3))
plot(q(:,4))
legend('q1','q2','q3','q4');

% Angular velocities
subplot(3,1,2)
plot(om(:,1))
hold on
plot(om(:,2))
plot(om(:,3))
legend('omega1','omega2','omega3');
