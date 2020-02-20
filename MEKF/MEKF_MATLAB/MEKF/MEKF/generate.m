clear all

dt = .1;
J = diag([1.1 1.2 2]);

%Noise on vector sensors
V1 = ((1*pi/180)^2)*eye(3); %sun sensor
V2 = ((5*pi/180)^2)*eye(3); %magnetometer
V3 = ((.01*pi/180)^2)*eye(3); %star tracker

%Noise on gyro
Vw = ((.01*pi/180)^2)*eye(3); %gyro measurement noise
Vb = ((.001*pi/180)^2)*eye(3); %bias noise

%Random external torque
tau = 0.001*randn(3,1);

q0 = randn(4,1);
q0 = q0/norm(q0);

w0 = (5*pi/180)*randn(3,1); %~5 deg/sec initial rate
b0 = (10*pi/180)*randn(3,1); %~5 deg/sec bias

%Simulate
x0 = [q0; w0];
[t,xhist] = ode45(@(t,x)dynamics(t,x,J,tau),0:dt:150,x0);
xhist = xhist';

%Generate vector sensor measurements
[rN1, rN2, rB1hist, rB2hist, qSThist, rB1true, rB2true] = sensors(xhist,V1,V2,V3);

%Generate gyro data
[whist,bhist] = gyro(xhist,b0,Vw,Vb);

close all

figure(1)
subplot(3,1,1);
plot(xhist(5,:));
hold on
plot(whist(1,:));
title('Gyro');
subplot(3,1,2);
plot(xhist(6,:));
hold on
plot(whist(2,:));
subplot(3,1,3);
plot(xhist(7,:));
hold on
plot(whist(3,:));

figure(2);
subplot(4,1,1);
plot(xhist(1,:));
subplot(4,1,2);
plot(xhist(2,:));
subplot(4,1,3);
plot(xhist(3,:));
subplot(4,1,4);
plot(xhist(4,:));

figure(3);
subplot(2,1,1);
plot(rB1hist(1,:));
hold on
plot(rB1hist(2,:));
plot(rB1hist(3,:));
subplot(2,1,2);
plot(rB2hist(1,:));
hold on
plot(rB2hist(2,:));
plot(rB2hist(3,:));

V = blkdiag(V3,V1,V2);
W = blkdiag(Vw*dt*dt,Vb);
qtrue = xhist(1:4,:);
btrue = bhist;

save('mekf_inputs','dt','rN1','rN2','rB1hist','rB2hist','qSThist','whist','V','W');
save('mekf_truth','qtrue','btrue');