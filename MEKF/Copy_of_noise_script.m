%% Simulation

clear

%noise covariances 
omega_process_sigma = .000001*eye(3);
quat_process_sigma = .000001*eye(4);
omega_sensor_sigma = .0001*eye(3);
quat_sensor_sigma = .0001*eye(4);

%process and sensor noise matricees 
Q_process(1:3,1:3) = omega_process_sigma;
Q_process(4:7,4:7) = quat_process_sigma;
R_sensor(1:3,1:3) = omega_sensor_sigma;
R_sensor(4:7,4:7) = quat_sensor_sigma;

%initial conditions 
omega0 = [deg2rad(2);deg2rad(4); deg2rad(8)];
quat0 = [0 0 0 1]';
dt = .001;

%run the sim in simulink
[sim_results]=sim('noise_euler_eqn.slx')

%index the simulation results
omega = sim_results.omega;
t_vec = sim_results.tout;
STM_omega = sim_results.STM_omega;
quaternion = sim_results.quaternion;
STM_quaternion = sim_results.STM_quaternion;
omega_estimated = sim_results.omega_estimated;
quat_estimated = sim_results.quat_estimated;


%% EKF Implementation 

Ix = 197.28;
Iy = 222.835;
Iz = 277.6;

%set inertia tensor
I = [Ix 0 0; 0 Iy 0; 0 0 Iz];
dt = .001;

EKF_mu(1,:)=[omega0;quat0]';
EKF_P{1} = eye(7);
EKF_mu2(1,:)=[omega0;quat0]';
EKF_P2{1} = eye(7);

MEKF_mu(1,:)=[omega0;quat0]';
MEKF_P{1} = eye(7);
MEKF_mu2(1,:)=[omega0;quat0]';
MEKF_P2{1} = eye(7);

UKF_mu(:,1)=[omega0;quat0]';
UKF_P{1} = eye(7);
UKF_mu(:,1)=[omega0;quat0]';
UKF_P2{1} = eye(7);

%sensitivity Matrix
H = eye(7);
i = 0;

%first measurement 
 y(i+1,1:7) = [omega_estimated(i+1,:)';quat_estimated(i+1,:)']';
 
 tic
for i =1:(length(quat_estimated)-1)

    %here is measurement 
    y(i+1,1:7) = [omega_estimated(i+1,:)';quat_estimated(i+1,:)']';
    
    %EKF
    [EKF_mu(i+1,:),EKF_P{i+1}] = EKF(EKF_mu(i,:),EKF_P{i},omega_estimated(i+1,:),quat_estimated(i+1,:),Q_process,R_sensor);

end
toc

tic
for i =1:(length(quat_estimated)-1)

    %here is measurement 
    y(i+1,1:7) = [omega_estimated(i+1,:)';quat_estimated(i+1,:)']';
    
    %EKF
    [MEKF_mu(i+1,:),MEKF_P{i+1}] = MEKF(MEKF_mu(i,:),MEKF_P{i},omega_estimated(i+1,:),quat_estimated(i+1,:),Q_process,R_sensor);

end
toc

tic
y_in = [omega_estimated';quat_estimated'];
for ii = 1:length(quat_estimated)-1
    [UKF_mu(:,ii+1), UKF_P{ii+1}] = UKF(UKF_mu(:,ii),UKF_P{ii},Q_process,R_sensor,y_in(:,ii),dt);
end
toc
UKF_mu = UKF_mu';


%% plot errors EKF

% quat error 
figure
hold on 
sgtitle('Quaternion True State')
subplot(4,1,1)
hold on
plot(t_vec,quaternion(:,1))

ylabel('q_1')

subplot(4,1,2)
hold on
plot(t_vec,quaternion(:,2))

ylabel('q_2')

subplot(4,1,3)
hold on
plot(t_vec,quaternion(:,3))

ylabel('q_3')

subplot(4,1,4)
hold on
plot(t_vec,quaternion(:,4))

ylabel('q_4')
xlabel('Time (s)')

hold off

% quat error 
figure
hold on 
sgtitle('Angular Velocity True State')
subplot(3,1,1)
hold on
plot(t_vec,omega(:,1))

ylabel('\omega_x')

subplot(3,1,2)
hold on
plot(t_vec,omega(:,2))

ylabel('\omega_y')

subplot(3,1,3)
hold on
plot(t_vec,omega(:,3))

ylabel('\omega_z')

xlabel('Time (s)')

hold off

%plot EKF and Measurement Errors
figure
hold on
sgtitle('EKF vs Measurement Errors for \omega')
subplot(3,1,1)
hold on 
plot(t_vec,abs(omega(:,1) - y(:,1)))
plot(t_vec,abs(omega(:,1) - EKF_mu(:,1)))
ylabel('\Delta \omega_x (rad/s)')
legend('Measurement Error','EKF Error')

subplot(3,1,2)
hold on 
plot(t_vec,abs(omega(:,2) - y(:,2)))
plot(t_vec,abs(omega(:,2) - EKF_mu(:,2)))
ylabel('\Delta \omega_y (rad/s)')
legend('Measurement Error','EKF Error')

subplot(3,1,3)
hold on 
plot(t_vec,abs(omega(:,2) - y(:,2)))
plot(t_vec,abs(omega(:,2) - EKF_mu(:,2)))
ylabel('\Delta \omega_z (rad/s)')
legend('Measurement Error','EKF Error')
xlabel('Time (s)')



figure
hold on
sgtitle('EKF vs Measurement Errors for Quaternion')
subplot(4,1,1)
hold on 
plot(t_vec,abs(quaternion(:,1) - y(:,4)))
plot(t_vec,abs(quaternion(:,1) - EKF_mu(:,4)))
ylabel('\Delta q_1')
legend('Measurement Error','EKF Error')

subplot(4,1,2)
hold on 
plot(t_vec,abs(quaternion(:,2) - y(:,5)))
plot(t_vec,abs(quaternion(:,2) - EKF_mu(:,5)))
ylabel('\Delta q_2')
legend('Measurement Error','EKF Error')


subplot(4,1,3)
hold on 
plot(t_vec,abs(quaternion(:,3) - y(:,6)))
plot(t_vec,abs(quaternion(:,3) - EKF_mu(:,6)))
ylabel('\Delta q_3')
legend('Measurement Error','EKF Error')

subplot(4,1,4)
hold on 
plot(t_vec,abs(quaternion(:,4) - y(:,7)))
plot(t_vec,abs(quaternion(:,4) - EKF_mu(:,7)))
ylabel('\Delta q_4')
legend('Measurement Error','EKF Error')
xlabel('Time (s)')

%% plot errors MEKF

% quat error 
figure
hold on 
sgtitle('Quaternion True State')
subplot(4,1,1)
hold on
plot(t_vec,quaternion(:,1))

ylabel('q_1')

subplot(4,1,2)
hold on
plot(t_vec,quaternion(:,2))

ylabel('q_2')

subplot(4,1,3)
hold on
plot(t_vec,quaternion(:,3))

ylabel('q_3')

subplot(4,1,4)
hold on
plot(t_vec,quaternion(:,4))

ylabel('q_4')
xlabel('Time (s)')

hold off

% quat error 
figure
hold on 
sgtitle('Angular Velocity True State')
subplot(3,1,1)
hold on
plot(t_vec,omega(:,1))

ylabel('\omega_x')

subplot(3,1,2)
hold on
plot(t_vec,omega(:,2))

ylabel('\omega_y')

subplot(3,1,3)
hold on
plot(t_vec,omega(:,3))

ylabel('\omega_z')

xlabel('Time (s)')

hold off

%plot MEKF and Measurement Errors
figure
hold on
sgtitle('MEKF vs Measurement Errors for \omega')
subplot(3,1,1)
hold on 
plot(t_vec,abs(omega(:,1) - y(:,1)))
plot(t_vec,abs(omega(:,1) - MEKF_mu(:,1)))
ylabel('\Delta \omega_x (rad/s)')
legend('Measurement Error','MEKF Error')

subplot(3,1,2)
hold on 
plot(t_vec,abs(omega(:,2) - y(:,2)))
plot(t_vec,abs(omega(:,2) - MEKF_mu(:,2)))
ylabel('\Delta \omega_y (rad/s)')
legend('Measurement Error','MEKF Error')

subplot(3,1,3)
hold on 
plot(t_vec,abs(omega(:,2) - y(:,2)))
plot(t_vec,abs(omega(:,2) - MEKF_mu(:,2)))
ylabel('\Delta \omega_z (rad/s)')
legend('Measurement Error','MEKF Error')
xlabel('Time (s)')



figure
hold on
sgtitle('MEKF vs Measurement Errors for Quaternion')
subplot(4,1,1)
hold on 
plot(t_vec,abs(quaternion(:,1) - y(:,4)))
plot(t_vec,abs(quaternion(:,1) - MEKF_mu(:,4)))
ylabel('\Delta q_1')
legend('Measurement Error','MEKF Error')

subplot(4,1,2)
hold on 
plot(t_vec,abs(quaternion(:,2) - y(:,5)))
plot(t_vec,abs(quaternion(:,2) - MEKF_mu(:,5)))
ylabel('\Delta q_2')
legend('Measurement Error','MEKF Error')


subplot(4,1,3)
hold on 
plot(t_vec,abs(quaternion(:,3) - y(:,6)))
plot(t_vec,abs(quaternion(:,3) - MEKF_mu(:,6)))
ylabel('\Delta q_3')
legend('Measurement Error','MEKF Error')

subplot(4,1,4)
hold on 
plot(t_vec,abs(quaternion(:,4) - y(:,7)))
plot(t_vec,abs(quaternion(:,4) - MEKF_mu(:,7)))
ylabel('\Delta q_4')
legend('Measurement Error','MEKF Error')
xlabel('Time (s)')

%% plot errors UKF

% quat error 
figure
hold on 
sgtitle('Quaternion True State')
subplot(4,1,1)
hold on
plot(t_vec,quaternion(:,1))

ylabel('q_1')

subplot(4,1,2)
hold on
plot(t_vec,quaternion(:,2))

ylabel('q_2')

subplot(4,1,3)
hold on
plot(t_vec,quaternion(:,3))

ylabel('q_3')

subplot(4,1,4)
hold on
plot(t_vec,quaternion(:,4))

ylabel('q_4')
xlabel('Time (s)')

hold off

% quat error 
figure
hold on 
sgtitle('Angular Velocity True State')
subplot(3,1,1)
hold on
plot(t_vec,omega(:,1))

ylabel('\omega_x')

subplot(3,1,2)
hold on
plot(t_vec,omega(:,2))

ylabel('\omega_y')

subplot(3,1,3)
hold on
plot(t_vec,omega(:,3))

ylabel('\omega_z')

xlabel('Time (s)')

hold off

%plot UKF and Measurement Errors
figure
hold on
sgtitle('UKF vs Measurement Errors for \omega')
subplot(3,1,1)
hold on 
plot(t_vec,abs(omega(:,1) - y(:,1)))
plot(t_vec,abs(omega(:,1) - UKF_mu(:,1)))
ylabel('\Delta \omega_x (rad/s)')
legend('Measurement Error','UKF Error')

subplot(3,1,2)
hold on 
plot(t_vec,abs(omega(:,2) - y(:,2)))
plot(t_vec,abs(omega(:,2) - UKF_mu(:,2)))
ylabel('\Delta \omega_y (rad/s)')
legend('Measurement Error','UKF Error')

subplot(3,1,3)
hold on 
plot(t_vec,abs(omega(:,2) - y(:,2)))
plot(t_vec,abs(omega(:,2) - UKF_mu(:,2)))
ylabel('\Delta \omega_z (rad/s)')
legend('Measurement Error','UKF Error')
xlabel('Time (s)')



figure
hold on
sgtitle('UKF vs Measurement Errors for Quaternion')
subplot(4,1,1)
hold on 
plot(t_vec,abs(quaternion(:,1) - y(:,4)))
plot(t_vec,abs(quaternion(:,1) - UKF_mu(:,4)))
ylabel('\Delta q_1')
legend('Measurement Error','UKF Error')

subplot(4,1,2)
hold on 
plot(t_vec,abs(quaternion(:,2) - y(:,5)))
plot(t_vec,abs(quaternion(:,2) - UKF_mu(:,5)))
ylabel('\Delta q_2')
legend('Measurement Error','UKF Error')


subplot(4,1,3)
hold on 
plot(t_vec,abs(quaternion(:,3) - y(:,6)))
plot(t_vec,abs(quaternion(:,3) - UKF_mu(:,6)))
ylabel('\Delta q_3')
legend('Measurement Error','UKF Error')

subplot(4,1,4)
hold on 
plot(t_vec,abs(quaternion(:,4) - y(:,7)))
plot(t_vec,abs(quaternion(:,4) - UKF_mu(:,7)))
ylabel('\Delta q_4')
legend('Measurement Error','UKF Error')
xlabel('Time (s)')


%% plot errors together

error_norm_quat_EKF = zeros(size(EKF_mu,1),1);
for i = 1:size(EKF_mu,1)
    error_norm_quat_EKF(i) = quat2angle(quaternion(i,:))-quat2angle(EKF_mu(i,4:7));
end
RMSE_EKF = sqrt(sum(error_norm_quat_EKF.^2)/t_vec(end));

error_norm_quat_MEKF = zeros(size(MEKF_mu,1),1);
for i = 1:size(MEKF_mu,1)
    error_norm_quat_MEKF(i) = norm(quaternion(i,:)-MEKF_mu(i,4:7));
end
RMSE_MEKF = sqrt(sum(error_norm_quat_MEKF.^2)/t_vec(end));

error_norm_quat_UKF = zeros(size(UKF_mu,1),1);
for i = 1:size(UKF_mu,1)
    error_norm_quat_UKF(i) = norm(quaternion(i,:)-UKF_mu(i,4:7));
end
RMSE_UKF = sqrt(sum(error_norm_quat_UKF.^2)/t_vec(end));

error_norm_omega_EKF = zeros(size(EKF_mu,1),1);
for i = 1:size(EKF_mu,1)
    error_norm_omega_EKF(i) = norm(omega(i,:)-EKF_mu(i,1:3));
end
RMSE_EKF_omega = sqrt(sum(error_norm_omega_EKF.^2)/t_vec(end));

error_norm_omega_MEKF = zeros(size(MEKF_mu,1),1);
for i = 1:size(MEKF_mu,1)
    error_norm_omega_MEKF(i) = norm(omega(i,:)-MEKF_mu(i,1:3));
end
RMSE_MEKF_omega = sqrt(sum(error_norm_omega_MEKF.^2)/t_vec(end));

error_norm_omega_UKF = zeros(size(UKF_mu,1),1);
for i = 1:size(UKF_mu,1)
    error_norm_omega_UKF(i) = norm(omega(i,:)-UKF_mu(i,1:3));
end
RMSE_UKF_omega = sqrt(sum(error_norm_omega_UKF.^2)/t_vec(end));

figure
hold on
sgtitle('Measurement Errors for \omega')
subplot(3,1,1)
hold on 
plot(t_vec,abs(omega(:,1) - UKF_mu(:,1)))
plot(t_vec,abs(omega(:,1) - MEKF_mu(:,1)))
plot(t_vec,abs(omega(:,1) - EKF_mu(:,1)))
ylabel('\Delta \omega_x (rad/s)')
legend('UKF Error','MEKF Error','EKF Error')

subplot(3,1,2)
hold on 
plot(t_vec,abs(omega(:,2) - UKF_mu(:,2)))
plot(t_vec,abs(omega(:,2) - MEKF_mu(:,2)))
plot(t_vec,abs(omega(:,2) - EKF_mu(:,2)))
ylabel('\Delta \omega_y (rad/s)')
legend('UKF Error','MEKF Error','EKF Error')

subplot(3,1,3)
hold on 
plot(t_vec,abs(omega(:,3) - UKF_mu(:,3)))
plot(t_vec,abs(omega(:,3) - MEKF_mu(:,3)))
plot(t_vec,abs(omega(:,3) - EKF_mu(:,3)))
ylabel('\Delta \omega_z (rad/s)')
legend('UKF Error','MEKF Error','EKF Error')
xlabel('Time (s)')

figure
hold on
sgtitle('Measurement Errors for Quaternion')
subplot(4,1,1)
hold on 
plot(t_vec,abs(quaternion(:,1) - EKF_mu(:,4)),'red')
plot(t_vec,abs(quaternion(:,1) - MEKF_mu(:,4)),'green')
plot(t_vec,abs(quaternion(:,1) - UKF_mu(:,4)),'blue')
ylabel('\Delta q_1')
legend('EKF Error','MEKF Error','UKF Error')

subplot(4,1,2)
hold on 
plot(t_vec,abs(quaternion(:,2) - EKF_mu(:,5)),'red')
plot(t_vec,abs(quaternion(:,2) - MEKF_mu(:,5)),'green')
plot(t_vec,abs(quaternion(:,2) - UKF_mu(:,5)),'blue')
ylabel('\Delta q_2')
legend('EKF Error','MEKF Error','UKF Error')


subplot(4,1,3)
hold on 
plot(t_vec,abs(quaternion(:,3) - EKF_mu(:,6)),'red')
plot(t_vec,abs(quaternion(:,3) - MEKF_mu(:,6)),'green')
plot(t_vec,abs(quaternion(:,3) - UKF_mu(:,6)),'blue')
ylabel('\Delta q_3')
legend('EKF Error','MEKF Error','UKF Error')

subplot(4,1,4)
hold on 
plot(t_vec,abs(quaternion(:,4) - EKF_mu(:,7)),'red')
plot(t_vec,abs(quaternion(:,4) - MEKF_mu(:,7)),'green')
plot(t_vec,abs(quaternion(:,4) - UKF_mu(:,7)),'blue')
ylabel('\Delta q_4')
legend('EKF Error','MEKF Error','UKF Error')
xlabel('Time (s)')



%% EKF Function 
function [mu_updated,sigma_updated] = EKF(mu,sigma,omega_estimated,quat_estimated,Q_process,R_sensor)

%set inertia tensor
Ix = 197.28;
Iy = 222.835;
Iz = 277.6;

%dt
dt = .001;

%H matrix
H = eye(7);

%here is measurement 
y = [omega_estimated';quat_estimated']';

%split up the EKF mu into the individual parts

%angular velocities 
wx = mu(1);
wy = mu(2);
wz = mu(3);
w = [wx wy wz];

%quaternions
q1 = mu(4);
q2 = mu(5);
q3 = mu(6);
q4 = mu(7);

%omega/omega part 
A(1:3,1:3) = [1, (dt/Ix)*(Iy-Iz)*w(3), (dt/Ix)*(Iy-Iz)*w(2);...
     (dt/Iy)*(Iz-Ix)*w(3),1,(dt/Iy)*(Iz-Ix)*w(1);...
     (dt/Iz)*(Ix-Iy)*w(2),(dt/Iz)*(Ix-Iy)*w(1),1];

%omega/quaternions part
A(1:3,4:7) = zeros(3,4);

%quaternion/omega part
A(4:7,1:3) = .5*dt*[ q4, -q3,  q2;...
               q3,  q4, -q1;...
              -q2,  q1,  q4;...
              -q1, -q2, -q3];

%quaternion/quaternion part            
A(4:7,4:7) = eye(4) + .5*dt*[0 wz -wy wx;...
                            -wz 0 wx wy;...
                             wy -wx 0 wz;...
                            -wx -wy -wz 0];

%predict              
mu_bar = A*mu';
P_bar = A*sigma*A' + Q_process;

%update
K = P_bar*H'*(inv(H*P_bar*H' + R_sensor));
mu_updated = mu_bar + K*(y' - mu');
sigma_updated = (eye(7) - K*H)*P_bar*(eye(7) - K*H)' + K*R_sensor*K';

%this sigma updat step is the joseph form, not the form we have used in 273
end

%% MEKF Function
function [mu_updated,sigma_updated] = MEKF(mu,sigma,omega_estimate,quat_estimated,Q_process,R_sensor)

%set inertia tensor
Ix = 197.28;
Iy = 222.835;
Iz = 277.6;

% sigmas
o = zeros(1,3);
o(1) = (Iy - Iz)/Ix;
o(2) = (Iz - Ix)/Iy;
o(3) = (Ix - Iy)/Iz;


% time step
dt = .001;

% withdraw previous
omega_previous = mu(1:3)';
quat_previous = mu(4:7)';

% no input for now
torque = zeros(3,1);

% generate predictions
% omega_predict = omega_previous + dt*(cross(diag([Ix,Iy,Iz])*omega_previous,...
%     omega_previous)+torque);
% 
r_prev = omega_previous/norm(omega_previous);
theta_prev = dt*norm(omega_previous);
quat_predict = qmult(quat_previous,[r_prev*sin(theta_prev/2);...
    cos(theta_prev/2)]);

% withdraw the rotation vector from the given quaternion
theta_estimated = 2*acos(quat_estimated(4));
r_estimated = quat_estimated(1:3)'/sin(theta_estimated/2);
theta_predict = 2*acos(quat_predict(4));
r_predict = quat_predict(1:3)/sin(theta_predict/2);

% generate A matrix
A = zeros(6,6);
% A(1:3,1:3) = eye(3) + dt*...
%     [0 o(1)*omega_previous(3) o(1)*omega_previous(2);...
%     o(2)*omega_previous(3) 0 o(2)*omega_previous(1);...
%     o(3)*omega_previous(2) o(3)*omega_previous(1) 0];
A(1:3,1:3) = [1, (dt/Ix)*(Iy-Iz)*omega_previous(3), (dt/Ix)*(Iy-Iz)*omega_previous(2);...
     (dt/Iy)*(Iz-Ix)*omega_previous(3),1,(dt/Iy)*(Iz-Ix)*omega_previous(1);...
     (dt/Iz)*(Ix-Iy)*omega_previous(2),(dt/Iz)*(Ix-Iy)*omega_previous(1),1];

% A(4:6,1:3) = eye(3)*dt;
A(4:6,4:6) = expm(-hat(omega_previous)*dt);

% find measurement error
omega_predict = A(1:3,1:3)*omega_previous;

direction_estimated = ROT_quat(quat_estimated);
direction_predicted = ROT_quat(quat_predict);
z = [omega_estimate' - omega_predict;...
    (direction_estimated(1,:)' - direction_predicted(1,:)');...
    (direction_estimated(2,:)' - direction_predicted(2,:)');...
    (direction_estimated(3,:)' - direction_predicted(3,:)')];

% generate C matrix
C = zeros(12,6);
C(1:3,1:3) = eye(3);
C(4:6,4:6) = hat(direction_estimated(1,:));
C(7:9,4:6) = hat(direction_estimated(2,:));
C(10:12,4:6) = hat(direction_estimated(3,:));

% predict variance
Sigma_predict = A*sigma(1:6,1:6)*A'+Q_process(1:6,1:6);

% find Kalman Gain
K = Sigma_predict * C' * inv(C * Sigma_predict * C' + R_sensor(1,1)*eye(12));

% linear update
dx_update = K * z;

% extrapolate new mu
theta_update = norm(dx_update(4:6));
r_update = dx_update(4:6)/theta_update;
omega_update = omega_previous + dx_update(1:3);
quat_update = qmult(quat_predict,...
    [r_update * sin(theta_update/2);cos(theta_update/2)]);

 
mu_updated = [omega_update;quat_update]';

% extrapolate new sigma
sigma_updated = (eye(6) - K * C)*Sigma_predict*(eye(6) - K * C)'...
    + K*R_sensor(1,1)*eye(12)*K';

end

function hat_matrix = hat(x)
% this function takes 3x1 in a vector and outputs it as a hat matrix

    hat_matrix = [0 x(3) -x(2); -x(3) 0 x(1); x(2) -x(1) 0];

end

function [product] = qmult(q1,q2)
% this funciton multiplies two quaternions together

v1 = q1(1:3);
v2 = q2(1:3);
s1 = q1(4);
s2 = q2(4);

product = [s2*eye(3)-hat(v2) v2; -v2' s2]*q1;

end

function [ROT] = ROT_quat(q)
% this function takes in a quaternion and pushes out a rotation matrix

v = q(1:3);
s = q(4);

ROT = eye(3) + 2*hat(v)*(s*eye(3)+hat(v));

end

function inv = q_inv(q)

inv = zeros(1,4);
inv(1:3) = -q(1:3);
inv(4) = q(4);

end

