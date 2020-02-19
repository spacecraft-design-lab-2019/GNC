
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

