function [xdot, dxdot] = satellite_dynamics(t,x,u)

J = 0.01*eye(3); % kgm^2

% Angular velocity
w = x(5:7);

% Quaternion components
s = x(1);
v = x(2:4);

% Non-linear dynamics
qdot = 0.5*[-v'; s*eye(3) + skew_mat(v)]*w;
wdot = -inv(J)*(skew_mat(w)*J*w - u);
xdot = [qdot; wdot];

% Jacobians
A = [0, -w', -v';
    w, -skew_mat(w), s*eye(3)+skew_mat(v);
    zeros(3,4), -inv(J)*(skew_mat(w)*J - skew_mat(J*w))];

B = [zeros(4,3);
    inv(J)];

dxdot = [A, B];



function [x_skew] = skew_mat(x)
% Returns skew symmetric - cross porduct matrix of a vector

x_skew = [0 -x(3) x(2); x(3) 0 -x(1); -x(2) x(1) 0];

end

end
