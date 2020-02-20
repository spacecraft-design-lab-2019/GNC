function [xn,A] = prediction(xk,w,dt)

q = xk(1:4);
b = xk(5:7);

p = [0.5*dt*(w-b); 1-0.5*((0.5*dt*(w-b))'*(0.5*dt*(w-b)))];
p = p/sqrt(p'*p);

qn = [q(4)*p(1:3) + p(4)*q(1:3) + cross(q(1:3),p(1:3)); q(4)*p(4) - q(1:3)'*p(1:3)];

% theta = norm(w-b)*dt;
% r = dt*(w-b)/theta;
% p2 = [r*sin(theta/2); cos(theta/2)];
% qn2 = qmult(q,p2);

xn = [qn; b];
A = [qtodcm(p)', -dt*eye(3); zeros(3), eye(3)];

end

