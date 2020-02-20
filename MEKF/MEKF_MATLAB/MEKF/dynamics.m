function xdot = dynamics(t,x,J,tau)

q = x(1:4);
q = q/sqrt(q'*q);
w = x(5:7);

v = q(1:3);
s = q(4);

qdot = 0.5*[s*w + cross(v,w); -v'*w];

wdot = J\(tau - cross(w,J*w));

xdot = [qdot; wdot];

end

