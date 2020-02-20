function l = logq(q)
%LOGQ computes the logarithm of a quaternion

v = q(1:3);
nv = sqrt(v'*v);
nq = sqrt(q'*q);

l = [(v/nv)*acos(q(4)/nq); log(nq)];

end

