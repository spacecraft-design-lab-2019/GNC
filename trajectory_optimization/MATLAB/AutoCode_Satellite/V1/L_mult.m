function [L] = L_mult(q) %#codegen
% Forms the "left matrix" for quaternion multiplication
% where q*q1 = L(q)q1

% Input
% q - the quaternion to build the left matrix from

% Output
% L - the left matrix
% L = [s, -v';
%      v, sI+skew(v)]
%---------------------------------------------------

assert(size(q,1)==4);

L = [q(1),-q(2),-q(3),-q(4);
     q(2), q(1),-q(4), q(3);
     q(3), q(4), q(1),-q(2);
     q(4),-q(3), q(2), q(1)];
end