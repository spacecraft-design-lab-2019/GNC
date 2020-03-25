function [x_b] = quat_rotate(x_a, q)
% Rotates a vector in R3 from frame {A} to frame {B} using
% a unit quaternion encoding the transformation
% 
% Inputs
%=====================
% x_a - the vector in frame {A} (3x1)
% q   - unit quaternion encoding rotation from frame {A} to frame {B}
%
% Outputs
%=====================
% x_b - the vector expressed in frame {B} (3x1)

assert(size(x_a,1)==3);
assert(size(x_a,2)==1);
assert(size(q,1)==4);
assert(size(q,2)==1);

L = [q(1),-q(2),-q(3),-q(4);
     q(2), q(1),-q(4), q(3);
     q(3), q(4), q(1),-q(2);
     q(4),-q(3), q(2), q(1)];
 
R = [q(1),-q(2),-q(3),-q(4);
     q(2), q(1), q(4),-q(3);
     q(3),-q(4), q(1), q(2);
     q(4), q(3),-q(2), q(1)];

H = [zeros(1,3);eye(3)];

% Perform the rotation
x_b = H'*L*R'*H*x_a;

end

