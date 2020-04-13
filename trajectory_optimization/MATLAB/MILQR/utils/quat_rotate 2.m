Copyright (c) 2020 Robotic Exploration Lab

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

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

