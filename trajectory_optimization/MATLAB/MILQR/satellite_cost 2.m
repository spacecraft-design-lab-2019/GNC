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

function [cost, cx, cu, cxx, cuu] = satellite_cost(x, xg, u, terminal)
% Calculates the cost contribution of a given state and control
% Utilizes a standard quadratic cost function for the angular velocity
% and a geodesic cost for the attitude quaternion
% Also returns the 2nd order expansion of the cost function

% Inputs
%=====================================
% x        - [quaternion; omega] (7x1)
% u        - (3x1) (passed in as zeros for final time-step/terminal cost)
% terminal - integer 0 or 1
%---------------------------------------------------

R = 0.05*eye(3);   % control hessian

if terminal
    Qw = 1.0*eye(3);  % terminal angular velocity cost hessian
    w = 10;           % terminal geodesic cost weight
else
    Qw = 0.01*eye(3);  % cumulative angular velocity cost hessian
    w = 1;             % cumulative geodesic cost weighting
end

[quat_cost, sign] = calc_quat_cost(x, xg);
cost = w*quat_cost + (1/2)*(x(5:7)-xg(5:7))'*Qw*(x(5:7)-xg(5:7)) + (1/2)*u'*R*u;

% State cost Hessian
cxx = [-sign*w*eye(3)*(xg(1:4)'*x(1:4)), zeros(3,3);
       zeros(3,3), Qw];

% State cost Jacobian
Gq = [-x(2:4)'; x(1)*eye(3) + skew_mat(x(2:4))];
cx = [sign*w*Gq'*xg(1:4);
      Qw*(x(5:7)-xg(5:7))];

% Control cost Hessian & Jacobian
cuu = R;
cu = R*u;

end


function [quat_cost, sign] = calc_quat_cost(x, xg)
% Finds the geodesic quaternion-error cost
% quat_cost = min(1+qd'q, 1-qd'q)) where qd is the desired attitude quaternion 
% Also records the sign (+ or -) which minimizes quat_cost
% this is used when calculating the Jacobain and Hessian

quat_cost = xg(1:4)'*x(1:4);
if (1.0 + quat_cost) < (1.0 - quat_cost)
    quat_cost = (1.0 + quat_cost);
    sign = 1;
else
    quat_cost = (1.0 - quat_cost);
    sign = -1;
end
    
end

