% Copyright (c) 2020 Robotic Exploration Lab
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.

function [X] = two_body_position(x0,t)
% This function takes in an initial position and velocity of a satellite along 
% with a time series and outputs the positions over the time using a simple
% point mass two body approximation

% INPUTS:
% x0 - initial state
%     r0 - initial position in km
%     v0 - initial velocity in km/s
%     q0 - initial quaternion
%     w0 - initial rotation rate
% t - time series in seconds

% OUTPUTS:
% X - array of states with columns:
%     r - position in Earth radii
%     v - velocity in km/s

days2sec = 24*60*60;

% initialize empty array
X = zeros(6,length(t));
X(:,1) = x0(1:6);

% fill it
for i = 1:length(t)-1
    
    X(:,i+1) = rk4(X(:,i),@two_body,(t(i+1)-t(i))*days2sec);
    
end

end

function x_new = rk4(x,func,dt)

k1 = dt*func(x);
k2 = dt*func(x+k1/2);
k3 = dt*func(x+k2/2);
k4 = dt*func(x+k3);

x_new = x + 1/6*(k1+2*k2+2*k3+k4);

end

function xdot = two_body(x)

Earth = InitializeEarth();

r = x(1:3);
v = x(4:6);

xdot = [v;-Earth.mu*r/(norm(r)^3)];

end
