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

function [B_eci_vec] = get_magnetic_field_series(x0,t)
% This function takes in a series of positions in ECI (in km) as well as a 
% time series in MJD and outputs a series of ECI magnetic field vectors

% INPUT:
% x0 - [r0;v0;q0;w0];
% t - time vector in MJD

% OUTPUT:
% B_eci_vec - magnetic field in ECI in Teslas

X = two_body_position(x0,t);

B_eci_vec = zeros(3,length(t)-1);

for i = 1:length(t)-1
    B_eci_vec(:,i) = fake_IGRF(X(1:3,i),t(i))*1E-9;
end


end
