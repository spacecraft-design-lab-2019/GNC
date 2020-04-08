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

clear; clc;

addpath('sgp4')

deg2rad = pi/180;
rad2deg = 180/pi;

%     satn        - satellite number
%     bstar       - sgp4 type drag coefficient              kg/m2er
%     ecco        - eccentricity
%     epoch       - epoch time in days from jan 0, 1950. 0 hr
%     argpo       - argument of perigee (output if ds)
%     inclo       - inclination
%     mo          - mean anomaly (output if ds)
%     no          - mean motion
%     nodeo      - right ascension of ascending node

satn = 1;
bstar = -.11606E-4;
ecco = 0;
epoch = 25852; % days at Oct 12, 2020
ndot = 0;
nddot = 0;
argpo = 0; % deg
inclo = 96;  %deg
mo = 0; % deg
no = 0.001106816457 * rad2deg; %deg/s
nodeo = 0; % deg

fake_input = [satn,bstar,ecco,epoch,ndot,nddot,argpo,inclo,mo,no,nodeo];

sgp4_struct = sgp4_wrapper(fake_input);

% test to see if it works
[sgp4, r, v] = sgp4(sgp4_struct,10);