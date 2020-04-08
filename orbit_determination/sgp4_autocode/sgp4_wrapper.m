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

% This function is a wrapper to interact with the sgp4 model for autocoding

% create a new sgp4
% INPUTS:
% orbit -
%     satn        - satellite number
%     bstar       - sgp4 type drag coefficient              kg/m2er
%     ecco        - eccentricity
%     epoch       - epoch time in days from jan 0, 1950. 0 hr
%     argpo       - argument of perigee (output if ds)
%     inclo       - inclination
%     mo          - mean anomaly (output if ds)
%     no          - mean motion
%     nodeo      - right ascension of ascending node

% OUTPUTS:
% sgp4 - 

function [sgp4] = create_sgp4(orbit)

% parse out
satn = orbit(1);
bstar = orbit(2);
ndot = orbit(3);
nddot = orbit(4);
ecco = orbit(5);
epoch = orbit(6); % days at Oct 12, 2020
argpo = orbit(7); % rad
inclo = orbit(8);  %rad
mo = orbit(9); % rad
no = orbit(10); %rad/s
nodeo = orbit(11); % rad

% run sgp4init
sgp4 = [];
sgp4 = sgp4init(84,1,sgp4,epoch,bstar, ndot, nddot, ...
         ecco, argpo, inclo, mo, no, nodeo);

end