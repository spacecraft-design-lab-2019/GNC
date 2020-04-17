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

function [B_eci] = fake_IGRF(r,t)
% This function takes in a position in ECI (in km) and outputs
% a B vector in teslas that's in the ECI frame

% INPUTS:
%   r - position in km
%   t - time in MJD

% first, translate the position into a lat and long
[lat,lon,alt] = convert_to_lat_lon_alt(r,t);

% find B_NED
B_ned = get_magnetic_field(lat,lon,alt,round(MJD2year(t)),5);

% convert to eci
B_eci = ned2eci(B_ned,lat*pi/180,lon*pi/180);

end

function year = MJD2year(t_MJD)
% gets the year out of MJD

year = t_MJD/365.25 + 1858 + 321/365;

end

function vec_eci = ned2eci(vec_ned,lat,lon)
% this function converts a vector from ned to eci

% first find ENU
E = [-sin(lon);cos(lon);0];
N = [-sin(lat)*cos(lon);-sin(lat)*sin(lon);cos(lat)];
U = [cos(lat)*cos(lon);cos(lat)*sin(lon);sin(lat)];

R = [E N U];

vec_eci = R*vec_ned;

end

function [lat,lon,alt] = convert_to_lat_lon_alt(r,t)
% this function converts between position in ECI and time to lat, lon and
% alt

% INPUTS:
%   r - position in km
%   t - time in MJD

% OUTPUTS:
%   lat - latitude in deg
%   lon - longitude in deg
%   alt - altitude in km

% convert to ecef
r_ecef = eci_to_ecef(r,t);

% get lat and lon
lon = atan2(r_ecef(2),r_ecef(1))*180/pi;
lat = asin(r_ecef(3)/norm(r_ecef))*180/pi;
alt = norm(r_ecef)-6378;

end

function r_ecef = eci_to_ecef(r,t)
% converts a position from eci to ecef

% find days since 1/1/2000 at 12h
theta = 280.4606 + 360.9856473*MJD2J2000(t);
R = Rz(theta*pi/180); % get the rotation matrix

r_ecef = R*r;

end

function t_J2000 = MJD2J2000(t_MJD)
% converts between MJD and J2000 (days since 1/1/2000 at 12h)
% MJD epoch: 0h Nov 17, 1858
% J2000 epoch: 12h Jan 1, 2000

t_J2000 = t_MJD  - 51544.500000;

end

function R = Rz(theta)
% creates a z axis rotation matrix (rad)

R = [cos(theta),sin(theta),0; ...
    -sin(theta),cos(theta),0; ...
    0,0,1];

end
