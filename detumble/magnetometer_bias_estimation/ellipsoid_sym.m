% Script for analyzing ellipsoid stuff
syms x y z c_x c_y c_z P11 P12 P13 P22 P23 P33 real

vec = [x-c_x;y-c_y;z-c_z];
P = [P11, P12, P13; P12, P22, P23; P13, P23, P33];

ellipsoid = vec'*P*vec