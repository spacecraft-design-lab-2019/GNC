% Script for testing functionA() Matlab vs C
clear 
close all 
clc

% Testing MATLAB  coder for different matlab functionalities
% Will try to run the tests using both MATLAB and C (via a .mex file)
%
% What functionA tests:
%----------------------
% * Creating a matrix inside the function
% * Backlash (linear solvers)
% * Matrix Slicing
% * 3D matrices
% * Matrix multiplication
% * Matrix transposes
% * Returning multiple values from a function

% Define function inputs
% A (4x4)
A = eye(4);
A(2,2) =  3;
A(3,2) = 2;
A(1,4) = 5;
A(3,1) = 4;

% B (4x3x20)
B = zeros(4, 3, 20);
B(:, :, 5) = [9, 2, 1;
              0, 1, 1;
              4, 5, 6;
              1, 4, 0];
% C (4x1)
C = [5; 4; 3; 2];


% Calling the MATLAB version of functionA
[Z, Y, X, W, V] = functionA(A, B, C);

% Call the mex file (which uses the C version of functionA)
% [Zmex, Ymex, Xmex, Wmex, Vmex] = functionA_mex(A, B, C);


% Comparison (should all be the same)
% display(Z)
% display(Zmex);
% display(Y)
% display(Ymex);
% display(X)
% display(Xmex);
% display(W)
% display(Wmex);
% display(V)
% display(Vmex);







