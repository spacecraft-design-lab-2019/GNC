function [Z, Y, X, W, V] = functionA(A, B, C)
    %codegen
    
    % Testing MATLAB  coder for different matlab functionalities
    % Will try to run the tests using a .mex file
    %
    % * Creating a matrix inside the function
    % * Backlash (linear solvers)
    % * Matrix Slicing
    % * 3D matrices
    % * Matrix multiplication
    % * Matrix transposes
    % * Returning multiple values from a function
    
    % A = (4 x 4) matrix
    % B = (4 x 3 x 20) matrix
    % C = (4 x 1) vector
    
    % Test 0: Creating a matrix in the function
    Z = eye(4);
    Z(1, 2) = 3;
    Z(2, 4) = 5;
    
    % Test 1: Solving linear systems
    Y = A \ C;
    
    % Test 2: Slicing & 3D matrices
    X = B(:, :, 5);  % (4 x 3)
    
    % Test 3: Multiplication
    W = A * Z;
    
    % Test 3: Transpose
    V = C*C';  % Transpose
    

end