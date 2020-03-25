function [X] = chol_solve(L,B)
% Solves the linear system AX = B for the unknown X using the Cholesky
% decomposition L of the matrix A. Where LL' = A
% X can be a vector or a matrix of size n x m

% Solution is found in O(nm) time using back-substitution
% This implementation only works for lower triangular factorisations

n = size(L,1);
m = size(B,2);

% Check sizes match and L is lower-triangular
assert(size(L,2)==n,'L must be square');
assert(size(B,1)==n,'The first dimension of b must match the shape of L');
assert(istril(L),'L must be lower triangular');


% Solve LY = B for Y
%=======================
Y = zeros(n,m);
coeffs = zeros(n,1);
for i = 1:n
    if i ~= 1
        coeffs(1:i-1) = L(i,1:i-1);
    end
    
    for j = 1:m
        Y(i,j) = (B(i,j) - coeffs'*Y(:,j))/L(i,i);
    end
end

% Solve L'X = Y for X
%=======================
X = zeros(n,m);
LT = L';
coeffs = zeros(n,1);
for i = n:-1:1
    if i ~= n
        coeffs(i+1:n) = LT(i,i+1:n);
    end
    
    for j = 1:m
        X(i,j) = (Y(i,j) - coeffs'*X(:,j))/LT(i,i);
    end
end
end

