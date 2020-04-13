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

