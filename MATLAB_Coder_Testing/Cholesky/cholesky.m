function [R] = cholesky(A)
% Returns the Cholesky decomposition of a Hermitian
% positive-definite matrix A

R = chol(A, 'lower');
end

