function [P_inv_updated, q_updated, bias_updated] = UpdateEllipse(P_inv, q, X_i)
%This function updates the outer product matrix, P, RHS vector, q, using
%the data found in X_i
x = X_i(1);
y = X_i(2);
z = X_i(3);

% Form a vector (use Yuri Petrov's regularization, A+B+C = -3)
a = [ x .* x + y .* y - 2 * z .* z, ...
    x .* x + z .* z - 2 * y .* y, ...
    2 * x .* y, ...
    2 * x .* z, ...
    2 * y .* z, ...
    2 * x, ...
    2 * y, ...
    2 * z, ...
    1 + 0 * x ]';
y = x^2 + y^2 + z^2;

% Rank one update formula, from EE 263 slides
P_inv_updated = P_inv - 1/(1+a'*P_inv*a)*(P_inv*a)*(P_inv*a)';
q_updated = q + y*a;

u = P_inv_updated*q_updated;

% This is the matrix that makes things linear in x (from De Leeuw)
v = zeros(10,1);
v(1) = u(1) +     u(2) - 1;
v(2) = u(1) - 2 * u(2) - 1;
v(3) = u(2) - 2 * u(1) - 1;
v( 4 : 10 ) = u( 3 : 9 );

W = [ v(1) v(4) v(5) v(7); ...
      v(4) v(2) v(6) v(8); ...
      v(5) v(6) v(3) v(9); ...
      v(7) v(8) v(9) v(10) ];
  
V_inv = W(1:3,1:3);

bias_updated = -V_inv\v( 7:9 );

end

