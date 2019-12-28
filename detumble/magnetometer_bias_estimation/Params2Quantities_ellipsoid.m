function [center, radii, evecs] = Params2Quantities_ellipsoid(v)
%This function takes in a vector of ellipsoid parameters:
% Ax^2 + By^2 + Cz^2 + 2Dxy + 2Exz + 2Fyz + 2Gx + 2Hy + 2Iz + J = 0
% and returns relevant quantities such as the radii, directions, and center

% This is the matrix that makes things linear in x (from De Leeuw)
W = [ v(1) v(4) v(5) v(7); ...
      v(4) v(2) v(6) v(8); ...
      v(5) v(6) v(3) v(9); ...
      v(7) v(8) v(9) v(10) ];
  
V_inv = W(1:3,1:3);

center = -V_inv\v( 7:9 );

% Translate to origin
% form the corresponding translation matrix
T = eye( 4 );
T( 4, 1:3 ) = center';
% translate to the center
R = T * A * T';

% solve the eigenproblem
[ evecs, evals ] = eig( R( 1:3, 1:3 ) / -R( 4, 4 ) );
radii = sqrt( 1 ./ diag( abs( evals ) ) );
sgns = sign( diag( evals ) );
radii = radii .* sgns;
end

