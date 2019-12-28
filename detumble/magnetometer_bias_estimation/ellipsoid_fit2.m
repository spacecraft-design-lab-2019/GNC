% New function for fitting ellipse based on stuff I understand
function v = ellipsoid_fit2( X )
x = X( :, 1 );
y = X( :, 2 );
z = X( :, 3 );

% Form A matrix (use Yuri Petrov's regularization, A+B+C = -3)
    D = [ x .* x + y .* y - 2 * z .* z, ...
        x .* x + z .* z - 2 * y .* y, ...
        2 * x .* y, ...
        2 * x .* z, ...
        2 * y .* z, ...
        2 * x, ...
        2 * y, ...
        2 * z, ...
        1 + 0 * x ];

% solve the normal system of equations
d2 = x .* x + y .* y + z .* z; % the RHS of the llsq problem (y's)
u = ( D' * D ) \ ( D' * d2 );

% (Un)regularize the coefficients
v = zeros(10,1);
v(1) = u(1) +     u(2) - 1;
v(2) = u(1) - 2 * u(2) - 1;
v(3) = u(2) - 2 * u(1) - 1;
v( 4 : 10 ) = u( 3 : 9 );


end