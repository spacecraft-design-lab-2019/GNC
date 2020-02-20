function x = vee(xhat)
%VEE Computes the vee map (inverse hat map), mapping a 3X3 skew-symmetric matrix
%into a 3-vector

x = zeros(3,1);

if xhat(3,2) == -xhat(2,3)
    x(1) = xhat(3,2);
else
    x(1) = (xhat(3,2) - xhat(2,3))/2;
end

if xhat(1,3) == -xhat(3,1)
    x(2) = xhat(1,3);
else
    x(2) = (xhat(1,3) - xhat(3,1))/2;
end

if xhat(2,1) == -xhat(1,2)
    x(3) = xhat(2,1);
else
    x(3) = (xhat(2,1) - xhat(1,2))/2;
end

end

