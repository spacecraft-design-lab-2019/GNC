function [q,R] = triad(rN1,rN2,rB1,rB2)

N = [rN1, rN2, cross(rN1,rN2)];
B = [rB1, rB2, cross(rB1,rB2)];

R = N*B^-1;

q = dcmtoq(R);

end

