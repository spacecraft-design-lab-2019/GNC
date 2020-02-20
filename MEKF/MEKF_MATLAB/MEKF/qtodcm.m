function R = qtodcm(q)

    R = eye(3) + 2*hat(q(1:3))*(hat(q(1:3)) + eye(3)*q(4));
    
end