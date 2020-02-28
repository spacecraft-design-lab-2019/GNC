function [cost, cx, cu, cxx, cuu] = car_cost(x, xg, u, terminal)
% Calculates the cost component from the input state/control
% Also calculates the cost derivates

% Cost matrices
% Cumulative
Q = 0.1*eye(4);
R = 0.01*eye(2);

% Terminal
Qf = eye(4);
Qf(3,3) = 100;

if terminal
    % Terminal cost
    cost = (1/2)*(x-xg)'*Qf*(x-xg);

    % Final cost gradients
    cxx = Qf;
    cx = Qf*(x-xg);
else
    % Cumulative cost
    cost = (1/2)*(x-xg)'*Q*(x-xg) + (1/2)*u'*R*u;
    
    % Cost gradients
    cxx = Q;
    cx = Q*(x-xg);
end
% Control cost gradients
cuu = R;
cu = R*u;

end