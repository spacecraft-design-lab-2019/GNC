function [cost, cx, cu, cxx, cuu] = satellite_cost(x, xg, u, terminal)
% Calculates the cost component from the input state/control
% Also calculates the cost derivates

% Cost matrices
% Cumulative
Q = zeros(7, 7);
Q(5:7, 5:7) = 0.1*eye(3);
R = 0.5*eye(3);

% Terminal
Qf= zeros(7, 7);
Qf(5:7, 5:7) = 10*eye(3);

if terminal
    % Terminal cost
    [quat_cost, sign] = calc_quat_cost(x, xg);
    cost = quat_cost + (1/2)*(x-xg)'*Qf*(x-xg);

    % Final cost gradients
    cxx = Qf;
    cx = Qf*(x-xg);
    cx(1:4) = sign*xg(1:4);
    cuu = R;  % unused
    cu = R*u; % unused
      
else
    % Cumulative cost
    [quat_cost, sign] = calc_quat_cost(x, xg);
    cost = quat_cost + (1/2)*(x-xg)'*Q*(x-xg) + (1/2)*u'*R*u;
    
    % Cost gradients
    cxx = Q;
    cx = Q*(x-xg);
    cx(1:4) = sign*xg(1:4);
    cuu = R;
    cu = R*u;
end 

end


function [quat_cost, sign] = calc_quat_cost(x, xg)
    % Find the linear quaternion-error cost (qcost = min(1+qg'q, 1-qg'q))
    % Also record the sign for use in the backward pass
    
    quat_cost = xg(1:4)'*x(1:4);
    if (1.0 + quat_cost) < (1.0 - quat_cost)
        quat_cost = (1.0 + quat_cost);
        sign = 1;
    else
        quat_cost = (1.0 - quat_cost);
        sign = -1;
    end
    
end

