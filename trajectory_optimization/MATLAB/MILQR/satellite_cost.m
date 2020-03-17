function [cost, cx, cu, cxx, cuu] = satellite_cost(x, xg, u, terminal)
% Calculates the cost contribution of a given state and control 
% Also calculates the 2nd order expansion of the cost function

% Inputs
%=====================================
% x        - [quaternion; omega] (7x1)
% u        - (3x1) (passed in as zeros for final time-step/terminal cost)
% terminal - int 0 or 1
%---------------------------------------------------

R = 0.05*eye(3);   % control hessian

if terminal
    Qw = 1*eye(3);  % terminal angular velocity hessian 
else
    Qw = 0.01*eye(3);  % angular velocity hessian
end

[quat_cost, sign] = calc_quat_cost(x, xg);
cost = quat_cost + (1/2)*(x(5:7)-xg(5:7))'*Qw*(x(5:7)-xg(5:7)) + (1/2)*u'*R*u;

% State cost Hessian
cxx = [sign*eye(3)*(xg(1:4)'*x(1:4)), zeros(3,3);
    zeros(3,3), Qw];

% State cost Jacobian
Gq = [-x(2:4)'; x(1)*eye(3) + skew_mat(x(2:4))];
cx = [sign*Gq'*xg(1:4);
      Qw*(x(5:7)-xg(5:7))];

% Control cost Hessian & Jacobian
cuu = R;
cu = R*u;

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

