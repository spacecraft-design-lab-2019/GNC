function [cost, cx, cu, cxx, cuu] = car_cost2(x, xg, u, terminal)

cost = cost_func(x,u,terminal);
[cx,cu,cxx,cuu] = cost_derivs(x,u,terminal);
xg;

end


function c = cost_func(x, u, terminal)
% cost function for car-parking problem
% sum of 3 terms:
% lu: quadratic cost on controls
% lf: final cost on distance from target parking configuration
% lx: running cost on distance from origin to encourage tight turns

% final = isnan(u(1,:));
% u(:,final)  = 0;

cu  = 1e-2*[1 .01];         % control cost coefficients

cf  = [ .1  .1   1  .3];    % final cost coefficients
pf  = [.01 .01 .01  1]';    % smoothness scales for final cost

cx  = 1e-3*[1  1];          % running cost coefficients
px  = [.1 .1]';             % smoothness scales for running cost

% control cost
lu = cu*u.^2;

% final cost
if terminal
   lf = cf*sabs(x,pf);
else
   lf = 0;
end

% running cost
lx = cx*sabs(x(1:2,:),px);

% total cost
c = lu + lx + lf;

end

function y = sabs(x,p)
% smooth absolute-value function (a.k.a pseudo-Huber)
y = pp( sqrt(pp(x.^2,p.^2)), -p);
end


function [cx,cu,cxx,cuu] = cost_derivs(x,u, terminal)
% state and control indices
ix = 1:4;
iu = 5:6;

% cost first derivatives
xu_cost = @(xu, terminal) cost_func(xu(ix,:),xu(iu,:), terminal);
J       = squeeze(finite_difference(xu_cost, [x; u], terminal));
cx      = J(ix);
cu      = J(iu);

% cost second derivatives
xu_Jcst = @(xu, terminal) squeeze(finite_difference(xu_cost, xu, terminal));
JJ      = finite_difference(xu_Jcst, [x; u], terminal);
JJ      = 0.5*(JJ + permute(JJ,[2 1 3])); %symmetrize
cxx     = JJ(ix,ix);
cuu     = JJ(iu,iu);

end

function J = finite_difference(fun, x, terminal)
% simple finite-difference derivatives
% assumes the function fun() is vectorized

h = 2^-17;

[n, K]  = size(x);
H       = [zeros(n,1) h*eye(n)];
H       = permute(H, [1 3 2]);
X       = pp(x, H);
X       = reshape(X, n, K*(n+1));
Y       = fun(X, terminal);
m       = numel(Y)/(K*(n+1));
Y       = reshape(Y, m, K, n+1);
J       = pp(Y(:,:,2:end), -Y(:,:,1)) / h;
J       = permute(J, [1 3 2]);

end

function c = pp(a,b)
c = bsxfun(@plus,a,b);
end

function c = tt(a,b)
c = bsxfun(@times,a,b);
end
