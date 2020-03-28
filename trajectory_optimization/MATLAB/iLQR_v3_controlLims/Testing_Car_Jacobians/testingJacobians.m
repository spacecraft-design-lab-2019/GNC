clear
clc


dX = 1e-6*eye(4);
dU = 1e-6*eye(2);
x0 = [1; 1; 3*pi/2; 0.1];
u0 = [0.1; 0.1];
dt = 0.03;

[xnew1 ,fx,fu] = car_step(x0,u0,dt);
for k = 1:4
  [xp(:,k), ~, ~] = car_step(x0+dX(:,k),u0,dt);
  [xm(:,k), ~, ~] = car_step(x0-dX(:,k),u0,dt);
end

for k = 1:2
    [xp2(:,k),~,~] = car_step(x0,u0+dU(:,k),dt);
    [xm2(:,k),~,~] = car_step(x0,u0-dU(:,k),dt);
end


% dynamics first derivatives
ix = 1:4;
iu = 5:6;
xu_dyn  = @(xu) car_dynamics(xu(ix,:),xu(iu,:));
J       = finite_difference(xu_dyn, [x0; u0]);
fx2      = J(:,ix,:);
fu2      = J(:,iu,:);
xnew = car_dynamics(x0, u0);

A = (xp - xm)/2e-6;
fx;
fx2;

B = (xp2-xm2)/2e-6;
fu;
fu2;


xnew1;
xnew;