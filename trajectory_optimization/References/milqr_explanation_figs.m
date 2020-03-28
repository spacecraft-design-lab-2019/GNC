clear;
close all;
clc;


% The true value function
f = @(x)(sin(2*x) + cos(x.^2));
x = linspace(0.5,3.5,300);
y = f(x);

% The quaddratic approximate value function at x_bar = 1.25
xbar = 1.43;
fnom = f(xbar);
fx = 2*cos(2*xbar) - 2*xbar*sin(xbar^2);
fxx = -2*sin(2*xbar) - 4*(xbar^2)*cos(xbar^2);
fapprox = @(dx)(fnom + fx*dx + fxx*(dx.^2));
dx = linspace(-0.25,1.2,200);
yapprox = fapprox(dx);

% Approx min
syms xd real
eqn = fnom + fx*xd + fxx*(xd^2);
[ymin, idx] = min(yapprox);
xmin = dx(idx) +xbar;

% True min
xt = 1.87;
yt = -1.5;



% Plot result
plot(x, y, 'LineWidth',1.5)
hold on
plot(xbar+dx, yapprox, '--', 'LineWidth',1.5)
plot(xmin, ymin, '-o','MarkerSize',5, 'MarkerEdgeColor', [0.8,0.2,0.3],'MarkerFaceColor',[0.8,0.2,0.3])
plot(xbar, fnom, '*','MarkerSize',7, 'MarkerEdgeColor',[0,0.5,1],'MarkerFaceColor',[0,0.5,1])
plot(xt, yt, 'o', 'MarkerSize', 6, 'MarkerEdgeColor',[0,0.5,1],'MarkerFaceColor',[0,0.5,1])
p = text(xbar-0.5, fnom, '$Q_k(\bar{x_k}, \bar{u}_k$)');
set(p, 'interpreter', 'Latex', 'FontSize',11);
xlabel('Control, u_k');
ylabel('Value, Q_k');
h = legend('True value function', 'Quadratic approximation about $\bar{u}_k$=1.4');
set(h,'interpreter','Latex','FontSize',11)
