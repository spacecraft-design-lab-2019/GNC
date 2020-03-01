function y = car_dynamics(x,u)

% === states and controls:
% x = [x y t v]' = [x; y; car_angle; front_wheel_velocity]
% u = [w a]'     = [front_wheel_angle; acceleration]

% constants
d  = 2.0;      % d = distance between back and front axles
h  = 0.03;     % h = timestep (seconds)

% controls
w  = u(1,:,:); % w = front wheel angle
a  = u(2,:,:); % a = front wheel acceleration

o  = x(3,:,:); % o = car angle
               % z = unit_vector(o)
z  = [cos(o); sin(o)]; 

v  = x(4,:,:); % v = front wheel velocity
f  = h*v;      % f = front wheel rolling distance
               % b = back wheel rolling distance
b  = d + f.*cos(w) - sqrt(d^2 - (f.*sin(w)).^2);
               % do = change in car angle
do = asin(sin(w).*f/d);

dy = [tt(b, z); do; h*a];   % change in state
y  = x + dy;                % new state


function c = tt(a,b)
c = bsxfun(@times,a,b);