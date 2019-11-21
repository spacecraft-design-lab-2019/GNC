
clear all

%generate a sphere of points
[x,y,z] = sphere(20);
points = [];
for i = 1:length(x)^2
    points = [points;[x(i) y(i) z(i)]];
end

%stretch it
scale_factor = eye(3) + diag([.1 -.1 -.2]);
rotation = quaternion_rot([cosd(5)/2 sind(5)/2 0 0]);
mu = 0;
sigma = .05;
noise = normrnd(mu,sigma,size(points));

bias_values = normrnd(mu,sigma,[3,1]);
bias = zeros(size(points));
for i = 1:length(bias)
    bias(i,:) = bias_values;
end
meas = points*scale_factor*rotation+noise+bias;

%plot
figure;
scatter3(meas(:,1),meas(:,2),meas(:,3));
hold on
scatter3(points(:,1),points(:,2),points(:,3));
legend('Measured','Actual');

%now we have to solve it!
gains = (meas'*meas)^(-1)*meas'*points;
new = meas*gains;
figure
scatter3(new(:,1),new(:,2),new(:,3));
hold on
scatter3(points(:,1),points(:,2),points(:,3));
legend('Measured','Actual');

function ROT = quaternion_rot(q)
    % this function takes in a quaternion and spits out a rotation matrix
    w = q(1);
    x = q(2);
    y = q(3);
    z = q(4);
    ROT = [1-2*y^2-2*z^2 2*x*y-2*z*w 2*x*z+2*y*w;
       2*x*y+2*z*w 1-2*x^2-2*z^2 2*y*z-2*x*w;
       2*x*z-2*y*w 2*y*z+2*x*w 1-2*x^2-2*y^2];

end
