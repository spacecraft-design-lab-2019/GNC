function [x_skew] = skew_mat(x)
% Returns skew symmetric - cross product matrix of a 3x1 vector

x_skew = [0 -x(3) x(2); x(3) 0 -x(1); -x(2) x(1) 0];

end