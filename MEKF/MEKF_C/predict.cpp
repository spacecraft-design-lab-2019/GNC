#include <iostream>
#include "eigen-git-mirror/Eigen/Dense"


using Eigen::MatrixXd;

MatrixXd SkewSymmetric(double a1, double a2, double a3);

MatrixXd predict(MatrixXd xk, MatrixXd w, double dt){
    MatrixXd q(4,1);
    q(:,0) = xk(0:4);
    b(:,0) = xk(4:7);


    MatrixXd u(3,1);
    u = w + b;

    MatrixXd theta(3,1);
    theta = (sqrt(u(0,0)*b(0,0)+u(1,0)*b(1,0)+u(2,0)*b(2,0)))*dt;

    MatrixXd r(3,1);
    r = (u-b)/(sqrt(u(0,0)*b(0,0)+u(1,0)*b(1,0)+u(2,0)*b(2,0)));

    MatrixXd s(4,1); // q2
    s(0,0) = cos(theta/2);
    s(1,0) = r(0)*sin(theta/2);
    s(2,0) = r(1)*sin(theta/2);
    s(3,0) = r(2)*sin(theta/2);

    MatrixXd xn(7,1);
    xn(0,0) = q(0,0)*s(0,0)-q(1:4).transpose()*s(1:4);
    xn(1:4,0) = q(0,0)*s(1:4,0)+s(0,0)*q(1:4,0)+q(1:4,0).cross(s(1:4,0));

    MatrixXd V(3,4);
    V(0,0) = V(2,0) = V(3,0) = 0;
    V(1,0) = V(1,1) = V(1,3) = 0;
    V(2,0) = V(2,1) = V(2,2) = 0;
    V(0,1) = V(1,2) = V(2,3) = 1;

    // Identity Matrix - (probably a command for this?)
    MatrixXd I(3,3);
    I(0,0) = I(1,1) = I(2,2) = 1;
    I(0,1) = I(0,2) = I(1,2) = I(2,1) = I(1,0) = I(2,0) = 0;

    // left quaternion multiply
    MatrixXd Ls(4,4);
    Ls(0,0) = s(0,0);
    Ls(0,1:4) = -s(1:4,0).transpose();
    Ls(1:4,0) = s(1:4,0);
    Ls(1:4,1:4) = s(0,0) * I + SkewSymmetric(s(1,0),s(2,0),s(3,0));

    // right quaternion multiply
    MatrixXd Rs(4,4);
    Rs(0,0) = s(0,0);
    Rs(0,1:4) = -s(1:4,0).transpose();
    Rs(1:4,0) = s(1:4,0);
    Rs(1:4,1:4) = s(0,0) * I - SkewSymmetric(s(1,0),s(2,0),s(3,0));

    MatrixXd A(6,6);
    A(0:3,0:3) = V*Ls.tranpose()*Rs*V.transpose();
    A(0:3,3:6) = -0.5*dt*I;
    A(3:6,3:6) = I;

    return xn,A
}


MatrixXd SkewSymmetric(double a1, double a2, double a3){
    MatrixXd a(3,3);
    a(1,1) = 0;
    a(2,2) = 0;
    a(3,3) = 0;
    a(1,2) = -a3;
    a(1,3) = a2;
    a(2,1) = a3;
    a(2,3) = -a1;
    a(3,1) = -a2;
    a(3,2) = a1;
    return a;
}


