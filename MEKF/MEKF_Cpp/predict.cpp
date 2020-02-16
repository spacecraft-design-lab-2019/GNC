#include <iostream>
#include "MEKF.hpp"
#include <cmath>
#include <stdio.h>

using namespace std;
using Eigen::MatrixXd;

MatrixXd SkewSymmetric(double a1, double a2, double a3);

void predict(MatrixXd xk, MatrixXd w, double dt, MatrixXd &xn, MatrixXd &A){
    MatrixXd q(4,1), b(3,1);
    q(all,0) = xk(seq(0,3),0);
    b(all,0) = xk(seq(4,6),0);


    MatrixXd u(3,1);
    u = w + b;

    double theta;
    theta = dt*sqrt(u(0,0)*b(0,0)+u(1,0)*b(1,0)+u(2,0)*b(2,0));

    MatrixXd r(3,1);
    r = (u-b);
    r(0,0) = r(0,0)/(sqrt(u(0,0)*b(0,0)+u(1,0)*b(1,0)+u(2,0)*b(2,0)));
    r(1,0) = r(1,0)/(sqrt(u(0,0)*b(0,0)+u(1,0)*b(1,0)+u(2,0)*b(2,0)));
    r(2,0) = r(2,0)/(sqrt(u(0,0)*b(0,0)+u(1,0)*b(1,0)+u(2,0)*b(2,0)));
    MatrixXd s(4,1); // q2
    s(0,0) = cos(theta/2);
    s(1,0) = r(0)*sin(theta/2);
    s(2,0) = r(1)*sin(theta/2);
    s(3,0) = r(2)*sin(theta/2);
    xn(0,0) = q(0,0)*s(0,0)-(q(1,0)*s(1,0)+q(2,0)*s(2,0)+q(3,0)*s(3,0));


    // cross product of s and q
    MatrixXd cross_sq(3,1);
    cross_sq(0,0) = q(2,0)*s(3,0)-s(2,0)*q(3,0);
    cross_sq(1,0) = q(3,0)*s(1,0)-s(3,0)*q(1,0);
    cross_sq(2,0) = q(1,0)*s(2,0)-s(1,0)*q(2,0);

    xn(seq(1,3),0) = q(0,0)*s(seq(1,3),0)+s(0,0)*q(seq(1,3),0)+cross_sq;
    xn(seq(4,6),0) = b;



    MatrixXd V(3,4);
    V(0,1) = 1.0;
    V(1,2) = 1.0;
    V(2,3) = 1.0;
    // left quaternion multiply
    MatrixXd Ls(4,4);
    Ls(0,0) = s(0,0);
    Ls(0,seq(1,3)) = -s(seq(1,3),0).transpose();
    Ls(seq(1,3),0) = s(seq(1,3),0);
    Ls(seq(1,3),seq(1,3)) = s(0,0) * MatrixXd::Identity(3,3) + SkewSymmetric(s(1,0),s(2,0),s(3,0));
    // right quaternion multiply
    MatrixXd Rs(4,4);
    Rs(0,0) = s(0,0);
    Rs(0,seq(1,3)) = -s(seq(1,3),0).transpose();
    Rs(seq(1,3),0) = s(seq(1,3),0);
    Rs(seq(1,3),seq(1,3)) = s(0,0) * MatrixXd::Identity(3,3) - SkewSymmetric(s(1,0),s(2,0),s(3,0));

    A(seq(0,2),seq(0,2)) = V*Ls.transpose()*Rs*V.transpose();
    A(seq(0,2),seq(3,5)) = -0.5*dt*MatrixXd::Identity(3,3);
    A(seq(3,5),seq(3,5)) = MatrixXd::Identity(3,3);
}

MatrixXd SkewSymmetric(double a1, double a2, double a3){
    MatrixXd a(3,3);
    a(0,0) = 0.0;
    a(1,1) = 0.0;
    a(2,2) = 0.0;
    a(0,1) = -a3;
    a(0,2) = a2;
    a(1,0) = a3;
    a(1,2) = -a1;
    a(2,0) = -a2;
    a(2,1) = a1;
    return a;
}


