#include <iostream>
#include "MEKF.hpp"
#include <cmath>
#include <stdio.h>

using namespace std;
using Eigen::MatrixXd;

MatrixXd SkewSymmetric2(double a1, double a2, double a3);
MatrixXd quat2dcm(MatrixXd q);

void measurement(MatrixXd q, MatrixXd rN, MatrixXd &y, MatrixXd &R, MatrixXd &C){
    R = quat2dcm(q);
    MatrixXd rB1(3,1);
    rB1 = R*rN(seq(0,2),0);
    MatrixXd rB2(3,1);
    rB2 = R*rN(seq(3,5),0);

    y(seq(0,2),0) = rB1;
    y(seq(3,5),0) = rB2;

    C(seq(0,2),seq(0,2)) = 2*SkewSymmetric2(rB1(0),rB1(1),rB1(2));
    C(seq(3,5),seq(0,2)) = 2*SkewSymmetric2(rB2(0),rB2(1),rB2(2));
    C(all,seq(3,5)) = MatrixXd::Zero(6,3);
}


MatrixXd quat2dcm(MatrixXd q){
    double q1 = q(0,0); // scalar
    double q2 = q(1,0); // start vec
    double q3 = q(2,0);
    double q4 = q(3,0);
    MatrixXd q_vec(3,1);
    q_vec(0,0) = q2;
    q_vec(1,0) = q3;
    q_vec(2,0) = q4;
// calculate skew symmetric matrix Q
    MatrixXd Q(3,3);
    Q = SkewSymmetric2(q2,q3,q4);
// calculate DCM
    MatrixXd DCM(3,3);
    DCM = (pow(q1,2.0)-q2*q2-q3*q3-q4*q4)*MatrixXd::Identity(3,3)-2.0*q1*Q+2.0*q_vec*q_vec.transpose();
    return DCM;
}
MatrixXd SkewSymmetric2(double a1, double a2, double a3){
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

