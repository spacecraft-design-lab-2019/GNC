#include <iostream>
#include "eigen-git-mirror/Eigen/Dense"


using Eigen::MatrixXd;

MatrixXd SkewSymmetric(double a1, double a2, double a3);
MatrixXd quat2dcm(MatrixXd q);

MatrixXd measurement(MatrixXd q, MatrixXd rN){
    MatrixXd R(3,3) = quat2dcm(q);
    MatrixXd rB1(3,1) = R*rN(0:3,0);
    MatrixXd rB2(3,1) = R*rN(3:6,0);
    MatrixXd y(6,1);
    y(0:3,0) = rB1;
    y(3:6,0) = rB2;
    MatrixXd C(6,6);
    C(0:3,0:3) = 2*SkewSymmetric(rB1(0),rB1(1),rB1(2));
    C(3:6,0:3) = 2*SkewSymmetric(rB2(0),rB2(1),rB2(2));
    return y,R,C;
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
    Q = SkewSymmetric(q2,q3,q4);
// calculate DCM
    MatrixXd DCM(3,3);
    DCM = (q1^2)-q_vec.transpose()*q_vec*identity(3)-2*q1*Q+2*q_vec*q_vec.transpose();
    return DCM
}

