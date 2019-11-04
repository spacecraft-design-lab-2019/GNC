#include <iostream>
#include "eigen-git-mirror/Eigen/Dense"

using Eigen::MatrixXd;

MatrixXd quatmult(MatrixXd q1, MatrixXd q2);
MatrixXd cross(MatrixXd a, MatrixXd b);

MatrixXd quatmult(MatrixXd q1, MatrixXd q2){
    MatrixXd qn_temp(4,1);
    qn_temp(0,0) = q1(0,0)*q2(0,0)-q1(1:4,0).transpose()*q2(1:4,0);
    qn_temp(1:4,0) = q1(0,0)*q2(1:4,0)+q2(0,0)*q1(1:4,0)+cross(q1(1:4,0),q2(1:4,0));
    MatrixXd qn(4,1);
    qn = qn_temp / sqrt(qn_temp(0,0)^2+qn_temp(1,0)^2+qn_temp(2,0)^2+qn_temp(3,0)^2);
    return qn;
}

MatrixXd cross(MatrixXd a, MatrixXd b){
    MatrixXd c(3,1);
    c(0,0) = a(1,0)*b(2,0)-a(2,0)*b(1,0);
    c(1,0) = -a(0,0)*b(2,0)+a(2,0)*b(0,0);
    c(2,0) = a(0,0)*b(1,0)-a(1,0)*b(0,0);
    return c;
}

MatrixXd update(MatrixXd L, MatrixXd z, MatrixXd xn, MatrixXd Pn, MatrixXd V, MatrixXd C){
    MatrixXd dx(6,1) = L * z;
    MatrixXd dphi(3,1) = dx(0:2,0);
    MatrixXd dbeta(3,1) = dx(3:5,0);
    double temp = sqrt(1-dphi.transpose()*dphi);
    MatrixXd temp_q(4,1);
    temp_q(0,0) = temp;
    temp_q(1:4,0) = dphi;
    MatrixXd q(4,1);
    q = quatmult(xn(0:4,0),temp_q);
    MatrixXd b(3,1);
    b = xn(4:7,0) + dbeta;
    MatrixXd x(7,1);
    x(0:4,0) = q;
    x(4:,0) = b;
    MatrixXd P(6,6);
    P = (identity(6)-L*C)*Pn*(identity(6)-L*C).transpose() + L*V*L.transpose();
    return x, P;
}
