#include <iostream>
#include "MEKF.hpp"
#include <cmath>
#include <stdio.h>

using namespace std;
using Eigen::MatrixXd;

MatrixXd quatmult(MatrixXd q1, MatrixXd q2);
MatrixXd cross(MatrixXd a, MatrixXd b);

void update(MatrixXd L, MatrixXd z, MatrixXd xn, MatrixXd Pn, MatrixXd V, MatrixXd C, MatrixXd &x, MatrixXd &P){
    MatrixXd dx(6,1);
    dx = L * z;
    MatrixXd dphi(3,1);
    dphi = dx(seq(0,2),0);
    MatrixXd dbeta(3,1);
    dbeta = dx(seq(3,5),0);
    double temp;
    temp = sqrt(1.0-pow(dphi(0,0),2)+pow(dphi(1,0),2)+pow(dphi(2,0),2));
    MatrixXd temp_q(4,1);
    temp_q(0,0) = temp;
    temp_q(seq(1,3),0) = dphi;
    MatrixXd q(4,1);
    q = quatmult(xn(seq(0,3),0),temp_q);
    MatrixXd b(3,1);
    b = xn(seq(4,6),0) + dbeta;
    x(seq(0,3),0) = q;
    x(seq(4,6),0) = b;
    P = (MatrixXd::Identity(6,6)-L*C)*Pn*(MatrixXd::Identity(6,6)-L*C).transpose() + L*V*L.transpose();
}

MatrixXd quatmult(MatrixXd q1, MatrixXd q2){
    MatrixXd qn_temp(4,1);
    qn_temp(0,0) = q1(0,0)*q2(0,0)-q1(seq(1,3),0).transpose()*q2(seq(1,3),0);
    qn_temp(seq(1,3),0) = q1(0,0)*q2(seq(1,3),0)+q2(0,0)*q1(seq(1,3),0)+cross(q1(seq(1,3),0),q2(seq(1,3),0));
    MatrixXd qn(4,1);
    qn = qn_temp / sqrt(pow(qn_temp(0,0),2)+pow(qn_temp(1,0),2)+pow(qn_temp(2,0),2)+pow(qn_temp(3,0),2));
    return qn;
}

MatrixXd cross(MatrixXd a, MatrixXd b){
    MatrixXd c(3,1);
    c(0,0) = a(1,0)*b(2,0)-a(2,0)*b(1,0);
    c(1,0) = -a(0,0)*b(2,0)+a(2,0)*b(0,0);
    c(2,0) = a(0,0)*b(1,0)-a(1,0)*b(0,0);
    return c;
}

