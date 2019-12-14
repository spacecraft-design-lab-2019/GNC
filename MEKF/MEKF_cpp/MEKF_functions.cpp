#include <iostream>
#include "MEKF.hpp"
#include <cmath>
#include <stdio.h>
#include "../../eigen-git-mirror/Eigen/Dense"


using namespace Eigen;
using namespace std;
/* Define Utilities for MEKF */

MatrixXd quatmult(MatrixXd q1, MatrixXd q2);
MatrixXd cross(MatrixXd a, MatrixXd b);
MatrixXd SkewSymmetric(double a1, double a2, double a3);
MatrixXd quat2dcm(MatrixXd q);

/* Predict Step */

MatrixXd predict_xn(MatrixXd xk, MatrixXd Pk, MatrixXd w, double dt, MatrixXd W){
    MatrixXd q(4,1), b(3,1);
    q(seq(0,3),0) = xk(seq(0,3),0);
    b(seq(0,2),0) = xk(seq(4,6),0);

    MatrixXd u(3,1);
    u = w + b;
    double theta;
    theta = dt*sqrt(w(0,0)*w(0,0)+w(1,0)*w(1,0)+w(2,0)*w(2,0));
    MatrixXd r0(3,1), r(3,1);
    r0 = (u-b);
    r(0,0) = r0(0,0)/(sqrt(r0(0,0)*r0(0,0)+r0(1,0)*r0(1,0)+r0(2,0)*r0(2,0)));
    r(1,0) = r0(1,0)/(sqrt(r0(0,0)*r0(0,0)+r0(1,0)*r0(1,0)+r0(2,0)*r0(2,0)));
    r(2,0) = r0(2,0)/(sqrt(r0(0,0)*r0(0,0)+r0(1,0)*r0(1,0)+r0(2,0)*r0(2,0)));
    MatrixXd s(4,1); // q2
    s(0,0) = cos(theta/2.0);
    s(1,0) = r(0,0)*sin(theta/2.0);
    s(2,0) = r(1,0)*sin(theta/2.0);
    s(3,0) = r(2,0)*sin(theta/2.0);

    MatrixXd xn(7,1);
    xn(0,0) = q(0,0)*s(0,0)-(q(1,0)*s(1,0)+q(2,0)*s(2,0)+q(3,0)*s(3,0));


    // cross product of s and q
    MatrixXd cross_sq(3,1);
    cross_sq(0,0) = q(2,0)*s(3,0)-s(2,0)*q(3,0);
    cross_sq(1,0) = q(3,0)*s(1,0)-s(3,0)*q(1,0);
    cross_sq(2,0) = q(1,0)*s(2,0)-s(1,0)*q(2,0);

    xn(1,0) = q(0,0)*s(1,0)+s(0,0)*q(1,0)+cross_sq(0,0);
    xn(2,0) = q(0,0)*s(2,0)+s(0,0)*q(2,0)+cross_sq(1,0);
    xn(3,0) = q(0,0)*s(3,0)+s(0,0)*q(3,0)+cross_sq(2,0);
    xn(seq(4,6),0) = b;

    return xn;
}

MatrixXd predict_Pn(MatrixXd xk, MatrixXd Pk, MatrixXd w, double dt,  MatrixXd W){
    MatrixXd q(4,1), b(3,1);
    q(seq(0,3),0) = xk(seq(0,3),0);
    b(seq(0,2),0) = xk(seq(4,6),0);

    MatrixXd u(3,1);
    u = w + b;

    double theta;
    theta = dt*sqrt(w(0,0)*w(0,0)+w(1,0)*w(1,0)+w(2,0)*w(2,0));
    
    MatrixXd r0(3,1), r(3,1);
    r0 = (u-b);

    r(0,0) = r0(0,0)/(sqrt(r0(0,0)*r0(0,0)+r0(1,0)*r0(1,0)+r0(2,0)*r0(2,0)));
    r(1,0) = r0(1,0)/(sqrt(r0(0,0)*r0(0,0)+r0(1,0)*r0(1,0)+r0(2,0)*r0(2,0)));
    r(2,0) = r0(2,0)/(sqrt(r0(0,0)*r0(0,0)+r0(1,0)*r0(1,0)+r0(2,0)*r0(2,0)));
    
    MatrixXd s(4,1); // q2
    s(0,0) = cos(theta/2.0);
    s(1,0) = r(0,0)*sin(theta/2.0);
    s(2,0) = r(1,0)*sin(theta/2.0);
    s(3,0) = r(2,0)*sin(theta/2.0);

    MatrixXd V(3,4);
    V = MatrixXd::Zero(3,4);
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

    MatrixXd A(6,6);
    A(seq(0,2),seq(0,2)) = V*Ls.transpose()*Rs*V.transpose();
    A(seq(0,2),seq(3,5)) = -0.5*dt*MatrixXd::Identity(3,3);
    A(seq(3,5),seq(3,5)) = MatrixXd::Identity(3,3);
    A(seq(3,5),seq(0,2)) = MatrixXd::Zero(3,3);
    MatrixXd Pn(6,6);
    Pn = A*Pk*A.transpose()+W;
    return Pn;
}

/* Measurement Step */

MatrixXd measurement(MatrixXd q, MatrixXd rN){
    MatrixXd R(3,3);
    R = quat2dcm(q);
    MatrixXd rB1(3,1);
    rB1 = R*rN(seq(0,2),0);
    MatrixXd rB2(3,1);
    rB2 = R*rN(seq(3,5),0);

    MatrixXd C(6,6);

    C(seq(0,2),seq(0,2)) = 2*SkewSymmetric(rB1(0),rB1(1),rB1(2));
    C(seq(3,5),seq(0,2)) = 2*SkewSymmetric(rB2(0),rB2(1),rB2(2));
    C(all,seq(3,5)) = MatrixXd::Zero(6,3);
    return C;
}
/* Innovation Step */

MatrixXd innovation(MatrixXd R, MatrixXd rN, MatrixXd rB){
    MatrixXd Q(6,6);
    Q(seq(0,2),seq(0,2)) = R;
    Q(seq(3,5),seq(3,5)) = R;
    MatrixXd z(6,1);
    z = rB - Q*rN;
    return z;
}

/* Update Step */


MatrixXd update_xk(MatrixXd L, MatrixXd z, MatrixXd xn, MatrixXd Pn, MatrixXd V, MatrixXd C){
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
    MatrixXd xk(7,1);
    xk(seq(0,3),0) = q;
    xk(seq(4,6),0) = b;
    return xk;
    
}

MatrixXd update_Pk(MatrixXd L, MatrixXd z, MatrixXd xn, MatrixXd Pn, MatrixXd V, MatrixXd C){
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
    MatrixXd xk(7,1);
    xk(seq(0,3),0) = q;
    xk(seq(4,6),0) = b;
    MatrixXd Pk(6,6);
    Pk = (MatrixXd::Identity(6,6)-L*C)*Pn*(MatrixXd::Identity(6,6)-L*C).transpose() + L*V*L.transpose();
    return Pk;
}

/* Other utilities for MEKF */

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
    DCM = (pow(q1,2.0)-q2*q2-q3*q3-q4*q4)*MatrixXd::Identity(3,3)-2.0*q1*Q+2.0*q_vec*q_vec.transpose();
    return DCM;
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



