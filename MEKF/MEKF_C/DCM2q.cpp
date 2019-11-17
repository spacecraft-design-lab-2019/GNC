#include <iostream>
#include "MEKF.hpp"

using Eigen::MatrixXd;
double trace(MatrixXd A);

MatrixXd DCM2q(MatrixXd A){
    double B1, B2, B3, B4;
    B1 = sqrt(1/4*(1+2*A(0,0)-trace(A)));
    B2 = sqrt(1/4*(1+2*A(1,1)-trace(A)));
    B3 = sqrt(1/4*(1+2*A(2,2)-trace(A)));
    B4 = sqrt(1/4*(1+trace(A)));
    MatrixXd B(4,1);
    B(0,0) = B1;
    B(1,0) = B2;
    B(2,0) = B3;
    B(3,0) = B4;
    double qs,q1,q2,q3;
    if (B1 == max(B)) {
        q1 = B1;
        q2 = (A(0, 1) + A(1, 0)) / (4 * B1);
        q3 = (A(0, 2) + A(2, 0)) / (4 * B1);
        qs = (A(1, 2) - A(2, 1)) / (4 * B1);
    }
    else if (B2 == max(B)){
        q1 = (A(0,1)+A(1,0))/(4*B2);
        q2 = B2;
        q3 = (A(1,2)+A(2,1))/(4*B2);
        qs = (A(2,0)-A(0,2))/(4*B2);
    }
    else if (B3 == max(B)){
        q1 = (A(2,0)+A(0,2))/(4*B3);
        q2 = (A(2,1)+A(1,2))/(4*B3);
        q3 = B3;
        qs = (A(0,1)-A(1,0))/(4*B3);
    }
    else{
        q1 = (A(1,2)-A(2,1))/(4*B4);
        q2 = (A(2,0)-A(0,2))/(4*B4);
        q3 = (A(0,1)-A(1,0))/(4*B4);
        qs = B4;
    }
    MatrixXd q(4,1);
    q(0,0) = qs;
    q(1,0) = q1;
    q(2,0) = q2;
    q(3,0) = q3;
    double qnorm;
    qnorm = sqrt(pow(qs,2)+pow(q1,2)+pow(q2,2)+pow(q3,2));
    for (int i = 0; i < 4; i++) {
        q(i) = q(i) / qnorm;
    }
    return q;
}

double trace(MatrixXd A){
    a = A(0,0) + A(1,1) + A(2,2);
    return a;
}
