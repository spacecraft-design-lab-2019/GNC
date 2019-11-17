#include <iostream>
#include "MEKF.hpp"
#include <cmath>
#include <stdio.h>

using namespace std;
using Eigen::MatrixXd;

MatrixXd DCM2q(MatrixXd A) {
    double B1, B2, B3, B4;
    B1 = sqrt((1.0 / 4.0) * (1.0 + 2.0 * A(0, 0) - trace(A)));
    B2 = sqrt((1.0 / 4.0) * (1.0 + 2.0 * A(1, 1) - trace(A)));
    B3 = sqrt((1.0 / 4.0) * (1.0 + 2.0 * A(2, 2) - trace(A)));
    B4 = sqrt((1.0 / 4.0) * (1.0 + trace(A)));

    double qs, q1, q2, q3;

    if ((B1 > B2) and (B1 > B3) and (B1 > B4)) {
        q1 = B1;
        q2 = (A(0, 1) + A(1, 0)) / (4 * B1);
        q3 = (A(0, 2) + A(2, 0)) / (4 * B1);
        qs = (A(1, 2) - A(2, 1)) / (4 * B1);
    }else if ((B2 > B1) and (B2 > B3) and (B2 > B4)) {
        q1 = (A(0, 1) + A(1, 0)) / (4 * B2);
        q2 = B2;
        q3 = (A(1, 2) + A(2, 1)) / (4 * B2);
        qs = (A(2, 0) - A(0, 2)) / (4 * B2);
    } else if ((B3 > B1) and (B3 > B2) and (B3>B4)){
        q1 = (A(2,0)+A(0,2))/(4*B3);
        q2 = (A(2,1)+A(1,2))/(4*B3);
        q3 = B3;
        qs = (A(0,1)-A(1,0))/(4*B3);
    } else{
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