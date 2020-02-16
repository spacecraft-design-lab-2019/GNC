#include <iostream>
#include "MEKF.hpp"
#include <cmath>
#include <stdio.h>

using namespace std;
using Eigen::MatrixXd;

void Innovation(MatrixXd R, MatrixXd rN, MatrixXd rB, MatrixXd P, MatrixXd V, MatrixXd C, MatrixXd &z, MatrixXd &S){
    MatrixXd Q(6,6);
    Q(seq(0,2),seq(0,2)) = R;
    Q(seq(3,5),seq(3,5)) = R;
    z = rB - Q*rN;
    S = C*P*C.transpose()+V;

}


