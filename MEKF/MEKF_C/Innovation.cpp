#include <iostream>
#include "eigen-git-mirror/Eigen/Dense"


using Eigen::MatrixXd;

MatrixXd innovation(MatrixXd R, MatrixXd rN, MatrixXd rB, MatrixXd P, MatrixXd V, MatrixXd C){
    MatrixXd Q(6,6);
    Q(0:3,0:3) = R;
    Q(3:6,3:6) = R;
    z = rB - Q*rN;
    S = C*P*C.transpose()+V;
    return z,S;
}


