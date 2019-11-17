//
// Created by jayde on 11/16/2019.
//

#ifndef GNC_MEKF_HPP
#define GNC_MEKF_HPP

#include "../../eigen-git-mirror/Eigen/Dense"

using namespace Eigen;
void Innovation(MatrixXd R, MatrixXd rN, MatrixXd rB, MatrixXd P, MatrixXd V, MatrixXd C, MatrixXd &z, MatrixXd &S);
void measurement(MatrixXd q, MatrixXd rN, MatrixXd &y, MatrixXd &R, MatrixXd &C);
void update(MatrixXd L, MatrixXd z, MatrixXd xn, MatrixXd Pn, MatrixXd V, MatrixXd C, MatrixXd &x, MatrixXd &P);
void predict(MatrixXd xk, MatrixXd w, double dt, MatrixXd &xn, MatrixXd &A);
MatrixXd DCM2q(MatrixXd A);
double trace(MatrixXd A);
MatrixXd triad_ad(MatrixXd M, MatrixXd V);

#endif //GNC_MEKF_HPP
