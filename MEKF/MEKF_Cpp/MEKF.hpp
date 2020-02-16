

#ifndef GNC_MEKF_HPP
#define GNC_MEKF_HPP

#include "../../eigen-git-mirror/Eigen/Dense"

using namespace Eigen;

MatrixXd predict_xn(MatrixXd xk, MatrixXd Pk, MatrixXd w, double dt, MatrixXd W);
MatrixXd predict_Pn(MatrixXd xk, MatrixXd Pk, MatrixXd w, double dt, MatrixXd W);
MatrixXd measurement(MatrixXd q, MatrixXd rN);
MatrixXd innovation(MatrixXd R, MatrixXd rN, MatrixXd rB);
MatrixXd update_xk(MatrixXd L, MatrixXd z, MatrixXd xn, MatrixXd Pn, MatrixXd V, MatrixXd C);
MatrixXd update_Pk(MatrixXd L, MatrixXd z, MatrixXd xn, MatrixXd Pn, MatrixXd V, MatrixXd C);


MatrixXd DCM2q(MatrixXd A);
double trace(MatrixXd A);
MatrixXd triad_ad(MatrixXd M, MatrixXd V);

#endif //GNC_MEKF_HPP
