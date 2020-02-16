//
// Created by Paul on 11/5/2019.
//

#ifdef __cplusplus

#include "../../eigen-git-mirror/Eigen/Dense"
using namespace Eigen;

Vector3d detumble_B_cross(Vector3d omega, Vector3d B, double k);
Vector3d detumble_B_cross_bang_bang(Vector3d omega, Vector3d B, double k, Vector3d max_dipoles);
Vector3d detumble_B_cross_directional(Vector3d omega, Vector3d B, double k, Vector3d max_dipoles);
Vector3d detumble_B_dot(Vector3d B, Vector3d B_dot, double k);
Vector3d detumble_B_dot_bang_bang(Vector3d B_dot, Vector3d max_dipoles);
Vector3d get_bias_estimate(MatrixXd B_mat);
Vector3d get_B_dot(Vector3d B1, Vector3d B2, double dt);

extern "C" {
#endif /* __cplusplus */

void detumble_B_dot_C(double* B_dot, double* max_dipoles, double* commanded_dipole);

#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */