#include <iostream>
#include "MEKF.hpp"
#include <cmath>
#include <stdio.h>

using namespace std;
using Eigen::MatrixXd;

// taken from utility functions
MatrixXd triad_ad(MatrixXd M, MatrixXd V) {
    /*
    Gives rotation matrix from inertial to body frame
    Inputs :
    M - Matrix where each column is a measurement vector in the body frame - Note, most accurate measurement should be in the first column
    V - Matrix where each column is a modeled vector of the corresponding measurement vector from M, in the inertial frame
    Outputs:
    R - rotation matrix from inertial to body frame
    */
    MatrixXd R(3,3);
    if (M.outerSize() == 2 && V.outerSize() == 2)
    {
        Vector3d m1 = M.col(0);
        Vector3d mtemp = M.col(1);
        Vector3d m2 = m1.cross(mtemp);
        Vector3d m3 = m1.cross(m2);

        Vector3d v1 = V.col(0);
        Vector3d vtemp = V.col(1);
        Vector3d v2 = v1.cross(vtemp);
        Vector3d v3 = v1.cross(v2);

        MatrixXd Rtemp(3,3);
        MatrixXd Vtemp(3, 3);
        Rtemp << m1, m2, m3;
        Vtemp << v1, v2, v3;
        R = Rtemp * Vtemp.inverse();
    }
    else
    {
        R = M * V.completeOrthogonalDecomposition().pseudoInverse();
    }
    return R;
}
