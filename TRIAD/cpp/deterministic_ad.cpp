//
// Created by Ethan on 10/16/2019.
//

#include "deterministic_ad.h"
#include <iostream>
#include "../../eigen-git-mirror/Eigen/Dense"
#include <../../pybind11/include/pybind11/pybind11.h>
#include <../../pybind11/include/pybind11/eigen.h>
namespace py = pybind11;
using namespace Eigen;
using namespace std;

MatrixXd triad_ad(MatrixXd M, MatrixXd V);

int main(){
    MatrixXd M(3,2);
    M << -0.0561111, 0.47183533,
    -0.84110217,  0.39004599,
    0.53795788, -0.79071838;
    MatrixXd V(3,2);
    V << -0.84639796, 0.39004599,
    0.06513016, -0.47183533,
    0.52855326, -0.79071838;
    MatrixXd ret(3,3);
    ret = triad_ad(M, V);
    cout << ret << endl;
    return 0;
}

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

PYBIND11_MODULE(triad_cpp, m) {
    m.doc() = "AD Methods"; // optional module docstring

    m.def("triad_ad", &triad_ad, "Performs TRIAD algorithm");
}

