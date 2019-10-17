//
// Created by Ethan on 10/11/2019.
//

#include "sun_utils.h"
#include <iostream>
#include <math.h>
#include "../../eigen-git-mirror/Eigen/Dense"
#include <../../pybind11/include/pybind11/pybind11.h>
using namespace Eigen;
using namespace std;
namespace py = pybind11;
VectorXd sun_position(double MJD);

int main(){

    double MJD = 58827.53750000009;
    VectorXd ret = sun_position(MJD);
    cout << ret << endl;
    return 0;
}
VectorXd sun_position(double MJD) {
    /*
    Gives the Greenwich Mean Sidereal Time
    Inputs :
    MJD - Modified Julian Day
    Outputs :
    GMST - Greenwich Mean Sidereal Time
    */
    double JD = MJD + 2400000.5;
    double OplusW = 282.94;
    double T = (JD - 2451545.0) / 36525;

    double M = 357.5256 + 35999.049 * T;

    double lon = OplusW + M + 6892 / 3600 * sin(M) + 72 / 3600 * sin(2*M);
    double r_mag = (149.619 - 2.499 * cos(M) - 0.021 * cos(2*M)) * pow(10, 6);

    double epsilon = 23.43929111;
    VectorXd r_vec(3);
    r_vec(0) = r_mag * cos(lon);
    r_vec(1) = r_mag * sin(lon);
    r_vec(2) = r_mag * sin(lon) * sin(epsilon);

    return r_vec;
}

PYBIND11_MODULE(sun_utils_cpp, m) {
    m.doc() = "pybind11 example plugin"; // optional module docstring

    m.def("sun_position", &sun_position, "A function which returns the sun position");
}