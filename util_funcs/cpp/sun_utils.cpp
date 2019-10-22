//
// Created by Ethan on 10/11/2019.
//

#include "sun_utils.h"
#define _USE_MATH_DEFINES
#include <iostream>
#include <math.h>
#include "../../eigen-git-mirror/Eigen/Dense"
#include <../../pybind11/include/pybind11/pybind11.h>
#include <../../pybind11/include/pybind11/eigen.h>
using namespace Eigen;
using namespace std;
namespace py = pybind11;
VectorXd sun_position(double MJD);
VectorXd sat_sun_vect(VectorXd r, double MJD);
int main(){

    double MJD = 51622.0;
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
    const double deg2rad = M_PI / 180.0;
    const double rad2deg = 180.0 / M_PI;
    double JD = MJD + 2400000.5;
    double OplusW = 282.94;
    double T = (JD - 2451545.0) / 36525;

    double M = (357.5256 + 35999.049 * T) * deg2rad;

    double lon = (OplusW + rad2deg * M + 6892 / 3600 * sin(M) + 72 / 3600 * sin(2*M)) * deg2rad;
    double r_mag = (149.619 - 2.499 * cos(M) - 0.021 * cos(2*M)) * pow(10, 6);
    cout << lon;
    cout << r_mag;
    double epsilon = deg2rad * 23.43929111;
    VectorXd r_vec(3);
    r_vec(0) = r_mag * cos(lon);
    r_vec(1) = r_mag * sin(lon);
    r_vec(2) = r_mag * sin(lon) * sin(epsilon);

    return r_vec;
}

VectorXd sat_sun_vect(VectorXd r, double MJD){
    /*
    Returns the unit vector from the satellite to the Sun in ECI coordinates
    Inputs:
    r - ECI position of the satellite
    MJD - Julian Day (J2000) as a Real Number
    Outputs:
    r_sat_sun - numpy array giving the unit direction to the Sun from the satellite
     */
    VectorXd r_sun = sun_position(MJD);
    VectorXd r_sat_sun = r_sun - r;
    VectorXd r_sat_sun_unit = r_sat_sun / r_sat_sun.norm();

    return r_sat_sun_unit;
}


PYBIND11_MODULE(sun_utils_cpp, m) {
    m.doc() = "Sun utilities"; // optional module docstring

    m.def("sun_position", &sun_position, "A function which returns the sun position");
    m.def("sat_sun_vect", &sat_sun_vect, "Returns the unit vector to the Sun from the satellite position in the inertial frame");
}