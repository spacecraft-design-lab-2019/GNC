//
// Created by Ethan on 10/23/2019.
//

#include "frame_conversions.h"
#include "../../eigen-git-mirror/Eigen/Dense"
#include <../../pybind11/include/pybind11/pybind11.h>
#include <../../pybind11/include/pybind11/eigen.h>
namespace py = pybind11;
using namespace Eigen;
using namespace std;

MatrixXd eci2ecef(double GMST);
std::tuple<double, double, double> ecef2lla(Vector3d);
MatrixXd ecef2enu(double lat, double lon);

int main(){
    return 0;
}

MatrixXd eci2ecef(double GMST){
    /*
    Rotation matrix from ECI to ECEF coordinates
    Inputs:
    GMST - Greenwich Mean Sidereal Time
    Outputs:
    R - Rotation matrix from ECI to ECEF
    */
    MatrixXd R(3,3);
    R <<  cos(GMST), sin(GMST), 0,
            -sin(GMST), cos(GMST), 0,
            0, 0 ,1;

    return R;
}

MatrixXd ecef2enu(double lat, double lon){
    /*
    Rotation matrix from ECEF to ENU coordinates
    Inputs:
    lat - Latitude in radians
    lon - Longitude in radians
    Outputs:
    R - Rotation matrix from ECEF to ENU
    */
    Vector3d ehat, nhat, uhat;
    ehat << -sin(lon), cos(lon), 0;
    nhat << -sin(lat) * cos(lon), -sin(lat) * sin(lon), cos(lat);
    uhat << cos(lat)*cos(lon), cos(lat)*sin(lon), sin(lat);

    MatrixXd R(3, 3);
    R << ehat, nhat, uhat;
    MatrixXd R_t(3, 3);
    R_t = R.transpose();

    return R_t;
}

std::tuple<double, double, double> ecef2lla(Vector3d r){
    /*
     * Gives longitude, latitude, and altitude from ECEF position vector
    Inputs:
    r - Position vector in ECEF
    Outputs:
    lat - geocentric latitude in radians
    lon - longitude in radians
    alt - altitude
     */
    double R_earth = 6378.1;
    double lat = asin(r(2) / r.norm());
    double lon = atan2(r(1), r(0));
    double alt = r.norm() - R_earth;
    return std::make_tuple(lat, lon, alt);
}



PYBIND11_MODULE(frame_conversions_cpp, m) {
    m.doc() = "Frame Conversions"; // optional module docstring

    m.def("eci2ecef", &eci2ecef, "Gives rotation matrix from ECI2ECEF");
    m.def("ecef2lla", &ecef2lla, "Converts position in ECEF to lat, long, alt");
    m.def("ecef2enu", &ecef2enu,  "Gives rotation matrix from ECEF2enu using long and lat");
}