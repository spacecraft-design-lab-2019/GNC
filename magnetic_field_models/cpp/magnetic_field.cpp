//
// Created by ayotundedemuren on 10/28/19.
//

#include "magnetic_field.h"
#include <math.h>
#include <iostream>
#include "../../eigen-git-mirror/Eigen/Dense"

using namespace Eigen;

// input coordinates


int main() {
    // radius of earth
    double R_E = 6378.1;
    // gauss coefficients (decrement n by 1 because of zero-indexing)
    MatrixXd g(3,4);
    g(0,0) = -29442;
    g(0,1) = -1501.0;
    g(0,2) = 0;
    g(0,3) = 0 ;
    g(1,0) = -2445.1;
    g(1,1) = 3012.9;
    g(1,2) = 1676.7;
    g(1,3) = 0;
    g(2,0) = 1350.7;
    g(2,1) = -2352.3;
    g(2,2) = 1225.6;
    g(2,3) = 582.0;
    g(3,0) = 0;
    g(3,1) = 0;
    g(3,2) = 0;
    g(3,3) = 0;

    MatrixXd h(3,4);
    h(0,0) = 0;
    h(0,1) = 4797.1;
    h(0,2) = 0;
    h(0,3) = 0;
    h(1,0) = 0;
    h(1,1) = -2845.6;
    h(1,2) = -641.9;
    h(1,3) = 0;
    h(2,0) = 0;
    h(2,1) = -115.3;
    h(2,2) = 244.9;
    h(2,3) = -538.4;
    h(3,0) = 0;
    h(3,1) = 0;
    h(3,2) = 0;
    h(3,3) = 0;
    return 1;
}

VectorXd get_magnetic_field(double lat, double lon, double alt)
{
    VectorXd B_vec(3);
    B_vec(0) = 0;
    B_vec(1) = 0;
    B_vec(2) = 0;
    return B_vec;
}



















