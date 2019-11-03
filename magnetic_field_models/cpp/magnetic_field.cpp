//
// Created by ayotundedemuren on 10/28/19.
//

#include "magnetic_field.h"
#include <math.h>
#include <iostream>
#include "../../eigen-git-mirror/Eigen/Dense"

using namespace Eigen;


int main() {
        return 0;
}

VectorXd get_magnetic_field(double lat, double lon, double alt, int year)
{
    // gauss coefficients (decrement n by 1 because of zero-indexing)
    // IGRF 2015 values
    MatrixXd g(3,4);
    g(0,0) = -29442; g(0,1) = -1501.0; g(0,2) = 0; g(0,3) = 0;
    g(1,0) = -2445.1; g(1,1) = 3012.9; g(1,2) = 1676.7; g(1,3) = 0;
    g(2,0) = 1350.7; g(2,1) = -2352.3; g(2,2) = 1225.6; g(2,3) = 582.0;

    MatrixXd h(3,4);
    h(0,0) = 0; h(0,1) = 4797.1; h(0,2) = 0; h(0,3) = 0;
    h(1,0) = 0; h(1,1) = -2845.6; h(1,2) = -641.9; h(1,3) = 0;
    h(2,0) = 0; h(2,1) = -115.3; h(2,2) = 244.9; h(2,3) = -538.4;

    //2015-2020 secular variation values
    MatrixXd g_sv(3,4);
    g_sv(0,0) = 10.3; g_sv(0,1) = 18.1; g_sv(0,2) = 0; g_sv(0,3) = 0 ;
    g_sv(1,0) = -8.7; g_sv(1,1) = -3.3; g_sv(1,2) = 2.1; g_sv(1,3) = 0;
    g_sv(2,0) = 3.4; g_sv(2,1) = -5.5; g_sv(2,2) = -0.7; g_sv(2,3) = -10.1;

    MatrixXd h_sv(3,4);
    h_sv(0,0) = 0; h_sv(0,1) = -26.6; h_sv(0,2) = 0; h_sv(0,3) = 0;
    h_sv(1,0) = 0; h_sv(1,1) = -27.4; h_sv(1,2) = -14.1; h_sv(1,3) = 0;
    h_sv(2,0) = 0; h_sv(2,1) = 8.2; h_sv(2,2) = -0.4; h_sv(2,3) = 1.8;
    const double deg2rad = M_PI / 180.0;
    lat = lat*deg2rad;
    lon = lon*deg2rad;
    // radius of earth
    double a = 6378.1;
    double r = a + alt;
    // year since 2015 for secular variation
    int dt = year - 2015;
    // magnetic field components
    double B_x = 0; double B_y = 0; double B_z = 0;

    // Associated Legendre functions and derivatives (wrt lat)
    MatrixXd P(3,4);
    P(0,0) = cos(lat);
    P(0,1) = -sqrt(1 - pow(cos(lat),2));
    P(0,2) = 0;
    P(0,3) = 0;
    P(1,0) = (1/2)*(3*pow(cos(lat), 2) - 1);
    P(1,1) = -3*cos(lat)*sqrt(1 - pow(cos(lat), 2));
    P(1,2) = 3*(1 - pow(cos(lat), 2));
    P(1,3) = 0;
    P(2,0) = (1/2)*(5*pow(cos(lat), 3) - 3*cos(lat));
    P(2,1) = (-3/2)*(5*pow(cos(lat), 2) - 1)*sqrt(1 - pow(cos(lat), 2));
    P(2,2) = 15*cos(lat)*(1 - pow(cos(lat), 2));
    P(2,3) = -15*pow(1 - pow(cos(lat), 2), (3/2));

    MatrixXd P_d(3,4);
    P_d(0,0) = -sin(lat);
    P_d(0,1) = (-sin(lat)*cos(lat))/(sqrt(1 - pow(cos(lat),2)));
    P_d(0,2) = 0;
    P_d(0,3) = 0;
    P_d(1,0) = -3*sin(lat)*cos(lat);
    P_d(1,1) = (3*sin(lat)*(1 - 2*pow(cos(lat), 2)))/(sqrt(1 - pow(cos(lat),2)));
    P_d(1,2) = 6*sin(lat)*cos(lat);
    P_d(1,3) = 0;
    P_d(2,0) = (-3/2)*sin(lat)*(5*pow(cos(lat), 2) - 1);
    P_d(2,1) = (3*sin(lat)*cos(lat)*(11 - 15*pow(cos(lat), 2)))/(2*sqrt(1 - pow(cos(lat),2)));
    P_d(2,2) = 15*sin(lat)*(3*pow(cos(lat), 2) - 1);
    P_d(2,3) = -45*sin(lat)*cos(lat)*sqrt(1 - pow(cos(lat), 2));

    // IGRF model B field calculation, n/m sets order (3)
    for(int n = 0; n < 3; n++){
        for(int m = 0; m < 4; m++){
            B_x += pow((a/r), n + 1)*((g(n,m) + dt*g_sv(n,m))*cos(m*lon) +
                    (h(n,m) + dt*h_sv(n,m))*sin(m*lon))*P_d(n,m);
            B_y += pow((a/r), n + 1)*(-m*(g(n,m) + dt*g_sv(n,m))*sin(m*lon) +
                    m*(h(n,m) + dt*h_sv(n,m))*cos(m*lon))*P(n,m);
            B_z += (-a*(n+1)*(pow((a/r), n))/(pow(r, 2)))*((g(n,m) + dt*g_sv(n,m))*cos(m*lon) +
                     (h(n,m) + dt*h_sv(n,m))*sin(m*lon))*P(n,m);
        }
    }
    B_x = B_x*(a/r);
    B_y = B_y*(-a/(r*sin(lat)));
    B_z = B_z*a;

    VectorXd B_vec(3);
    B_vec(0) = B_x;
    B_vec(1) = B_y;
    B_vec(2) = B_z;


    return B_vec;
}



















