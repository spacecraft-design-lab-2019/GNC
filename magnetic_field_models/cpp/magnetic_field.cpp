//
// Created by ayotundedemuren on 10/28/19.
//

#include "magnetic_field.h"
#include <math.h>
#include <iostream>
#include "../../eigen-git-mirror/Eigen/Dense"
#include <../../pybind11/include/pybind11/pybind11.h>
#include <../../pybind11/include/pybind11/eigen.h>

using namespace std;
using namespace Eigen;
namespace py = pybind11;

//function declaration
//VectorXd get_magnetic_field(double lat, double lon, double alt, int year);


VectorXd get_magnetic_field(double lat, double lon, double alt, double year, int order);
MatrixXd get_P_coefficients(double x, int order);
MatrixXd get_g_coefficients();
MatrixXd get_h_coefficients();
MatrixXd get_g_sv_coefficients();
MatrixXd get_h_sv_coefficients();
MatrixXd get_Pd_coefficients(MatrixXd P, double x, int order);


const MatrixXd g = get_g_coefficients();
const MatrixXd h = get_h_coefficients();
const MatrixXd g_sv = get_g_sv_coefficients();
const MatrixXd h_sv = get_h_sv_coefficients();

int main() {
    // double lat = 45; double lon = 45; double alt = 400; double year = 2015;
    // MatrixXd P = get_P_coefficients(cos((90.0 - lat) * M_PI / 180.0));
    // MatrixXd Pd = get_Pd_coefficients(P, cos(M_PI/2.0 - lat*M_PI/180.0));
    // VectorXd B = get_magnetic_field(lat, lon, alt, year);
    // cout << P << endl;
    // cout << Pd << endl;
    // cout << B << endl;
    return 0;

}


VectorXd get_magnetic_field(double lat, double lon, double alt, double year, int order){

    /*
    lat is geocentric latitude in degrees
    lon is longitude in degrees
    alt is altitude in km
    year is the fractional year (include months/days essentially)

    outputs B vector in NED
    */
    const double deg2rad = M_PI / 180.0;
    lat = 90 - lat;
    lat = lat*deg2rad;
    lon = lon*deg2rad;
    // radius of earth
    double a = 6371.2;
    double r = a + alt;
    // year since 2015 for secular variation
    double dt = year - 2015.0;
    // magnetic field components
    double B_r = 0; double B_lat = 0; double B_lon = 0;

    MatrixXd P = get_P_coefficients(cos(lat), order);
    MatrixXd Pd = get_Pd_coefficients(P, cos(lat), order);

    // IGRF model B field calculation, n/m sets order (5)

    for(int n = 1; n < order+1; n++){
        for(int m = 0; m < order+1; m++){

            double coef = pow((a/r), n + 2) * ((g(n, m) + dt * g_sv(n, m)) * cos(m * lon) +
                    (h(n, m) + dt * h_sv(n, m)) * sin(m * lon));
            // Radial component
            B_r +=  (n + 1) * coef * P(n, m);
            // Colatitudinal component
            B_lat -= coef * Pd(n, m);
            // Address singularity at colatitude of 0
            if(sin(lat) == 0){
                // Longitudinal component
                B_lon += -cos(lat) * pow((a/r), n + 2) * (- (g(n, m) + dt * g_sv(n, m)) * sin(m * lon) +
                        (h(n, m) + dt * h_sv(n, m)) * cos(m * lon)) * Pd(n, m);
            }else{
                B_lon += (-1/sin(lat)) * pow((a/r), n + 2) * m * (-(g(n, m) + dt * g_sv(n, m)) * sin(m * lon) +
                        (h(n, m) + dt * h_sv(n, m)) * cos(m * lon)) * P(n, m);
            }
        }
    }


    VectorXd B_vec(3);
    // NED (North, East, Down) coordinate frame
    B_vec(0) = -B_lat;
    B_vec(1) = B_lon;
    B_vec(2) = -B_r;


    return B_vec;
}


MatrixXd get_g_coefficients(){
    Matrix<double, 11, 11> g;

    g <<      0.0,     0.0,    0.0,    0.0,    0.0,  0.0,   0.0,   0.0,  0.0,   0.0,  0.0,
         -29442.0, -1501.0,    0.0,    0.0,    0.0,  0.0,   0.0,   0.0,  0.0,   0.0,  0.0,
          -2445.1,  3012.9, 1676.7,    0.0,    0.0,  0.0,   0.0,   0.0,  0.0,   0.0,  0.0,
           1350.7, -2352.3, 1225.6,  582.0,    0.0,  0.0,   0.0,   0.0,  0.0,   0.0,  0.0,
            907.6,   813.7,  120.4, -334.9,   70.4,  0.0,   0.0,   0.0,  0.0,   0.0,  0.0,
           -232.6,   360.1,  192.4, -140.9, -157.5,  4.1,   0.0,   0.0,  0.0,   0.0,  0.0,
             70.0,    67.7,   72.7, -129.9,  -28.9, 13.2, -70.9,   0.0,  0.0,   0.0,  0.0,
             81.6,   -76.1,   -6.8,   51.8,   15.0,  9.4,  -2.8,   6.8,  0.0,   0.0,  0.0,
             24.2,     8.8,  -16.9,   -3.2,  -20.6, 13.4,  11.7, -15.9, -2.0,   0.0,  0.0,
             5.4,      8.8,    3.1,   -3.3,    0.7, -13.3, -0.1,   8.7,  9.1, -10.5,  0.0,
            -1.9,     -6.3,    0.1,    0.5,   -0.5,   1.8, -0.7,   2.1,  2.4,  -1.8, -3.6;

    return g;
}

MatrixXd get_h_coefficients(){
    Matrix<double, 11, 11> h;

    h << 0.0,     0.0,    0.0,    0.0,    0.0,   0.0,   0.0,  0.0,  0.0,  0.0,  0.0,
         0.0,  4797.1,    0.0,    0.0,    0.0,   0.0,   0.0,  0.0,  0.0,  0.0,  0.0,
         0.0, -2845.6, -641.9,    0.0,    0.0,   0.0,   0.0,  0.0,  0.0,  0.0,  0.0,
         0.0,  -115.3,  244.9, -538.4,    0.0,   0.0,   0.0,  0.0,  0.0,  0.0,  0.0,
         0.0,   283.3, -188.7,  180.9, -329.5,   0.0,   0.0,  0.0,  0.0,  0.0,  0.0,
         0.0,    47.3,  197.0, -119.3,   16.0, 100.2,   0.0,  0.0,  0.0,  0.0,  0.0,
         0.0,   -20.8,   33.2,   58.9,  -66.7,   7.3,  62.6,  0.0,  0.0,  0.0,  0.0,
         0.0,   -54.1,  -19.5,    5.7,   24.4,   3.4, -27.4, -2.2,  0.0,  0.0,  0.0,
         0.0,    10.1,  -18.3,   13.3,  -14.6,  16.2,   5.7, -9.1,  2.1,  0.0,  0.0,
         0.0,   -21.6,   10.8,   11.8,   -6.8,  -6.9,   7.8,  1.0, -4.0,  8.4,  0.0,
         0.0,     3.2,   -0.4,    4.6,    4.4,  -7.9,  -0.6, -4.2, -2.8, -1.2, -8.7;

    return h;
}

// 2015-2020 secular variation values
MatrixXd get_g_sv_coefficients() {
    Matrix<double, 11, 11> g_sv;

    g_sv <<  0.0,   0.0,  0.0,   0.0,  0.0,  0.0,  0.0,  0.0, 0.0, 0.0, 0.0,
            10.3,  18.1,  0.0,   0.0,  0.0,  0.0,  0.0,  0.0, 0.0, 0.0, 0.0,
            -8.7,  -3.3,  2.1,   0.0,  0.0,  0.0,  0.0,  0.0, 0.0, 0.0, 0.0,
             3.4,  -5.5, -0.7, -10.1,  0.0,  0.0,  0.0,  0.0, 0.0, 0.0, 0.0,
            -0.7,   0.2, -9.1,   4.1, -4.3,  0.0,  0.0,  0.0, 0.0, 0.0, 0.0,
            -0.2,   0.5, -1.3,  -0.1,  1.4,  3.9,  0.0,  0.0, 0.0, 0.0, 0.0,
            -0.3,  -0.1, -0.7,   2.1, -1.2,  0.3,  1.6,  0.0, 0.0, 0.0, 0.0,
             0.3,  -0.2, -0.5,   1.3,  0.1, -0.6, -0.8,  0.2, 0.0, 0.0, 0.0,
             0.2,   0.0, -0.6,   0.5, -0.2,  0.4,  0.1, -0.4, 0.3, 0.0, 0.0,
             0.0,   0.0,  0.0,   0.0,  0.0,  0.0,  0.0,  0.0, 0.0, 0.0, 0.0,
             0.0,   0.0,  0.0,   0.0,  0.0,  0.0,  0.0,  0.0, 0.0, 0.0, 0.0;

    return g_sv;
}

MatrixXd get_h_sv_coefficients() {
    Matrix<double, 11, 11> h_sv;

    h_sv <<  0.0,   0.0,   0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.0, 0.0, 0.0,
             0.0, -26.6,   0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.0, 0.0, 0.0,
             0.0, -27.4, -14.1,  0.0,  0.0,  0.0,  0.0,  0.0, 0.0, 0.0, 0.0,
             0.0,   8.2,  -0.4,  1.8,  0.0,  0.0,  0.0,  0.0, 0.0, 0.0, 0.0,
             0.0,  -1.3,   5.3,  2.9, -5.2,  0.0,  0.0,  0.0, 0.0, 0.0, 0.0,
             0.0,   0.6,   1.7, -1.2,  3.4,  0.0,  0.0,  0.0, 0.0, 0.0, 0.0,
             0.0,   0.0,  -2.1, -0.7,  0.2,  0.9,  1.0,  0.0, 0.0, 0.0, 0.0,
             0.0,   0.8,   0.4, -0.2, -0.3, -0.6,  0.1, -0.2, 0.0, 0.0, 0.0,
             0.0,  -0.3,   0.3,  0.1,  0.5, -0.2, -0.3,  0.3, 0.0, 0.0, 0.0,
             0.0,   0.0,   0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.0, 0.0, 0.0,
             0.0,   0.0,   0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.0, 0.0, 0.0;

    return h_sv;
}

double doub_fact(int x){
    double val = 1.0;
    for(int k = 0; k <= ceil(x/2 - 1); k++){
        val = val * (x - 2 * k);
    }
    return val;
}

MatrixXd get_P_coefficients(double x, int order){

    MatrixXd P(order+2, order+2);
    P = MatrixXd::Zero(order+2, order+2);
    P(0,0) = 1.0;
    P(1, 1) = sqrt(1 - pow(x, 2));
    for(int n = 1; n < order+2; n++){
        for(int m = 0; m < order+2; m++){

            if(n != 1 || m != 1){
                int prev2 = n - 2;
                if(prev2 < 0){
                    prev2 = 0;
                }
                if(m < n){
                    P(n, m) = (2 * n - 1) / sqrt(pow(n, 2) - pow(m, 2)) * x * P(n - 1, m) - sqrt((pow(n-1, 2) - pow(m, 2)) / (pow(n, 2) - pow(m, 2))) * P(prev2, m);
                }else{
                    P(n, m) = sqrt(1.0 - 1.0/(2.0*m)) * sqrt(1 - pow(x, 2)) * P(n - 1, m - 1);
                }
            }
        }
    }
    return P;
}

MatrixXd get_Pd_coefficients(MatrixXd P, double x, int order){

    MatrixXd Pd(order+1, order+1);
    Pd = MatrixXd::Zero(order+1, order+1);

    Pd(1,1) = x;
    for(int n = 1; n < order+1; n++){
        for(int m = 0; m < order+1; m++){
            if(n != 1 || m != 1){
                if(m < n){
                    int prev2 = n - 2;
                    if(prev2 < 0){
                        prev2 = 0;
                    }
                    Pd(n, m) = (2 * n - 1) / sqrt(pow(n, 2) - pow(m, 2)) * (x * Pd(n-1, m) -
                                                                            sqrt(1 - pow(x, 2)) * P(n-1, m)) - sqrt((pow(n-1, 2) - pow(m, 2))/(pow(n, 2) - pow(m, 2))) * Pd(prev2, m);
                }else{
                    Pd(n, m) = sqrt(1.0 - 1.0 / (2.0 * m)) * (sqrt(1 - pow(x, 2)) * Pd(n - 1, m - 1) + x * P(n - 1, m - 1));
                }
            }
        }
    }
    return Pd;
}


PYBIND11_MODULE(magnetic_field_cpp, m) {
    m.doc() = "Magnetic Field"; // optional module docstring

    m.def("get_magnetic_field", &get_magnetic_field, "Gives mag field in NED at the given lat lon alt and year, use geocentric");
    m.def("get_P_coefficients", &get_P_coefficients);
    m.def("get_Pd_coefficients", &get_Pd_coefficients);
    m.def("get_g_coefficients", &get_g_coefficients);
    m.def("get_h_coefficients", &get_h_coefficients);
    m.def("get_g_sv_coefficients", &get_g_sv_coefficients);
    m.def("get_h_sv_coefficients", &get_h_sv_coefficients);
}










