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


VectorXd get_magnetic_field(double lat, double lon, double alt, int year);
MatrixXd get_P_coefficients(double x);
MatrixXd get_g_coefficients();
MatrixXd get_h_coefficients();
MatrixXd get_Pd_coefficients(MatrixXd P, double x);
int main() {
    double lat = 45; double lon = 45; double alt = 400; int year = 2015;
    MatrixXd P = get_P_coefficients(cos((90.0 - lat) * M_PI / 180.0));
    MatrixXd Pd = get_Pd_coefficients(P, cos(M_PI/2.0 - lat*M_PI/180.0));
    VectorXd B = get_magnetic_field(lat, lon, alt, year);
    cout << P << endl;
    cout << Pd << endl;
    cout << B << endl;
    return 0;

}

VectorXd get_magnetic_field(double lat, double lon, double alt, int year){


    //2015-2020 secular variation values
    // MatrixXd g_sv(3,4);
    // g_sv(0,0) = 10.3; g_sv(0,1) = 18.1; g_sv(0,2) = 0.0; g_sv(0,3) = 0.0 ;
    // g_sv(1,0) = -8.7; g_sv(1,1) = -3.3; g_sv(1,2) = 2.1; g_sv(1,3) = 0.0;
    // g_sv(2,0) = 3.4; g_sv(2,1) = -5.5; g_sv(2,2) = -0.7; g_sv(2,3) = -10.1;

    // MatrixXd h_sv(3,4);
    // h_sv(0,0) = 0.0; h_sv(0,1) = -26.6; h_sv(0,2) = 0.0; h_sv(0,3) = 0.0;
    // h_sv(1,0) = 0.0; h_sv(1,1) = -27.4; h_sv(1,2) = -14.1; h_sv(1,3) = 0.0;
    // h_sv(2,0) = 0.0; h_sv(2,1) = 8.2; h_sv(2,2) = -0.4; h_sv(2,3) = 1.8;
    MatrixXd g_sv(6, 6);
    g_sv = MatrixXd::Zero(6, 6);
    MatrixXd h_sv(6, 6);
    h_sv = MatrixXd::Zero(6, 6);

    const double deg2rad = M_PI / 180.0;
    lat = 90 - lat;
    lat = lat*deg2rad;
    lon = lon*deg2rad;
    // radius of earth
    double a = 6371.2;
    double r = a + alt;
    // year since 2015 for secular variation
    int dt = year - 2015;
    // magnetic field components
    double B_r = 0; double B_lat = 0; double B_lon = 0;


    // MatrixXd P_d(3,4);
    // P_d(0,0) = -sin(lat);
    // P_d(0,1) = (-sin(lat)*cos(lat))/(sqrt(1 - pow(cos(lat),2)));
    // P_d(0,2) = 0;
    // P_d(0,3) = 0;
    // P_d(1,0) = -3*sin(lat)*cos(lat);
    // P_d(1,1) = (3*sin(lat)*(1 - 2*pow(cos(lat), 2)))/(sqrt(1 - pow(cos(lat),2)));
    // P_d(1,2) = 6*sin(lat)*cos(lat);
    // P_d(1,3) = 0;
    // P_d(2,0) = (-3/2)*sin(lat)*(5*pow(cos(lat), 2) - 1);
    // P_d(2,1) = (3*sin(lat)*cos(lat)*(11 - 15*pow(cos(lat), 2)))/(2*sqrt(1 - pow(cos(lat),2)));
    // P_d(2,2) = 15*sin(lat)*(3*pow(cos(lat), 2) - 1);
    // P_d(2,3) = -45*sin(lat)*cos(lat)*sqrt(1 - pow(cos(lat), 2));

    MatrixXd g = get_g_coefficients();
    MatrixXd h = get_h_coefficients();
    MatrixXd P = get_P_coefficients(cos(lat));
    MatrixXd Pd = get_Pd_coefficients(P, cos(lat));
    // IGRF model B field calculation, n/m sets order (3)

    for(int n = 1; n < 6; n++){
        for(int m = 0; m < 6; m++){

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

//PYBIND11_MODULE(magnetic_field_cpp, m) {
//    m.doc() = "Magnetic Field Model"; // optional module docstring
//    m.def("get_magnetic_field", &get_magnetic_field, "Get magnetic field (ENU) based on lla and year");
//}

MatrixXd get_g_coefficients(){
    MatrixXd g(6, 6);
    g.row(0) << 0, 0, 0 ,0, 0, 0;
    g(1, 0) = -29442.0; g(1, 1) = -1501.0;  g(1, 2) = 0.0;      g(1, 3) = 0.0;      g(1, 4) = 0.0;      g(1, 5) = 0.0;
    g(2, 0) = -2445.1;  g(2, 1) = 3012.9;   g(2, 2) = 1676.7;   g(2, 3) = 0.0;      g(2, 4) = 0.0;      g(2, 5) = 0.0;
    g(3, 0) = 1350.7;   g(3, 1) = -2352.3;  g(3, 2) = 1225.6;   g(3, 3) = 582.0;    g(3, 4) = 0.0;      g(3, 5) = 0.0;
    g(4, 0) = 907.6;    g(4, 1) = 813.7;    g(4, 2) = 120.4;    g(4, 3) = -334.9;   g(4, 4) = 70.4;     g(4, 5) = 0.0;
    g(5, 0) = -232.6;   g(5, 1) = 360.1;    g(5, 2) = 192.4;    g(5, 3) = -140.9;   g(5, 4) = -157.5;   g(5, 5) = 4.1;
    return g;
}

MatrixXd get_h_coefficients(){
    MatrixXd h(6, 6);
    h.row(0) << 0, 0, 0 ,0, 0, 0;
    h(1, 1) = 4797.1;
    h(2, 1) = -2845.6;  h(2, 2) = -641.9;
    h(3, 1) = -115.3;   h(3, 2) = 244.9;    h(3, 3) = -538.4;
    h(4, 1) = 283.3;    h(4, 2) = -188.7;   h(4, 3) = 180.9;    h(4, 4) = -329.5;
    h(5, 1) = 47.3;     h(5, 2) = 197.0;    h(5, 3) = -119.3;   h(5, 4) = 16.0;     h(5, 5) = 100.2;
    return h;
}

double doub_fact(int x){
    double val = 1.0;
    for(int k = 0; k <= ceil(x/2 - 1); k++){
        val = val * (x - 2 * k);
    }
    return val;
}

MatrixXd get_P_coefficients(double x){

    MatrixXd P(7, 7);
    P.row(0) << 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
    P(1, 1) = sqrt(1 - pow(x, 2));
    for(int n = 1; n < 7; n++){
        for(int m = 0; m < 7; m++){

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

MatrixXd get_Pd_coefficients(const MatrixXd P, double x){

    MatrixXd Pd(6, 6);
    Pd.row(0) << 0, 0, 0 ,0, 0, 0;
    Pd(1,1) = x;
    for(int n = 1; n < 6; n++){
        for(int m = 0; m < 6; m++){
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

    m.def("get_magnetic_field", &get_magnetic_field, "Gives mag field in NED at the given lat lon alt and year(geocentric)");

}










