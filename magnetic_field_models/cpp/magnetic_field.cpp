//
// Created by ayotundedemuren on 10/28/19.
//

#include "magnetic_field.h"
#include <math.h>
#include <iostream>
#include "../../eigen-git-mirror/Eigen/Dense"
using namespace std;
using namespace Eigen;

VectorXd get_magnetic_field(double lat, double lon, double alt, int year);
MatrixXd get_P_coefficients(double x);
MatrixXd get_g_coefficients();
MatrixXd get_h_coefficients();
MatrixXd get_Pd_coefficients(MatrixXd P, double x);
int main() {
    double lat = 45; double lon = 45; double alt = 400; int year = 2015;
    // VectorXd B = get_magnetic_field(lat, lon, alt, year);
    MatrixXd P = get_P_coefficients(cos(lat * M_PI / 180.0));
    MatrixXd Pd = get_Pd_coefficients(P, cos(M_PI/2 - lat*M_PI/180.0));
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
    MatrixXd g_sv(5, 6);
    g_sv = MatrixXd::Zero(5, 6);
    MatrixXd h_sv(5, 6);
    h_sv = MatrixXd::Zero(5, 6);

    const double deg2rad = M_PI / 180.0;
    lat = 90 - lat;
    lat = lat*deg2rad;
    lon = lon*deg2rad;
    // radius of earth
    double a = 6378.1;
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

    for(int n = 0; n < 5; n++){
        for(int m = 0; m < 6; m++){
            B_r += pow((a/r), n + 2) * (n+1) * ((g(n, m) + dt * g_sv(n, m)) * cos(m * lon) +
                    (h(n, m) + dt * h_sv(n, m)) * sin(m * lon)) * P(n + 1, m);

            B_lat += - pow((a/r), n + 2) * ((g(n, m) + dt * g_sv(n, m)) * sin(m * lon) +
                    m * (h(n, m) + dt * h_sv(n, m)) * cos(m * lon)) * Pd(n + 1,m);

            B_lon += (-1/sin(lat)) * pow((a/r), n + 2) * m * (-(g(n, m) + dt * g_sv(n, m)) * sin(m * lon) + 
                (h(n, m) + dt * h_sv(n, m)) * cos(m * lon)) * P(n + 1, m);

        }
    }
    

    VectorXd B_vec(3);
    // NED coordinate frame
    B_vec(0) = -B_lat;
    B_vec(1) = B_lon;
    B_vec(2) = -B_r;


    return B_vec;
}

MatrixXd get_g_coefficients(){
    MatrixXd g(5, 6);
    g(0, 0) = -29442.0; g(0, 1) = -1501.0;  g(0, 2) = 0.0;      g(0, 3) = 0.0;      g(0, 4) = 0.0;      g(0, 5) = 0.0;
    g(1, 0) = -2445.1;  g(1, 1) = 3012.9;   g(1, 2) = 1676.7;   g(1, 3) = 0.0;      g(1, 4) = 0.0;      g(1, 5) = 0.0;
    g(2, 0) = 1350.7;   g(2, 1) = -2352.3;  g(2, 2) = 1225.6;   g(2, 3) = 582.0;    g(2, 4) = 0.0;      g(2, 5) = 0.0;
    g(3, 0) = 907.6;    g(3, 1) = 813.7;    g(3, 2) = 120.4;    g(3, 3) = -334.9;   g(3, 4) = 70.4;     g(3, 5) = 0.0;
    g(4, 0) - -232.6;   g(4, 1) = 360.1;    g(4, 2) = 192.4;    g(4, 3) = -140.9;   g(4, 4) = -157.5;   g(4, 5) = 4.1;
    return g;
}

MatrixXd get_h_coefficients(){
    MatrixXd h(5, 6);
    h(0, 1) = 4797.1;
    h(1, 1) = -2845.6;  h(1, 2) = -641.9;
    h(2, 1) = -115.3;   h(2, 2) = 244.9;    h(2, 3) = -538.4;
    h(3, 1) = 283.3;    h(3, 2) = -188.7;   h(3, 3) = 180.9;    h(3, 4) = -329.5;
    h(4, 1) = 47.3;     h(4, 2) = 197.0;    h(4, 3) = -119.3;   h(4, 4) = 16.0;     h(4, 5) = 100.2;
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
    for(int n = 1; n < 7; n++){
        for(int m = 0; m < 7; m++){

            if(n == m){
                P(n, m) = pow((-1), n) * doub_fact(2 * n - 1) * pow(1 - pow(x, 2), (n / 2.0));

            }else if(n < m){
                P(n, m) = 0.0;

            }else if(m == 0 && n == 1){
                P(n, m) = x * (2 * m + 1) * P(n - 1, m);

            }else if(n > 1){
                P(n, m) = ((2 * (n-1) + 1) * x * P(n - 1, m) - (n-1 + m) * P(n - 2, m)) / (n - m);

            }else{
                P(n, m) = -(2.0 * (n-1) + 1.0) * sqrt(1 - pow(x, 2)) * P(n-1, m-1);

            }
            
        }
    }
    return P;
}

MatrixXd get_Pd_coefficients(const MatrixXd P, double x){

    MatrixXd Pd(6, 6);
    Pd.row(0) << 0, 0, 0 ,0, 0, 0;
    for(int n = 1; n < 6; n++){
        for(int m = 0; m < 6; m++){
            // cout << (pow(x, 2.0) - 1.0 ) << endl;
            if(n < m){
                Pd(n, m) = 0.0;  
            }else{
                Pd(n, m) = 1.0/(2.0 * n + 1.0) * ((n + 1.0) * (n + m) * P(n - 1, m) - n * (n - m + 1.0) * P(n + 1, m)) / (1.0 - pow(x, 2));
            }
            
                    
        }
    }
    return Pd;
}













