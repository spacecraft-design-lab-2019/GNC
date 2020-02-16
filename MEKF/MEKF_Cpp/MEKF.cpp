#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "MEKF.hpp"

using namespace std;
using Eigen::MatrixXd;
//using Eigen::

MatrixXd readMatrix(const char *filename);
double trace(MatrixXd A);

#define MAXBUFSIZE  ((int) 1e6)

int main() {
    // initialize W and V based on MEKF template data
    MatrixXd W(6,6);
    MatrixXd V(6,6);
    for (int i = 0; i<6; i++){
        for (int j = 0; j<6; j++){
            if (i==j){
                if (i < 3){
                    W(i,j) = 0.003;
                    V(i,j) = 0.003;
                }
                else{
                    W(i,j) = 0.0076;
                    V(i,j) = 0.0076;
                }
            }

        }
    }


    // Measurement History
    MatrixXd rB1hist(3,1501);
    MatrixXd rB2hist(3,1501);
    MatrixXd whist(3,1501);
    rB1hist = readMatrix("rB1hist.txt");
    rB2hist = readMatrix("rB2hist.txt");
    whist = readMatrix("whist.txt");

    // two measurements for TRIAD
    MatrixXd rN1(3,1);
    MatrixXd rN2(3,1);
    rN1(0,0) = 0.6693;
    rN1(1,0) = -0.6818;
    rN1(2,0) = 0.2952;
    rN2(0,0) = 0.2627;
    rN2(1,0) = -0.7082;
    rN2(2,0) = -0.6553;
    // vert stack rN
    MatrixXd rN(6,1);
    rN(seq(0,2),0) = rN1;
    rN(seq(3,5),0) = rN2;

    double dt = 0.1;

    // initialize covariance matrix
    MatrixXd Pk(6,6);
    Pk = pow((10.0*3.14159/180.0),2.0)*MatrixXd::Identity(6,6);
    MatrixXd rN_full(3,2);
    MatrixXd rB_full(3,2);
    for (int i = 0; i<3;i++){
        rN_full(i,0) = rN1(i,0);
        rB_full(i,0) = rB1hist(i,0);
        rN_full(i,1) = rN2(i,0);
        rB_full(i,1) = rB2hist(i,0);
    }
    // initialize quaternion using TRIAD
    MatrixXd R0(3,3);
    R0 = triad_ad(rN_full, rB_full);
    MatrixXd q0(4,1);
    q0 = DCM2q(R0);

    // gyro bias
    MatrixXd Beta0(3,1);
    Beta0(0,0) = 0.1;
    Beta0(1,0) = 0.1;
    Beta0(2,0) = 0.1;

    MatrixXd xk(7,1);
    xk(seq(0,3),0) = q0;
    xk(seq(4,6),0) = Beta0;

    // pre-allocate
    MatrixXd q_mekf(4,1501), b_mekf(3,1501);
    q_mekf(all,0) = q0;
    b_mekf(all,0) = Beta0;
    // start sim
   for (int i = 0; i < 1500; i++){
        /*
         * xn - predicted state (mu_k+1|k)
         * Pn - predicted covariance (sigma_k+1|k)
         * A - linearized state transition matrix
         * W - process noise covariance
         * V - measurement noise covariance
         * rN - vector measurement in newtonian frame
         * rB - vector measurement in body frame
         */
        // get measurement vector
        MatrixXd rB1(3,1), rB2(3,1), rB(6,1);
        for (int j = 0; j < 3; j++){
            rB1(j,0) = rB1hist(j,i+1);
            rB2(j,0) = rB2hist(j,i+1);
            rB(j,0) = rB1(j,0);
            rB(j+3,0) = rB2(j,0);
        }

        // predict step
        MatrixXd xn(7,1), A(6,6);
        predict(xk,whist(all,i+1),dt,xn,A);
        MatrixXd Pn(6,6);
        Pn = A*Pk*A.transpose()+W;

        // run measurement step
        MatrixXd y(6,1), R(3,3), C(6,6);
        measurement(xn(seq(0,3),0),rN,y,R,C);

        // run innovation step to find z
        MatrixXd z(6,1), S(6,6);
        Innovation(R, rN, rB, Pn, V, C, z, S);

        // find Kalman Gain
        MatrixXd L(6,6);
        L = Pn * C.transpose() * S.inverse();

        // update step
        update(L,z,xn,Pn,V,C,xk,Pk);

        // update state
        q_mekf(all,i+1) = xk(seq(0,3),0);
        b_mekf(all,i+1) = xk(seq(4,6),0);
    }
    ofstream myfile;
    myfile.open("output_C.txt");
    myfile << q_mekf.transpose() << endl;
    myfile.close();
    return 0;
}

// function to read in measurement data
MatrixXd readMatrix(const char *filename)
{
    int cols = 0, rows = 0;
    double buff[MAXBUFSIZE];

    // Read numbers from file into buffer.
    ifstream infile;
    infile.open(filename);
    while (! infile.eof())
    {
        string line;
        getline(infile, line);

        int temp_cols = 0;
        stringstream stream(line);
        while(! stream.eof())
            stream >> buff[cols*rows+temp_cols++];

        if (temp_cols == 0)
            continue;

        if (cols == 0)
            cols = temp_cols;

        rows++;
    }

    infile.close();

    rows--;

    // Populate matrix with numbers.
    MatrixXd result(rows,cols);
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            result(i,j) = buff[ cols*i+j ];

    return result;
}

double trace(MatrixXd A){
    double a = A(0,0) + A(1,1) + A(2,2);
    return a;
}


