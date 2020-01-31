/*
 * Author: Nick Goodson
 * Jan 30th 2020
 *
*/

#include "iLQR.h"

using namespace Eigen;


/**
 * Test the iLQR algorithm on the pendulum
 */
int main() {
    // Define sizes and sim params
    const int Nx = 7;  //(q, w)
    const int Nu = 3;
    const double dt = 0.01;
    const double tol = 0.001;
    int N = 300;

    // Cost matrices
    MatrixXd Qf = MatrixXd::Identity(Nx, Nx) * 100;
    MatrixXd Q = MatrixXd::Identity(Nx, Nx) * 0.01;
    MatrixXd R(Nu, Nu);
    R << 0.3;

    // Initial and final states
    MatrixXd x0 = MatrixXd::Zero(Nx, 1);
    MatrixXd xg(Nx, 1);
    xg << M_PI, 0;

    // Outputs from iLQR (intitialized)
    MatrixXd xtraj = MatrixXd::Zero(Nx, N);
    xtraj(all, 0) = x0;
    MatrixXd utraj = MatrixXd::Zero(Nu, N-1);
    MatrixXd K = MatrixXd::Zero(Nu, Nx*(N-1));
    vector<double> Jhist;  // Size depends on how many iterations of while loop run. Use for testing only, don't implement on MCU

    // Call to iLQR function
    bool success = iLQRsimple(xg, Q, R, Qf, dt, tol, xtraj, utraj, K, Jhist);
    cout << "Results: " << success << endl;

    // Write results to a csv file for comparison with MATLAB or python
    writeToFile(xtraj, utraj, Jhist);

    return 0;
}