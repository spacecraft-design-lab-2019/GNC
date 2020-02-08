/*
 * Author: Nick Goodson
 * Jan 30th 2020
 *
*/

#include "iLQR.h"

using namespace Eigen;
using namespace std;


/**
 * Test the iLQR algorithm on a desired attitude change
 */
int main() {
    // Define sizes and sim params
    const int Nx = 7;  //(q, w)
    const int Nw = 3; // dimension of rotational dynamics
    const int Nu = 3;
    const double dt = 0.01;
    const double tol = 0.001;
    int N = 300;

    // Cost matrices
    MatrixXd Qwf = MatrixXd::Identity(Nw, Nw) * 100;  // Omega terminal cost
    MatrixXd Qw = MatrixXd::Identity(Nw, Nw) * 0.01;  // Omega cumulative cost
    MatrixXd R = MatrixXd::Identity(Nu, Nu) * 0.3;  // Control cumulative cost
    const double Qqf = 30;  // Attitude (quaternion) terminal cost

    // Initial and final states
    double theta0 = M_PI / 2;
    double thetaF = 3 * M_PI / 2;
    MatrixXd x0(Nx, 1);
    MatrixXd xg(Nx, 1);
    x0 << sin(theta0/2), 0, 0, cos(theta0/2), 0, 0, 0;
    xg << sin(thetaF/2), 0, 0, cos(thetaF/2), 0, 0, 0;

    // Outputs from iLQR (intitialized)
    MatrixXd xtraj = MatrixXd::Zero(Nx, N);
    xtraj(all, 0) = x0;
    MatrixXd utraj = MatrixXd::Zero(Nu, N-1);
    MatrixXd K = MatrixXd::Zero(Nu, Nx*(N-1));
    vector<double> Jhist;  // Size depends on how many iterations of while loop run. Use for testing only, don't implement on MCU

    // Call to iLQR function
    bool success = iLQRsimple(xg, Qw, R, Qwf, Qqf, dt, tol, xtraj, utraj, K, Jhist);
    cout << "Results: " << success << endl;

    // Write results to a csv file for comparison with MATLAB or python
    writeToFile(xtraj, utraj, Jhist);

    return 0;
}