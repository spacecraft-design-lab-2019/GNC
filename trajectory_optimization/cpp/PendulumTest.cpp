/*
 * Author: Nick Goodson 
 * Jan 15th 2020
 *
*/

#include "iLQRsimple.h"
using namespace Eigen;
using namespace std;


/**
  * Simulates the pendulum's dynamics. Used for forward step with runge-kutta integrator.
  * 
  @ t, current simulation time
  @ x, current state vector (Nx, 1)
  @ u, control input (Nu, 1)
  @ xdot, state vector derivative (return value)
  @ dxdot, error-state vector jacobian [A, B] (return value)
  */
void pendulumDynamics(double t, const MatrixXd& x, const MatrixXd& u, MatrixXd& xdot, MatrixXd& dxdot) {

	// parameters
	const double m = 1.0; // kg
	const double l = .5;  // m
	const double b = 0.1; // kg m^2 /s
	const double lc = 0.5; // m
	const double I = 0.25; // (m*l^2)  kg*m^2
	const double g = 9.81; // m/s^2

	// Non-linear EOM's  (Returning xdot vector)
	xdot(0, 0) = x(1, 0);
	xdot(1, 0) = (u(0, 0) - m * g * lc * sin(x(0, 0)) - b * x(1, 0)) / I;

	// Returning concatenated matrices of linearized dynamics (jacobians)
	// dxdot = [A, B]
	dxdot(0, 0) = 0;
	dxdot(0, 1) = 1;
	dxdot(0, 2) = 0;
	dxdot(1, 0) = -m * g * lc * cos(x(0, 0)) / I;
	dxdot(1, 1) = -b / I;
	dxdot(1, 2) = 1.0 / I;

}



// Test the iLQR algorithm on the pendulum
int main() {

	// Define sizes and sim params
	const int Nx = 2;
	const int Nu = 1;
	const double dt = 0.01;
	const double tol = 0.001;
	int N = 249;

	// Cost matrices
	MatrixXd Qf = MatrixXd::Identity(Nx, Nx) * 30;
	MatrixXd Q = MatrixXd::Identity(Nx, Nx) * 0.01;
	MatrixXd R(Nu, Nu);
	R << 3.0;

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
	iLQRsimple(xg, Q, R, Qf, dt, tol, xtraj, utraj, K, Jhist);

	return 0;
}



