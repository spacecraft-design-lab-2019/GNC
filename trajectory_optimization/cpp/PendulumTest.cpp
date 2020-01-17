/*
 *
 * Author: Nick Goodson 
 * Jan 15th 2020
 *
*/

#include "iLQR.h"
using namespace Eigen;
using namespace std;


int main() {

	// Type definition for pointer to dynamics function
	//typedef void (*dynamics)(double, const MatrixXd&, const MatrixXd&, MatrixXd&, MatrixXd&);
	dynamicsFunc pendDynPtr = &pendulumDynamics;

	// Define sizes and sim params
	int Nx = 2;
	int Nu = 1;
	int N = 249;
	double dt = 0.01;
	double tol = 0.001;

	// Cost matrices
	MatrixXd Qf = MatrixXd::Identity(Nx, Nx) * 30;
	MatrixXd Q = MatrixXd::Identity(Nx, Nx) * 0.01;
	double R = 3.0;

	// Initial and final states
	MatrixXd x0 = MatrixXd::Zero(Nx, 1);

	MatrixXd xg(Nx, 1);
	xg << M_PI, 0;

	// Outputs from iLQR (intiialized)
	MatrixXd xtraj = MatrixXd::Zero(Nx, N); 
	MatrixXd utraj = MatrixXd::Zero(Nu, N-1);  // Initial control trajectory
	MatrixXd K = MatrixXd::Zero(Nu, Nx*(N-1));
	vector<double> Jhist;  // Size depends on how many iterations of while loop run. Use for testing only, don't implement on MCU

	// Call to iLQR function
	iLQRsimple(pendDynPtr, x0, xg, Q, R, Qf, dt, tol, xtraj, utraj, K, Jhist);
}



/**
  * Simulates the pendulum's dynamics. Used for forward step with runge-kutta integrator.
  * 
  @ t, current simulation time
  @ x, state vector
  @ u, control input
  @ xdot, state vector derivative (return value)
  @ dxdot, error-state vector jacobian [A, B] (return value)
  */
void pendulumDynamics(double t, const MatrixXd& x, const MatrixXd& u, MatrixXd& xdot, MatrixXd& dxdot) {

	// parameters
	double m = 1.0; // kg
	double l = .5;  // m
	double b = 0.1; // kg m^2 /s
	double lc = 0.5; // m
	double I = 0.25; // (m*l^2)  kg*m^2
	double g = 9.81; // m/s^2

	// Non-linear EOM's
	double q = x(0, 0);
	double qd = x(1, 0);
	double qdd = (u(0, 0) - m * g * lc * sin(q) - b * qd) / I;

	// Returning xdot vector
	xdot(0, 0) = qd;
	xdot(1, 0) = qdd;

	// Returning concatenated matrices of linearized dynamics (jacobians)
	// dxdot = [A, B]
	dxdot(0, 0) = 0;
	dxdot(0, 1) = 1;
	dxdot(0, 2) = 0;
	dxdot(1, 0) = -m * g * lc * cos(q) / I;
	dxdot(1, 1) = -b / I;
	dxdot(1, 2) = 1.0 / I;

}



