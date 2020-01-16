/*
 *
 * Author: Nick Goodson 
 * Jan 15th 2020
 *
*/

#include "iLQR.h"
#include <iostream>
#include <cmath>
#include "../../eigen-git-mirror/Eigen/Dense"

using namespace Eigen;
using namespace std;



int main() {

	// Type definition for pointer to dynamics function (for clarity)
	typedef void (*dynamics)(double, const MatrixXd&, const MatrixXd&, VectorXd&, MatrixXd&);
	dynamics pendDynPtr = &pendulumDynamics;

	// Define Sizes
	int Nx = 2;
	int Nu = 1;
	int N = 249;

	// Cost matrices
	MatrixXd Qf = MatrixXd::Identity(Nx, Nx) * 30;
	MatrixXd Q = MatrixXd::Identity(Nx, Nx) * 0.01;
	double R = 3.0;

	// Initial and final states
	VectorXd x0 = VectorXd::Zero(Nx, 1)

	VectorXd xg;
	xg << M_PI, 0;

	// Outputs from iLQR
	MatrixXd xtraj = MatrixXd::Zero(Nx, N);
	MatrixXd utraj = Matrix::Zero(1, N-1); // Initial control trajectory
	MatrixXd K = MatrixXd::Zero(Nx, Nu*N);



	// Call to iLQR function
	//iLQRsimple(pendDynPtr, )
}



/**
  * Simulates the pendulum's dynamics
  * Could split this into two functions (linearized, non-linear) each returning a single value
  * rather than passing in the outputs by reference
  * 
  @ t, current simulation time
  @ x, state vector
  @ x
  */
void pendulumDynamics(double t, const MatrixXd& x, const MatrixXd& u, VectorXd& xdot, MatrixXd& dxdot) {

	// parameters
	double m = 1.0; // kg
	double l = .5;  // m
	double b = 0.1; // kg m^2 /s
	double lc = 0.5; // m
	double I = 0.25; // (m*l^2)  kg*m^2
	double g = 9.81; // m/s^2

	// Non-linear EOM's
	q = x(0);
	qd = x(1);
	qdd = (u - m * g * lc * sin(q) - b * qd) / I;

	// Returning xdot vector
	xdot(0) = qd;
	xdot(1) = qdd;

	// Returning concatenated matrices of linearized dynamics (jacobians)
	// dxdot = [A, B]
	dxdot(0, 0) = 0;
	dxdot(0, 1) = 1;
	dxdot(0, 2) = 0;
	dxdot(1, 0) = -m * g * lc * cos(q) / I;
	dxdot(1, 1) = -b / I;
	dxdot(1, 2) = 1.0 / I;

}



