/*
 *
 * Author: Nick Goodson
 * Jan 14th 2020
*/

#include <iostream>
#include <cmath>
#include iLQR.h
#include "../../eigen-git-mirror/Eigen/Dense"

using namespace Eigen
using namespace std


/**
  * Simulates the pendulum's dynamics
  * 
  @ t, current simulation time
  @ x, state vector
  @ x
  */

void pendulumDynamics(double t, const MatrixXd& x, const MatrixXd& u, MatrixXd& x, VectorXd& xdot, MatrixXd& dxdot) {

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

	// Retrning concatenated matrix of linearized dynamics
	// dxdot = [A, B]
	dxdot(0, 0) = 0;
	dxdot(0, 1) = 1;
	dxdot(0, 2) = 0;
	dxdot(1, 0) = -m * g * lc * cos(q) / I;
	dxdot(1, 1) = -b / I;
	dxdot(1, 2) = 1.0 / I;

}

