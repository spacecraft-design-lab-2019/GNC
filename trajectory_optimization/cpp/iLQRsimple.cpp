/*
 *
 * Author: Nick Goodson 
 * Jan 14th 2020
 *
*/


#include "iLQR.h"
using namespace Eigen;
using namespace std;


// Type definition for pointer to dynamics function (for clarity)
//typedef void (*dynamicsFunc)(double, const MatrixXd&, const MatrixXd&, MatrixXd&, MatrixXd&);	


/**
  * Simple iLQR implementation
  *
  *
  */
void iLQRsimple(dynamicsFunc pendDynPtr,
				MatrixXd& x0, 
				MatrixXd& xg,  
				MatrixXd& Q, 
				double R, 
				MatrixXd& Qf, 
				double dt, 
				double tol,
				MatrixXd& xtraj,
				MatrixXd& utraj,
				MatrixXd& K,
				vector<double>& Jhist) {

	// Define sizes
	double Nx = (unsigned int) x0.rows();
	double Nu = (unsigned int) utraj.rows();
	double N = (unsigned int) utraj.cols() + 1;

	MatrixXd A = MatrixXd::Zero(Nx, Nx * N-1);
	MatrixXd B = MatrixXd::Zero(Nx, Nu * N-1);

	// Forward simulate with initial controls utraj0
	double J = 0;
	for (int k = 0; k < N-1; k++) {
		auto Jk = 0.5*((xtraj(all, k) - xg).transpose()) * Q * (xtraj(all, k) - xg) + 0.5*(utraj(all, k).transpose()) * R * utraj(all, k);
		J += Jk(0);
		// Perform rk step
	}
	// Add terminal cost
	Jn = 0.5*((xtraj(all, N) - xg).transpose()) * Qf * ((xtraj(all, N) - xg);
	J += Jn(0);
	Jhist[0] = J;

	MatrixXd S = MatrixXd::Zero(Nx, Nx);
	MatrixXd s = MatrixXd::Zero(Nx, 1);
	MatrixXd l = MatrixXd::Zero(Nu, N); 

}

