/*
 *
 * Author: Nick Goodson 
 * Jan 14th 2020
 *
*/


#include "iLQR.h"
#include <iostream>
#include <cmath>
#include "../../eigen-git-mirror/Eigen/Dense"

using namespace Eigen;
using namespace std;


// Type definition for pointer to dynamics function (for clarity)
typedef void (*dynamics)(double, const MatrixXd&, const MatrixXd&, VectorXd&, MatrixXd&);	


/**
  * Simple iLQR implementation
  *
  *
  */
void iLQRsimple(dynamics pendDynPtr,
				VectorXd& x0, 
				VectorXd& xg, 
				MatrixXd& utraj0, 
				MatrixXd& Q, 
				double R, 
				MatrixXd& Qf, 
				double dt, 
				double tol
				MatrixXd& xtraj
				MatrixXd& utraj
				MatrixXd& K
				MatrixXd& Jhist) {

	// Define sizes
	double Nx = (unsigned int) x0.rows();
	double Nu = (unsigned int) utraj0.rows();
	double N = (unsigned int) utraj0.cols() + 1;





}



