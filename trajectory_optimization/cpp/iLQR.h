

#ifndef GNC_ILQR_H
#define GNC_ILQR_H


#include <iostream>
#include <cmath>
#include <vector>
#include "../../eigen-git-mirror/Eigen/Dense"
using namespace Eigen;
using namespace std;

// Type definition for pointer to dynamics function (for clarity)
typedef void (*dynamicsFunc)(double, const MatrixXd&, const MatrixXd&, MatrixXd&, MatrixXd&);	

void pendulumDynamics(double t, const MatrixXd& x, const MatrixXd& u, MatrixXd& xdot, MatrixXd& dxdot);

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
				vector<double>& Jhist);

void rkstep();

#endif  // GNC_ILQR_H