/*
 * Author: Nick Goodson 
 * Jan 15th 2020
 *
*/

#ifndef GNC_ILQR_H
#define GNC_ILQR_H


#include <iostream>
#include <cmath>
#include <vector>
#include "../../eigen-git-mirror/Eigen/Dense"

// Type definition for pointer to dynamics function (for clarity)
// typedef void (*dynamicsFunc)(double, const MatrixXd&, const MatrixXd&, MatrixXd&, MatrixXd&);	

bool iLQRsimple(Eigen::MatrixXd& x0, 
				Eigen::MatrixXd& xg,  
				Eigen::MatrixXd& Q, 
				Eigen::MatrixXd& R, 
				Eigen::MatrixXd& Qf, 
				double dt, 
				double tol,
				Eigen::MatrixXd& xtraj,
				Eigen::MatrixXd& utraj,
				Eigen::MatrixXd& K,
				std::vector<double>& Jhist);

// void rkstep(const Eigen::MatrixXd& x0, const Eigen::MatrixXd& u0, double dt, Eigen::MatrixXd& x1, Eigen::MatrixXd& A, Eigen::MatrixXd& B);
void rkstep(const Eigen::MatrixXd& u0, double dt, int k, Eigen::MatrixXd& x, Eigen::MatrixXd& A, Eigen::MatrixXd& B);

void pendulumDynamics(double t, const Eigen::MatrixXd& x, const Eigen::MatrixXd& u, Eigen::MatrixXd& xdot, Eigen::MatrixXd& dxdot);


#endif  // GNC_ILQR_H