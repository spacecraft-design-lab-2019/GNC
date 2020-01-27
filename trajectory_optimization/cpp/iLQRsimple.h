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
#include <sstream>
#include <fstream>
#include "../../eigen-git-mirror/Eigen/Dense"


/* iLQRsimple.cpp */
bool iLQRsimple(Eigen::MatrixXd& xg,  
				Eigen::MatrixXd& Q, 
				Eigen::MatrixXd& R, 
				Eigen::MatrixXd& Qf, 
				double dt, 
				double tol,
				Eigen::MatrixXd& xtraj,
				Eigen::MatrixXd& utraj,
				Eigen::MatrixXd& K,
				std::vector<double>& Jhist);

void rkstep(const Eigen::MatrixXd& u0, double dt, int k, Eigen::MatrixXd& x, Eigen::MatrixXd& A, Eigen::MatrixXd& B);


/* PendulumTest.cpp */
void pendulumDynamics(double t, const Eigen::MatrixXd& x, const Eigen::MatrixXd& u, Eigen::MatrixXd& xdot, Eigen::MatrixXd& dxdot);


/* utils.cpp */
template<typename T>
std::string toString(T const& value);

void writeToFile(const Eigen::MatrixXd& xtraj, const Eigen::MatrixXd& utraj, const std::vector<double>& Jhist);

void printMatrix(const Eigen::MatrixXd& mat);


#endif  // GNC_ILQR_H