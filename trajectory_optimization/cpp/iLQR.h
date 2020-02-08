/*
 * Author: Nick Goodson
 * Jan 30th 2020
 *
*/


#ifndef ILQR_H
#define ILQR_H

#include <iostream>
#include <cmath>
#include <vector>
#include <sstream>
#include <fstream>
#include "../../eigen-git-mirror/Eigen/Dense"


/* iLQRsimple.cpp */
bool iLQR(const Eigen::MatrixXd& xg,
                const Eigen::MatrixXd& Qw,
                const Eigen::MatrixXd& R,
                const Eigen::MatrixXd& Qwf,
                const Eigen::MatrixXd& Qqf,
                const double dt,
                const double tol,
                Eigen::MatrixXd& xtraj,
                Eigen::MatrixXd& utraj,
                Eigen::MatrixXd& K,
                std::vector<double>& Jhist);

void rkstep(const Eigen::MatrixXd& u0, double dt, int k, Eigen::MatrixXd& x, Eigen::MatrixXd& A, Eigen::MatrixXd& B);
void satelliteDynamics(double t, const Eigen::MatrixXd& x, const Eigen::MatrixXd& u, Eigen::MatrixXd& xdot, Eigen::MatrixXd& dxdot);


/* utils.cpp */
template<typename T>
std::string toString(T const& value);

void writeToFile(const Eigen::MatrixXd& xtraj, const Eigen::MatrixXd& utraj, const std::vector<double>& Jhist);

void printMatrix(const Eigen::MatrixXd& mat);


#endif  // ILQR_H

