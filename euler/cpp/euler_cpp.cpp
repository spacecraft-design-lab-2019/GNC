//
// Created by Ethan on 10/23/2019.
//

#include "euler_cpp.h"
#include "../../eigen-git-mirror/Eigen/Dense"
#include <../../eigen-git-mirror/Eigen/Geometry> 
#include <../../pybind11/include/pybind11/pybind11.h>
#include <../../pybind11/include/pybind11/eigen.h>
#include <iostream>

namespace py = pybind11;
using namespace Eigen;
using namespace std;


int main(){
    return 0;
}

MatrixXd hat(Vector3d w);
MatrixXd Lq(Vector4d q);
MatrixXd Rq(Vector4d q);

Vector3d get_w_dot(Vector3d w, Vector3d M, MatrixXd I){
    /*
    Takes in angular rate, net torque (3x1, principal frame), and the
    principal moment of inertia matrix (3x3, principal frame). Returns the rates of change
    of the angular rates as a 3x1.

    Implements the Euler equations.
    Inputs:
        w - angular velocity vector, 3x1, principal frame, rad/s
        M - vector of moments, 3x1, principal frame, N-m
        I - principal moment of inertia matrix, 3x3, kg-m^2
    Outputs:
        w_dot - angular acceleration vector, 3x1, principal frame, rad/s^2
   */

	Vector3d w_dot;
	Vector3d Iw = I * w;
	w_dot = I.inverse() * (M - w.cross(Iw));


    return w_dot;

}

Vector4d get_q_dot(Vector4d q, Vector3d w){
    /*
     Takes in a quaternion and the rotation rate vector,
    and returns the time derivative of the quaternion 
    Inputs:
        q -  4x1, scalar first, normalized
        w -  3x1, rad/s
    Outputs:
        q_dot - 4x1, scalar first, 1/sec
    */

	Vector4d q_dot;
	Vector4d wvec;

	wvec << 0, w(0), w(1), w(2);
	q_dot = .5 * Lq(q) * wvec;

    return q_dot;

}

VectorXd get_attitude_derivative(double t, VectorXd x, Vector3d M, MatrixXd I){
    /*
    Takes in an attitude state parametrized by a quaternion and an angular rate and returns a state derivative.

    Inputs:
        t - time, scalar, [sec]
        x - state, [q;w], 7x1 vector of stacked q and w.
            q - quaternion, 4x1 vector, scalar first, represents rotation from body coordinates to ECI coordinates
            w - angular rate, 3x1 vector, [rad/s], expressed in body frame of the spacecraft
        M - Torques on spacecraft, 3x1 vector, [Newton-meters], expressed in body frame of the spacecraft
        I - Moment of Inertia matrix, 3x3 matrix, [kg-m^2]

    Outputs:
        x_dot - rate of change of state, [q_dot;w_dot]
    */

    VectorXd x_dot(7,1);
    Vector4d q_dot;
    Vector3d w_dot;
    Vector4d q = x(seq(0,3));
    Vector3d w = x(seq(4,6),0);

    q_dot = get_q_dot(q, w);
    w_dot = get_w_dot(w, M, I);
    x_dot << q_dot, w_dot;

    return x_dot;
}

MatrixXd Lq(Vector4d q){
    /*
    TODO: add documentation
    */
	MatrixXd Lq(4, 4);
	double s = q(0);
	Vector3d v;
	v << q(1), q(2), q(3);
	Lq.row(0) << q(0), -q(1), -q(2), -q(3);
	Lq.col(0).tail(3) << q(1), q(2), q(3);
	Lq.block(1, 1, 3, 3) << s * MatrixXd::Identity(3, 3) + hat(v);

	return Lq;
}

MatrixXd Rq(Vector4d q){
    /*
    TODO: add documentation
    */
	MatrixXd Rq(4, 4);
	double s = q(0);
	Vector3d v;
	v << q(1), q(2), q(3);
	Rq.row(0) << q(0), -q(1), -q(2), -q(3);
	Rq.col(0).tail(3) << q(1), q(2), q(3);
	Rq.block(1, 1, 3, 3) << s * MatrixXd::Identity(3, 3) - hat(v);

	return Rq;
}

MatrixXd hat(Vector3d w){
    /*
    This function takes in a vector, w,  and outputs a skew-symmetric matrix, w_hat, that represents
    the vector cross product of that vector. Premultiplying a vector by this matrix is the same as taking the cross product
    of that vector with w ( cross(w,vector) ).

    Inputs:
        w -3x1 vector. (Often an angular velocity [rad/s]).
    Outputs:
        w_hat - 3x3 matrix
    */
	MatrixXd w_hat(3,3);
	w_hat.row(0) << 0, -w(2), w(1);
	w_hat.row(1) << w(2), 0, -w(0);
	w_hat.row(2) << -w(1), w(0), 0;


	return w_hat;

}

Vector3d rotate_vec(Vector3d xb, Vector4d q){
	/*
	Goes from BODY to INERTIAL IFF q represents the BODY to INERTIAL attitude
	*/
	Vector4d xtemp_sol, xtemp;
	xtemp << 0, xb;
	MatrixXd lq = Lq(q);
	MatrixXd rq = Rq(q);
	xtemp_sol = lq * rq.transpose() * xtemp;
	return xtemp_sol.tail(3);
}

Vector4d get_inverse_quaternion(Vector4d q){
    /*
    This function takes in a quaternion and outputs its inverse (i.e. if a quaternion describes a rotation from body to ECI,
    this function returns the quaternion describing rotation from ECI to body coordinates)

    Inputs
        q - quaternion, 4x1, scalar first
    Outputs
        q_inv - quaternion, 4x1, scalar first
    */
    Vector4d q_inv;
    q_inv << -q(0), q(1), q(2), q(3);
    return q_inv;

}

PYBIND11_MODULE(euler_cpp, m) {
    m.doc() = "Euler equations and Quat Functions"; // optional module docstring

    m.def("get_w_dot", &get_w_dot, "Gets angular velocity using Euler equations");
    m.def("get_q_dot", &get_q_dot, "Gets derivative of quaterion");
    m.def("Lq", &Lq, "Left quat multiply");
    m.def("Rq", &Rq, "Right quat multiply");
    m.def("hat", &hat, "Gets hat matrix");
    m.def("rotate_vec", &rotate_vec, "Rotates vector x from body to inertial");
    m.def("get_inverse_quaternion", &get_inverse_quaternion, "Returns inverse quaternion rotation");
    m.def("get_attitude_derivative", &get_attitude_derivative, "Returns derivative of attitude state");
    
}