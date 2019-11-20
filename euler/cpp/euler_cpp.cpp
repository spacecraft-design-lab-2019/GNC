//
// Created by Ethan on 10/23/2019.
//

#include "euler_cpp.h"
#include "../../eigen-git-mirror/Eigen/Dense"
#include <../../eigen-git-mirror/Eigen/Geometry> 
#include <../../pybind11/include/pybind11/pybind11.h>
#include <../../pybind11/include/pybind11/eigen.h>
namespace py = pybind11;
using namespace Eigen;
using namespace std;


int main(){
    return 0;
}

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

VectorXd get_q_dot(Vector4d q, Vector3d w){
    /*
     Takes in a quaternion and the rotation rate vector,
    and returns the time derivative of the quaternion 
    Inputs:
        q -  4x1, scalar first, normalized
        w -  3x1, rad/s
    Outputs:
        q_dot - 4x1, scalar first, 1/sec
    */

	Vector3d q_dot;
	MatrixXd What(4,4);
	// What.row(0) << 
	q_dot = 1/2 * What * q;


    return q_dot;

}
PYBIND11_MODULE(euler_cpp, m) {
    m.doc() = "Euler equations and Quat Functions"; // optional module docstring

    m.def("get_w_dot", &get_w_dot, "Gets angular velocity using Euler equations");
    m.def("get_q_dot", &get_q_dot, "Gets derivative of quaterion");
    
}