//
// Created by Paul DeTrempe on 11/5/2019
//

#include "frame_conversions.h"
#include "../../eigen-git-mirror/Eigen/Dense"
#include <../../pybind11/include/pybind11/pybind11.h>
#include <../../pybind11/include/pybind11/eigen.h>
#include <cmath>
namespace py = pybind11;
using namespace Eigen;
using namespace std;

// constants for default dipole moments
const double dipole_x = 8.8e-3;
const double dipole_y = 1.372e-2;
const double dipole_z = 8.2e-3;

Vector3d detumble_B_cross(Vector3d omega, Vector3d B, double k);
Vector3d detumble_B_dot(Vector3d B, Vector3d B_dot, double k);
Vector3d get_B_dot(Vector3d B1, Vector3d B2, double dt);
Vector3d get_B_dot_bang_bang( Vector3d B_dot, Vector3d max_dipoles( double dipole_x, double dipole_y, double dipole_z) );


int main(){
    return 0;
}

Vector3d detumble_B_cross(Vector3d omega, Vector3d B, double k){
    /*
    Takes in the angular rate [rad/s] and magnetic field [nanoTesla] vectors (in principal frame)
    (as 3x1 np arrays) and returns a 3x1 vector of control torque [N-m] to detumble.

    - based on the paper by AVANZINI AND GIULIETTI
    TODO: find optimal k for our system
    */

    Vector3d b_hat, M;   // unit B field, control moment

    b_hat = B/B.norm();
    M = -k*(MatrixXd::Identity(3, 3) - b_hat*b_hat.transpose())*omega;

    return M;
}

Vector3d detumble_B_dot(Vector3d B, Vector3d B_dot, double k){
    /*
    !-----------This function has been deprecated in favor of detumble_B_dot_bang_bang--------!

    Takes in magnetic field [nanoTesla], magnetic field rate [nanoTesla/s] (as 3x1 vectors, in principal frame), and control gain (scalar)
    and returns a 3x1 control moment [Amp-meters^2]. Torque [N-m] = M [Amp-meters^2] cross B [Tesla]

    TODO: find optimal k for our system
    */

    Vector3d M;

    M = -k*(B_dot/B.norm());

    return M;

}

Vector3d detumble_B_dot_bang_bang(Vector3d B_dot, Vector3d max_dipoles){
    /*
    !-----------------This function replaces detumble_B_dot (11/17/2019)-----------------------!

    Takes in the rate of change of magnetic field as a 3x1 vector, [nanoTesla/s](though it should be noted that units don't
    actually matter for this function, as it uses the sign of the dot product of the magnetic field rate of change with the
    principal axes of the spacecraft, not the magnitude), and an optional input for the maximum dipole as a 3x1 vector,
    [Amp-m^2].

    It returns a dipole acting to maximally oppose the rate of change of the magnetic field.
    */

//    Python code:
//    m = np.zeros((3,1))
//    m[0] = max_dipoles[0] * np.sign(np.dot(np.transpose(np.array([[1.0],[0.0],[0.0]])), np.transpose(B_dot)))
//    m[1] = max_dipoles[1] * np.sign(np.dot(np.transpose(np.array([[0.0],[1.0],[0.0]])) , np.transpose(B_dot)))
//    m[2] = max_dipoles[2] * np.sign(np.dot(np.transpose(np.array([[0.0],[0.0],[1.0]])) , np.transpose(B_dot)))
//    return m

    Vector3d M;             // dipole moment [Amp-m^2]

    if( signbit(B_dot(0)) ){        // signbit returns true is a number is negative
        M(0) = max_dipoles(0);
    }
    else{
        M(0) = -max_dipoles(0);
    }

    if( signbit(B_dot(1)) ){        // signbit returns true is a number is negative
        M(1) = max_dipoles(1);
    }
    else{
        M(1) = -max_dipoles(1);
    }

    if( signbit(B_dot(2)) ){        // signbit returns true is a number is negative
        M(2) = max_dipoles(2);
    }
    else{
        M(2) = -max_dipoles(2);
    }


    return M;

}

Vector3d get_B_dot(Vector3d B1, Vector3d B2, double dt){
    /*
    Takes in two magnetic field measurements (3x1) and the timestep between them and returns
    a first order approximation to the rate of change of the magnetic field (3x1)
    */
//    Reference Python code:
//    B_dot = (B2-B1)/dt
//    return B_dot
    Vector3d B_dot;

    B_dot = (B2-B1)/dt;

    return B_dot;
}


PYBIND11_MODULE(detumble_cpp, m) {
    m.doc() = "Detumble algorithms"; // optional module docstring

    m.def("detumble_B_cross", &detumble_B_cross, "Returns moment needed for detumble using B cross");
    m.def("detumble_B_dot", &detumble_B_dot, "Returns moment need for detumble using B dot");
    m.def("detumble_B_dot_bang_bang", &detumble_B_dot_bang_bang, "Bang bang controller for detumbling");
//    m.def("detumble_B_dot_bang_bang", &detumble_B_dot_bang_bang, "Bang bang controller for detumbling",
//    py::arg("B_dot") = (0.0, 0.0, 0.0), py::arg("max_dipoles") = ( dipole_x, dipole_y, dipole_z));
    m.def("get_B_dot", &get_B_dot,  "Performs simple forward/backward difference to get magnetic field rate");
}