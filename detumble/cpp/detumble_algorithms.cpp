//
// Created by Paul DeTrempe on 11/5/2019
//

#include "frame_conversions.h"
#include "../../eigen-git-mirror/Eigen/Dense"
#include <../../pybind11/include/pybind11/pybind11.h>
#include <../../pybind11/include/pybind11/eigen.h>
namespace py = pybind11;
using namespace Eigen;
using namespace std;

Vector3d detumble_B_cross(Vector3d omega, Vector3d B, double k);
Vector3d detumble_B_dot(Vector3d B, Vector3d B_dot, double k);
Vector3d get_B_dot(Vector3d B1, Vector3d B2, double dt);

int main(){
    return 0;
}

Vector3d detumble_B_cross(Vector3d omega, Vector3d B, double k){
    /*
    Takes in the angular rate and magnetic field vectors (in principal frame)
    (as 3x1 np arrays) and returns a 3x1 vector of control torque to detumble.

    - based on the paper by AVANZINI AND GIULIETTI
    TODO: find optimal k for our system
    */

//    Reference Python code:
//    b = B/np.linalg.norm(B) # normalize magnetic field vector
//    L = -k*np.matmul((np.identity(3) - np.matmul(b,np.transpose(b))),omega)
//    return L
    Vector3d b_hat, M;   // unit B field, control moment

    b_hat = B/B.norm();
    M = -k*(MatrixXd::Identity(3, 3) - b_hat*b_hat.transpose());

    return M;
}

Vector3d detumble_B_dot(Vector3d B, Vector3d B_dot, double k){
    /*
    Takes in magnetic field, magnetic field rate (as 3x1 vectors, in principal frame), and control gain (scalar)
    and returns a 3x1 control moment

    TODO: find optimal k for our system
    */
//    Reference Python code:
//    m = -k*B_dot/np.linalg.norm(B)
//    return m

    Vector3d M;

    M = -k*(B_dot/B.norm());

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
    m.def("get_B_dot", &get_B_dot,  "Performs simple forward/backward difference to get magnetic field rate");
}