//
// Created by Ethan on 10/23/2019.
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
Vector3d get_B_dot(Vector3d B1, Vector3d B2, dt);

int main(){
    return 0;
}

Vector3d detumble_B_cross(Vector3d omega, Vector3d B, double k){

}

Vector3d detumble_B_dot(B,B_dot,k=1.0){

}

Vector3d get_B_dot(B1,B2,dt){

}

PYBIND11_MODULE(frame_conversions_cpp, m) {
    m.doc() = "Frame Conversions"; // optional module docstring

    m.def("detumble_B_cross", &detumble_B_cross, "Returns moment needed for detumble using B cross");
    m.def("detumble_B_dot", &detumble_B_dot, "Returns moment need for detumble using B dot");
    m.def("get_B_dot", &get_B_dot  "Performs simple forward/backward difference to get magnetic field rate");
}