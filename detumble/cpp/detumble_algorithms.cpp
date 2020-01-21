//
// Created by Paul DeTrempe on 11/5/2019
//
#include "detumble_algorithms.h"
#include "../../eigen-git-mirror/Eigen/Dense"
#include "../../pybind11/include/pybind11/pybind11.h"
#include "../../pybind11/include/pybind11/eigen.h"
#include <cmath>
#include <iostream>
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
Vector3d get_bias_estimate( MatrixXd B_mat);
Vector3d detumble_B_cross_bang_bang(Vector3d omega, Vector3d B, double k, Vector3d max_dipoles);
Vector3d detumble_B_cross_directional(Vector3d omega, Vector3d B, double k, Vector3d max_dipoles);
void detumble_B_dot_C(double* B_dot, double* max_dipoles);

int main(){
    // MatrixXd B_mat_test = MatrixXd::Zero(10,3);
    // std::cout<< B_mat_test << std::endl;

    // Vector3d bias_estimated;
    // bias_estimated = get_bias_estimate(MatrixXd B_mat);

    return 0;
}

// C wrapper functions

extern "C" {

    void detumble_B_dot_C(double* B_dot, double* max_dipoles, double* commanded_dipole){
        // Takes in two pointers to C/C++ arrays and one pointer to where you would like the answer stored
        // Updates the contents of the answer array

  // ----------------------Try to add two vectors using Eigen (map to Eigen, add, map back)-------

        // Create Eigen vectors from C/C++ arrays
        Vector3d B_dot_vec = Map<Vector3d>(B_dot);
        Vector3d max_dipole_vec = Map<Vector3d>(max_dipoles);
        Vector3d commanded_dipole_vec = Map<Vector3d>(commanded_dipole);

        // Perform operation on vector/matrix (using Eigen notation)
        commanded_dipole_vec = detumble_B_dot_bang_bang(B_dot_vec, max_dipole_vec);


        // Copy contents back into array preallocated for holding answer used at higher level by C calling function
        for (int i = 0; i < 3; ++i)
        {
            commanded_dipole[i] = commanded_dipole_vec(i);
        }

    }
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

Vector3d detumble_B_cross_bang_bang(Vector3d omega, Vector3d B, double k, Vector3d max_dipoles){
    /*
    Takes in the angular rate [rad/s] and magnetic field [nanoTesla] vectors (in principal frame)
    (as 3x1 np arrays), and a vector specifying the maximum dipoles in each axis as a 3x1 vector.
    It returns a 3x1 vector specifying a magnetic dipole command.

    - based on Wertz (7.49) and created from my brain.
    TODO: find optimal k for our system
    */

    Vector3d b_hat, M_temp, M;   // unit B field, control moment

    b_hat = B/B.norm();
    M_temp = k*omega.cross(b_hat);

    if( signbit(M_temp(0))){        // signbit returns true is a number is negative
        M(0) = -max_dipoles(0);
    }
    else{
        M(0) = max_dipoles(0);
    }

    if( signbit(M_temp(1)) ){        // signbit returns true is a number is negative
        M(1) = -max_dipoles(1);
    }
    else{
        M(1) = max_dipoles(1);
    }

    if( signbit(M_temp(2)) ){        // signbit returns true is a number is negative
        M(2) = -max_dipoles(2);
    }
    else{
        M(2) = max_dipoles(2);
    }

    return M;
}

Vector3d detumble_B_cross_directional(Vector3d omega, Vector3d B, double k, Vector3d max_dipoles){
    /*
    Takes in the angular rate [rad/s] and magnetic field [nanoTesla] vectors (in principal frame)
    (as 3x1 np arrays), and a vector specifying the maximum dipoles in each axis as a 3x1 vector.
    It returns a 3x1 vector specifying a magnetic dipole command

    - based on Wertz (7.49) and created from my brain.
    TODO: find optimal k for our system
    */

    Vector3d b_hat, M_temp, M;   // unit B field, control moment
    double ratio, max_ratio;

    b_hat = B/B.norm();
    M_temp = k*omega.cross(b_hat);


    max_ratio = 1.0;
    for (int i = 0; i < 3; ++i) // 3 is number of dipole axes
    {
        ratio = abs(M_temp(i))/max_dipoles(i);
        if (ratio>max_ratio)
        {
            max_ratio = ratio;
        }
    }

    M = M_temp/max_ratio;

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

Vector3d get_bias_estimate(MatrixXd B_mat){
    /*
    This function takes in an matrix of magnetometer measurements and returns an estimate of the magnetometer bias.
    This implementation uses a least-squares estimate to fit parameters of an ellipsoid. The bias returned is the center
    of this estimated ellipsoid.

    Inputs:
        B_mat - matrix of magnetic field measurements, m x 3 matrix, [same as desired units as bias]

    Outputs:
        bias - vector of estimated bias, 3x1 vector, [same units as input]

    References:
        Ellipsoid fit, MATLAB stack exchange: https://www.mathworks.com/matlabcentral/fileexchange/24693-ellipsoid-fit
        LEAST SQUARES FITTING OF ELLIPSOID USING ORTHOGONAL DISTANCES: http://dx.doi.org/10.1590/S1982-21702015000200019
            (Note, this function uses the non-iterative method described in section 4 of the above paper).
    */

    // build matrix of ellipsoid parameters (LHS of lls problem)
    MatrixXd D(B_mat.rows(),9);
    VectorXd x(B_mat.rows());
    VectorXd y(B_mat.rows());
    VectorXd z(B_mat.rows());

    x = B_mat.col(0);
    y = B_mat.col(1);
    z = B_mat.col(2);

    D.col(0) = x.cwiseProduct(x) + y.cwiseProduct(y) - 2*z.cwiseProduct(z);
    D.col(1) = x.cwiseProduct(x) + z.cwiseProduct(z) - 2*y.cwiseProduct(y);
    D.col(2) = 2*x.cwiseProduct(y);
    D.col(3) = 2*x.cwiseProduct(z);
    D.col(4) = 2*y.cwiseProduct(z);
    D.col(5) = 2*x;
    D.col(6) = 2*y;
    D.col(7) = 2*z;
    D.col(8) = MatrixXd::Ones(B_mat.rows(),1);

    // build vector of distances from origin (RHS of lls problem)
    VectorXd dist(B_mat.rows(),1);
    dist = x.cwiseProduct(x) + y.cwiseProduct(y) + z.cwiseProduct(z);

    // solve the system using SVD for numerical stability (try SVD and QR)
    VectorXd u(9);
    u = D.bdcSvd(ComputeThinU | ComputeThinV).solve(dist);

    // compute center from solved system
    Vector3d bias;
    VectorXd v(10);
    v(0) = u(0) + u(1) -1;
    v(1) = u(0) - 2*u(1) -1;
    v(2) = u(1) - 2*u(0) -1;
    v.tail<7>() = u.tail<7>();

    // std::cout<<v<<std::endl;

    Matrix4d A;
    A.row(0) << v(0), v(3), v(4), v(6);
    A.row(1) << v(3), v(1), v(5), v(7);
    A.row(2) << v(4), v(5), v(2), v(8);
    A.row(3) << v(6), v(7), v(8), v(9);

    // std::cout<<A<<std::endl;


    Matrix3d A_concat;
    A_concat << A.block(0,0,3,3);


    // std::cout<<A_concat<<std::endl;

    bias = -A_concat.colPivHouseholderQr().solve(v.segment(6,3));


    // TODO: implement output for translating system, find evals and evecs about origin to find parameters of T matrix
    //       (scaling and orientation errors)

    return bias;
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
    m.def("detumble_B_cross_bang_bang", &detumble_B_cross_bang_bang, "Returns moment needed for detumble using bang bang B _cross");
    m.def("detumble_B_cross_directional", &detumble_B_cross_directional, "Returns moment needed for detumble using direction B_cross");
    m.def("detumble_B_dot", &detumble_B_dot, "Returns moment need for detumble using B dot");
    m.def("detumble_B_dot_bang_bang", &detumble_B_dot_bang_bang, "Bang bang controller for detumbling");
    m.def("get_bias_estimate", &get_bias_estimate, "Estimates magnetometer bias based on matrix of measurements");
    m.def("get_B_dot", &get_B_dot,  "Performs simple forward/backward difference to get magnetic field rate");
}