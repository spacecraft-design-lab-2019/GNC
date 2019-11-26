//
// Created by Andrew Gatherer on 11/25/19.
//

/* The purpose of this function is to perform an iLQR forwards
 * and backwards pass. For the CircuitPython implementation,
 * this should be run iteratively over so many times until
 * it coan be assumed that convergence has been reached
 *
 * The simulation does forward and then backward passes
 *
 * For info: https://homes.cs.washington.edu/~todorov/papers/LiICINCO04.pdf
 *
 * OUTLINE:
 * You know all the satellite dynamics, as well as the magnetic field
 * At the start, you have your initial condition and you propagate
 * forward according to the discrete dynamics you have established.
 * Then, you have to work backwards and generate the changes to the
 * given input that improve the ultimate solution at each point.
 *
 *
 *
 * INPUTS:
 *  - all inputs U_bar
 *  - Initial condition, final desired condition, dt, knot points
 *
 * OUTPUTS:
 *  - dU, the change in input you should have
 *  - some metric of accuracy (not sure what yet)
 *
 *
 */

#include <iostream>
#include <math.h>
#include <../eigen-git-mirror/Eigen/Dense>
//#include <../../pybind11/include/pybind11/pybind11.h>
//#include <../../pybind11/include/pybind11/eigen.h>
using namespace Eigen;
using namespace std;

VectorXd rk4(double t, double h,  VectorXd y, MatrixXd I, Vector3d B, Vector3d m);
VectorXd state_dot(VectorXd y, MatrixXd I, Vector3d B, Vector3d m);
Matrix3d skew(Vector3d x);

int main(){

    //initialize the MOI matrix
    MatrixXd I(3,3);
    I << 1, 0, 0,
         0, 1, 0,
         0, 0, 1;

    //initialize the time
    int N = 10; //create 100 simulation steps
    double dt = 0.1; //10 Hz
    double t0 = 0;
    VectorXd t(N);

    //initialize the control matrix
    Array<Vector3d, Dynamic, 1> U(N);    

    //initialize the statef
    VectorXd y0(7);
    Array<VectorXd, Dynamic, 1> y(N);
    y0 << 1, 0, 0, 1, 0, 0, 0;
    cout << "y0 = " << endl << y0 << endl; // print it just for shits
    y(0) = y0;

    //initialize the mag field
    Array<Vector3d, Dynamic, 1> B(N);

    //forward pass
    for (int i = 0; i < N-1; i++){
        //cout << i << endl;
        y(i+1) = rk4(t0, dt, y(i), I, 
            B(i), U(i));
        cout << endl << y(i+1) << endl;
    }

    return 0;
}

VectorXd rk4(double t, double h,  VectorXd y, MatrixXd I, Vector3d B, Vector3d m) {
    /* simple rk4 stepping given a time and a position
     *
     * INPUTS:
     * t - current time
     * h - time step
     * y - state vector
     * I - MOI matrix
     * B - magnetic field vector (TODO: change so that we can vary this based on time)
     * m - magnetic moment
     *
     * OUTPUTS:
     * y_new - new state vector at future time step t+h
     */

    VectorXd y_new(7);

    // calculate first step
    VectorXd k1(7);
    k1 = state_dot(y, I, B, m);

    // calculate second step
    VectorXd k2(7);
    k2 = state_dot(y + h*.5*k1, I, B, m);

    // calculate third step
    VectorXd k3(7);
    k3 = state_dot(y + h*.5*k2, I, B, m);

    // calculate forth step
    VectorXd k4(7);
    k4 = state_dot(y + h*k3, I, B, m);

    // put em together
    y_new = y + h*(1.0/6*k1 + 1.0/3*k2 + 1.0/3*k3 +1.0/6*k4);

    // renormalize quaternion
    //cout << y_new.segment<4>(3).norm() << endl;
    y_new.segment<4>(3) /= y_new.segment<4>(3).norm();

    return y_new;
}

VectorXd state_dot(VectorXd y, MatrixXd I, Vector3d B, Vector3d m){
    /*this just calculates the dot product of the state
     * INPUTS:
     * y - (omega, q)
     * I - MOI matrix
     * B - magnetic field vector
     * m - the magnetic moment input
     *
     * OUTPUTS:
     * y_dot - (omega_dot, q_dot)
     */

    VectorXd y_dot(7);

    Vector3d omega(0,0,0);
    omega = y.segment<3>(0); // specify rotation rate
    VectorXd q(4);
    q = y.segment<4>(3); // specify the quaternion

    // omega dot
    Vector3d c(0,0,0); //temporary cross product thingy
    c = I*omega;
    y_dot.segment<3>(0) = I.inverse() * (-skew(omega)*c + skew(m)*B); // generate omega_dot

    // q dot
    Matrix4d omega_expanded = Matrix4d::Zero();
    omega_expanded.block<1,3>(0,1) = -omega.transpose();
    omega_expanded.block<3,1>(1,0) = omega;
    omega_expanded.block<3,3>(1,1) = skew(omega); //generated R(omega)
    y_dot.segment<4>(3) = .5 * omega_expanded * q; // generate q_dot

    //cout << "y_dot is: " << y_dot << endl;
    return y_dot;
}

Matrix3d skew(Vector3d x){
    /* this function delivers the skew-symmetric matrix for a given vector
     * INPUTS:
     * x - 3 value vector
     * OUTPUTS:
     * x_hat - 3x3 matrix
     */

    Matrix3d x_hat = Matrix3d::Zero(); //initialize

    x_hat(0,1) = -x(2);
    x_hat(0,2) = x(1);
    x_hat(1,0) = x(2);
    x_hat(1,2) = -x(0);
    x_hat(2,0) = -x(1);
    x_hat(2,1) = x(0);

    return x_hat;
}
