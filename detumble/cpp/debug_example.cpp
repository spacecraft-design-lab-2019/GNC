//sample program for debugging
// Paul DeTrempe
#include<iostream>
#include "../../eigen-git-mirror/Eigen/Dense"

using namespace Eigen;
using namespace std;

Vector3d get_bias_estimate( MatrixXd B_mat);

int main(){
    MatrixXd B_mat_test = MatrixXd::Random(10,3);
    // std::cout<< B_mat_test << std::endl;

    Vector3d bias_estimated;
    bias_estimated = get_bias_estimate(B_mat_test);
    std::cout<<bias_estimated <<std::endl;

    return 0;
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

    // build matrix of ellipsoid parameters
    MatrixXd D(B_mat.rows(),9);
    VectorXd x(B_mat.rows());
    VectorXd y(B_mat.rows());
    VectorXd z(B_mat.rows());

    x = B_mat.col(0);
    y = B_mat.col(1);
    z = B_mat.col(2);

    D.col(0) = x.cwiseProduct(x) + y.cwiseProduct(y) - 2*z;
    D.col(1) = x.cwiseProduct(x) + z.cwiseProduct(z) - 2*x;
    D.col(2) = 2*x.cwiseProduct(y);
    D.col(3) = 2*x.cwiseProduct(z);
    D.col(4) = 2*y.cwiseProduct(z);
    D.col(5) = 2*x;
    D.col(6) = 2*y;
    D.col(7) = 2*z;
    D.col(8) = MatrixXd::Ones(B_mat.rows(),1);

    // build vector of distances from origin


    // solve the system using SVD for numerical stability (try SVD and QR)

    // compute center from solved system
    Vector3d bias;
    bias << 1, 2, 3;


    // TODO: implement output for translating system, find evals and evecs about origin to find parameters of T matrix
    //       (scaling and orientation errors)

    return bias;
}
