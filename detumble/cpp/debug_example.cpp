//sample program for debugging
// Paul DeTrempe
#include<iostream>
#include "../../eigen-git-mirror/Eigen/Dense"

using namespace Eigen;
using namespace std;

Vector3d get_bias_estimate( MatrixXd B_mat);

int main(){
	MatrixXd B_mat_test = MatrixXd::Random(100,3);	
	Vector3d bias_true;
	bias_true << 1, 1, 1;

	std::cout<< bias_true << std::endl;


	// Add in bias
	for (int i = 0; i < B_mat_test.rows(); i++)
	{
		B_mat_test.row(i) = B_mat_test.row(i) + bias_true.transpose();
	}

	std::cout<<B_mat_test<<std::endl;


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

    std::cout<<v<<std::endl;

    Matrix4d A;
    A.row(0) << v(0), v(3), v(4), v(6);
    A.row(1) << v(3), v(1), v(5), v(7);
    A.row(2) << v(4), v(5), v(2), v(8);
    A.row(3) << v(6), v(7), v(8), v(9);

    std::cout<<A<<std::endl;


    Matrix3d A_concat;
    A_concat << A.block(0,0,3,3);


    std::cout<<A_concat<<std::endl;

    bias = -A_concat.colPivHouseholderQr().solve(v.segment(6,3));


    // TODO: implement output for translating system, find evals and evecs about origin to find parameters of T matrix
    //       (scaling and orientation errors)

    return bias;
}
