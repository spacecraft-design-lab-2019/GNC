
#include <iostream>
#include <cmath>
#include <typeinfo>
#include "../../eigen-git-mirror/Eigen/Dense"

using namespace Eigen;
using namespace std;

/* Script for testing passing Matrices between functions, Mapping memory blocks into matrices, etc */

int main() {


	// Testing if output of x^TAx can be converted to double (Yes)
	MatrixXd m(2, 2);	
	m << 1, 2,
		3, 4;

	Vector2d v = Vector2d::Constant(2);
	cout << m << "\n" << v << endl;

	double x = v.transpose()*m*v;
	cout << x << endl;
	cout << typeid(x).name() << endl;


	// Testing passing slices of matrices to functions
	Matrix<double, 3, 5> A = Matrix<double, 3, 5>::Constant(0.2);
	MatrixXd B = MatrixXd::Constant(2, 6);
	MatrixXd C = MatrixXd::Constant(3, 6, 0.5);

	myTestLater1(A, B, C);

	return 0;
}


// Testing if slices of matrics can be passed in as arguments
void myTestLayer1(const MatrixXd& A, MatrixXd& B, MatrixXd C){

	cout << "A = " << A << endl;
	cout << "B = " << B << endl;
	cout << "C = " << C << endl;

	myTestLayer2(A(all, seq(1, 3)), C);

	return 20;
}


void myTestLayer2(const MatrixXd& D, MatrixXd& E) {

	cout << "D = " << D << endl;


	return 30;
}