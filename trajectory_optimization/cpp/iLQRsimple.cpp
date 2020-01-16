/*
 *
 * Author: Nick Goodson 
 * Jan 14th 2020
 *
*/

#include "iLQR.h"
using namespace Eigen;
using namespace std;


// Type definition for pointer to dynamics function (for clarity)
//typedef void (*dynamicsFunc)(double, const MatrixXd&, const MatrixXd&, MatrixXd&, MatrixXd&);	


/**
  * Simple iLQR implementation
  *
  *
  */
void iLQRsimple(dynamicsFunc pendDynPtr,
				MatrixXd& x0, 
				MatrixXd& xg,  
				MatrixXd& Q, 
				double R, 
				MatrixXd& Qf, 
				double dt, 
				double tol,
				MatrixXd& xtraj,
				MatrixXd& utraj,
				MatrixXd& K,
				vector<double>& Jhist) {

	// Define sizes
	double Nx = static_cast<unsigned int>( x0.rows() );
	double Nu = static_cast<unsigned int>( utraj.rows() );
	double N = static_cast<unsigned int>( xtraj.cols() );

	MatrixXd A = MatrixXd::Zero(Nx, Nx * (N-1));
	MatrixXd B = MatrixXd::Zero(Nx, Nu * (N-1));

	// Forward simulate with initial controls utraj0
	double J = 0;
	for (int k = 0; k < N-1; k++) {
		auto Jk = 0.5*((xtraj(all, k) - xg).transpose()) * Q * (xtraj(all, k) - xg) + 0.5*(utraj(all, k).transpose()) * R * utraj(all, k);
		J += Jk(0);
		// Perform rk step
	}
	// Add terminal cost
	Jn = 0.5*((xtraj(all, N) - xg).transpose()) * Qf * ((xtraj(all, N) - xg);
	J += Jn(0);
	Jhist[0] = J;


	// Intialize matrices for optimisation
	MatrixXd S = MatrixXd::Zero(Nx, Nx);
	MatrixXd s = MatrixXd::Zero(Nx, 1);
	MatrixXd l = MatrixXd::Constant(Nu, N, 1 + tol);

	MatrixXd Snew = S; // temp matrices for updating
	MatrixXd snew = s;
	MatrixXd q(Nx, 1);
	MatrixXd r(Nu, 1);
	MatrixXd Ak(Nx, Nx);
	MatrixXd Bk(Nx, Nu);
	MatrixXd Kk(Nu, Nx);
	MatrixXd LH;
	

	iter = 0;
	while (l.lpNorm<Infinity>() > tol) {

		iter += 1;

		// Initialize backwards pass
		S << Qf;
		s << Qf*(xtraj(all, N) - xg);
		for (int k = (N-1); k > 0; k--) {

			// Calculate cost gradients (for this time step)
			q = Q * (xtraj(all, k) - xg);
			r = R * utraj(all, k);

			// Calculate feed-forward correction and feedback gain
			// Need to check that these assignments don't slow the algorithm down too much. (Alternative code at botoom of file)
			Ak = A(all, seq(Nx*(k-1), Nx*k));
			Bk = B(all, seq(Nu*(k-1), Nu*k));
			Kk = K(all, seq(Nx*(k-1), Nx*k));

			LH = (R + Bk.transpose()*S*Bk);
			l(all, k) = LH.colPivHouseholderQr().solve((r + Bk.transpose()*s));
			K(all, seq(aIdxStrt, aIdxEnd)) = LH.colPivHouseholderQr().solve(Bk.transpose()*S*Ak);

			// Calculate new cost to go matrices (Sk, sk)
			Snew = Q + Kk.transpose()*R*Kk + (Ak - Bk*Kk).transpose()*S*(Ak - Bk*Kk);
			snew = q - Kk.transpose()*R + Kk.transpose()*R*l(all, k) + (Ak - Bk*Kk).transpose()*(s - S*Bk*l(all, k));
			S = Snew;
			s = snew;
			
		}
	}




}



void rkstep() {


}

/*
			// This code is less readable than assigning Ak, Bk, Kk as new MatrixXd variables
			// but indexing should be more compuationally efficient.
			int bIdxStrt = Nu * k - Nu;
			int bIdxEnd = Nu * k;
			int aIdxStrt = Nx * k - Nx;
			int aIdxEnd = Nx * k;

			MatrixXd LH = (R + B(all, seq(bIdxStrt, bIdxEnd)).transpose()*S*B(all, seq(bIdxStrt, bIdxEnd)));
			MatrixXd l_RH = (r + B(all, seq(bIdxStrt, bIdxEnd)).transpose()*s);
			MatrixXd K_RH = (B(all, seq(bIdxStrt, bIdxEnd)).transpose()*S*A(all, seq(aIdxStrt, aIdxEnd)));

			l(all, k) = LH.colPivHouseholderQr().solve(l_RH);  // Use solver to perform inversion
			K(all, seq(aIdxStrt, aIdxEnd)) = LH.colPivHouseholderQr().solve(K_RH); // K has same columns as A matrix

			// Calculate new cost to go matrices (Sk, sk)
			Snew = Q + K(all, seq(aIdxStrt, aIdxEnd)).transpose()*R*K(all, seq(aIdxStrt, aIdxEnd)) + 
					(A(all, seq(aIdxStrt, aIdxEnd)) - B(all, seq(bIdxStrt, bIdxEnd))*K(all, seq(aIdxStrt, aIdxEnd).transpose()) *
					S*((A(all, seq(aIdxStrt, aIdxEnd)) - B(all, seq(bIdxStrt, bIdxEnd))*K(all, seq(aIdxStrt, aIdxEnd));

*/
