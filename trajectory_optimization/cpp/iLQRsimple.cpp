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
	for ( int k = 0; k < N-1; k++ ) {
		auto Jk = 0.5*((xtraj(all, k) - xg).transpose()) * Q * (xtraj(all, k) - xg) + 0.5*(utraj(all, k).transpose()) * R * utraj(all, k);
		J += Jk(0);

		// Need to check if can pass an index to a matrix by reference (eg. xnew(all, k+1))
		// Otherwise have to pass the entire matrix and the step 'k'
		rkstep(xtraj(all, k), utraj(all, k), dt, xtraj(all, k+1), A(all, seq(Nx*k, Nx*(k+1)-1)), B(all, seq(Nu*k, Nu*(k+1)-1)));
	}
	// Add terminal cost
	auto Jn = 0.5*((xtraj(all, N) - xg).transpose()) * Qf * ((xtraj(all, N) - xg);
	J += Jn(0);
	Jhist[0] = J;


	// Intialize matrices for optimisation
	MatrixXd S = MatrixXd::Zero(Nx, Nx);
	MatrixXd s = MatrixXd::Zero(Nx, 1);
	MatrixXd l = MatrixXd::Constant(Nu, N, 1 + tol);
	MatrixXd q(Nx, 1);
	MatrixXd r(Nu, 1);

	MatrixXd Snew = S; // temporary matrices for updating
	MatrixXd snew = s;
	MatrixXd Ak(Nx, Nx);
	MatrixXd Bk(Nx, Nu);
	MatrixXd Kk(Nu, Nx);
	MatrixXd LH(Nu, Nu);

	// Line search temp matrices and params
	MatrixXd xnew = MatrixXd::Zero(Nx, N);
	MatrixXd unew = MatrixXd::Zero(Nu, N-1);
	double alpha;
	double Jnew;

	int iter = 0;
	while ( l.lpNorm<Infinity>() > tol ) {

		iter += 1;

		// Initialize backwards pass
		S << Qf;
		s << Qf*(xtraj(all, N) - xg);
		for ( int k = (N-1); k >= 0; k-- ) {

			// Calculate cost gradients (for this time step)
			q = Q * (xtraj(all, k) - xg);
			r = R * utraj(all, k);

			// Calculate feed-forward correction and feedback gain
			// Need to check that these assignments don't slow the algorithm down too much. (Alternative code at botoom of file)
			Ak = A(all, seq(Nx*k, Nx*(k+1)-1));
			Bk = B(all, seq(Nu*k, Nu*(k+1)-1));
			

			LH = (R + Bk.transpose()*S*Bk);
			l(all, k) = LH.colPivHouseholderQr().solve((r + Bk.transpose()*s));
			K(all, seq(Nx*k, Nx*(k+1)-1)) = LH.colPivHouseholderQr().solve(Bk.transpose()*S*Ak);

			// Calculate new cost to go matrices (Sk, sk)
			Kk = K(all, seq(Nx*k, Nx*(k+1)-1));
			Snew = Q + Kk.transpose()*R*Kk + (Ak - Bk*Kk).transpose()*S*(Ak - Bk*Kk);
			snew = q - Kk.transpose()*R + Kk.transpose()*R*l(all, k) + (Ak - Bk*Kk).transpose()*(s - S*Bk*l(all, k));
			S = Snew;
			s = snew;
		}

		// Forward pass line search with new l and K
		xnew(all, 0) = xtraj(all, 0);
		Jnew = J + 1;
		alpha = 1;

		while ( Jnew > J ) {
			Jnew = 0;
			for ( int k = 0; k < N; k++ ) {
				unew(all, k) = utraj(all, k) - alpha*l(all, k) - K(all, seq(Nx*k, Nx*(k+1)-1))*(xnew(all, k) - xtraj(all, k));
				rkstep(xnew(all, k), unew(all, k), dt, xnew(all, k+1), A(all, seq(Nx*k, Nx*(k+1)-1)), B(all, seq(Nu*k, Nu*(k+1)-1)));

				// Find new cost
				auto Jk = 0.5*((xnew(all, k) - xg).transpose())*Q*(xnew(all, k) - xg) + 0.5*(unew(all, k).transpose())*R*unew(all, k);
				J += Jk(0);
			}
			auto Jn = 0.5*((xtraj(all, N) - xg).transpose()) * Qf * ((xtraj(all, N) - xg);
			Jnew += Jn(0);
			alpha *= 0.5;
		}
		xtraj = xnew; // Make sure this assigns to the reference as desired
		utraj = unew;
		J = Jnew;
		Jhist[iter] = J;
	}

	return;
}


/**
  *
  * Performs a single midpoint (2nd order Runge-Kutta) step
  *
  */
void rkstep(const MatrixXd& x0, const MatrixXd& u0, double dt, MatrixXd& x1, MatrixXd& A, MatrixXd& B) {

	// Define sizes (could pass as parameters)
	double Nx = static_cast<unsigned int>( x0.rows() );
	double Nu = static_cast<unsigned int>( u0.rows() );

	// Initialize matrices to pass to dynamics
	// In final iplemenation, the size of these matrices should be specified
	MatrixXd xdot1(Nx, 1), xdot2(Nx, 1);
	MatrixXd dxdot1(Nx, Nx+Nu), dxdot2(Nx, Nx+Nu);

	//Initialize matrices for RK midpoint step
	MatrixXd A1, A2;
	MatrixXd B1, B2;

	// x is reference and is passed to pendulumDynamics by reference
	pendulumDynamics(0, x0, u0, xdot1, dxdot1);
	pendulumDynamics(0, x0 + dt*0.5*xdot1, u0, xdot2, dxdot2);

	x1 = x0 + dt * xdot2;

	A1 = dxdot1(all, seq(0, Nx));
	A2 = dxdot2(all, seq(0, Nx));

	B1 = dxdot1(all, seq(Nx+1, Nx+Nu));
	B2 = dxdot2(all, seq(Nx+1, Nx+Nu));

	A = MatrixXd::Identity(Nx, Nx) + dt*A2 + 0.5*dt*dt*A2*A1;
	B = dt*B2 + 0.5*dt*dt*A2*B1;

	return;
}


/*
			// This code is less readable than assigning Ak, Bk, Kk as new MatrixXd variables
			// but indexing should be more compuationally efficient.
			int bIdxStrt = Nu * k 
			int bIdxEnd = Nu * k + Nu -1;
			int aIdxStrt = Nx * k;
			int aIdxEnd = Nx * k + Nx -1;

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
