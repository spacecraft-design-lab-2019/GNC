/*
 * Author: Nick Goodson
 * Jan 30th 2020
 *
*/

#include "iLQR.h"
using namespace Eigen;
using namespace std;


// TODO: potentially should be a variable
#define MAX_ITERS 1000


/**
 * A simple implementation of the Iterative Linear-Quadratic-Regulator algorithm
 *
 * @return, xtraj, utraj, K
 */
bool iLQR(MatrixXd& xg,
                MatrixXd& Q,
                MatrixXd& R,
                MatrixXd& Qf,
                double dt,
                double tol,
                MatrixXd& xtraj,
                MatrixXd& utraj,
                MatrixXd& K,
                vector<double>& Jhist) {


    bool success_flag = true;  // Returns true if algorithm converged

    // Define sizes
    const auto Nx = static_cast<unsigned int>( xtraj.rows() );
    const auto Nu = static_cast<unsigned int>( utraj.rows() );
    const auto N = static_cast<unsigned int>( xtraj.cols() );

    MatrixXd A = MatrixXd::Zero(Nx, Nx * (N-1));
    MatrixXd B = MatrixXd::Zero(Nx, Nu * (N-1));


    // Forward simulate with initial controls utraj0
    double J = 0;
    for ( int k = 0; k < N-1; k++ ) {
        J = J + (0.5*((xtraj(all, k) - xg).transpose()) * Q * (xtraj(all, k) - xg) + 0.5*(utraj(all, k).transpose()) * R * utraj(all, k))(0);
        rkstep(utraj(all, k), dt, k, xtraj, A, B);
    }
    J = J + (0.5*((xtraj(all, N-1) - xg).transpose()) * Qf * ((xtraj(all, N-1) - xg)))(0); 	// Add terminal cost
    Jhist.push_back(J);

    // Intialize matrices for optimisation
    MatrixXd S = MatrixXd::Zero(Nx, Nx);
    MatrixXd s = MatrixXd::Zero(Nx, 1);
    MatrixXd l = MatrixXd::Constant(Nu, N-1, 1 + tol);
    MatrixXd q(Nx, 1);
    MatrixXd r(Nu, 1);

    MatrixXd Snew = S; // temporary matrices for update step
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
    while ( l.lpNorm<Infinity>() > tol) {
        cout << "lp-norm: " << l.lpNorm<Infinity>() << endl;

        iter += 1;
        if (iter > MAX_ITERS){
            // Break from the loop and don't perform the maneuver
            success_flag = false;
            return success_flag;
        }

        // Initialize backwards pass
        S << Qf;
        s << Qf*(xtraj(all, N-1) - xg);
        for ( int k = (int)N-2; k >= 0; k-- ) {

            // Calculate cost gradients (for this time step)
            q = Q * (xtraj(all, k) - xg);
            r = R * utraj(all, k);

            // Calculate feed-forward correction and feedback gain
            Ak = A(all, seq(Nx*k, Nx*(k+1)-1));
            Bk = B(all, seq(Nu*k, Nu*(k+1)-1));

            // Cholesky
            // TODO: Ensure matrix is positive definite before using Cholesky decomposition
            LH = (R + Bk.transpose()*S*Bk);
            l(all, k) = LH.llt().solve((r + Bk.transpose()*s));
            K(all, seq(Nx*k, Nx*(k+1)-1)) = LH.llt().solve(Bk.transpose()*S*Ak);

            // Calculate new cost to go matrices (Sk, sk)
            Kk = K(all, seq(Nx*k, Nx*(k+1)-1));
            Snew = Q + Kk.transpose()*R*Kk + (Ak - Bk*Kk).transpose()*S*(Ak - Bk*Kk);
            snew = q - Kk.transpose()*r + Kk.transpose()*R*l(all, k) + (Ak - Bk*Kk).transpose()*(s - S*Bk*l(all, k));
            S = Snew;
            s = snew;
        }

        // Forward pass line search with new l and K
        xnew(all, 0) = xtraj(all, 0);
        Jnew = J + 1;
        alpha = 1;

        while ( Jnew > J ) {
            Jnew = 0;
            for (int k = 0; k < N - 1; k++) {
                unew(all, k) = utraj(all, k) - alpha * l(all, k) -
                               K(all, seq(Nx * k, Nx * (k + 1) - 1)) * (xnew(all, k) - xtraj(all, k));
                rkstep(unew(all, k), dt, k, xnew, A, B);
                Jnew = Jnew + (0.5 * ((xnew(all, k) - xg).transpose()) * Q * (xnew(all, k) - xg) +
                               0.5 * (unew(all, k).transpose()) * R * unew(all, k))(0);
            }
            Jnew = Jnew + (0.5 * ((xnew(all, N - 1) - xg).transpose()) * Qf * ((xnew(all, N - 1) - xg)))(0);
            alpha *= 0.5;
        }
        xtraj = xnew;
        utraj = unew;
        J = Jnew;
        Jhist.push_back(J);
    }
    return success_flag;
}


/**
  * Performs a single midpoint (2nd order Runge-Kutta) step
  * Can modify this function to include the satellites dynamics
  * in the final implemenation to prevent having to pass matrices to a dynamics function
  *
  */
void rkstep(const MatrixXd& u0, double dt, int k, MatrixXd& xtraj, MatrixXd& A, MatrixXd& B) {

    // Define sizes (hard code in final version)
    auto Nx = static_cast<unsigned int>( xtraj.rows() );
    auto Nu = static_cast<unsigned int>( u0.rows() );

    // Extract current state
    MatrixXd x0 = xtraj(all, k);

    // Initialize matrices
    MatrixXd xdot1(Nx, 1), xdot2(Nx, 1);
    MatrixXd dxdot1(Nx, Nx+Nu), dxdot2(Nx, Nx+Nu);  // [A B]
    MatrixXd A1, A2;
    MatrixXd B1, B2;

    // Step Dynamics (update xtraj, A and B at timestep k+1)
    satelliteDynamics(0, x0, u0, xdot1, dxdot1);
    satelliteDynamics(0, x0 + dt*0.5*xdot1, u0, xdot2, dxdot2);
    xtraj(all, k+1) = x0 + dt * xdot2;

    A1 = dxdot1(all, seq(0, Nx-1));
    A2 = dxdot2(all, seq(0, Nx-1));

    B1 = dxdot1(all, seq(Nx, Nx+Nu-1));
    B2 = dxdot2(all, seq(Nx, Nx+Nu-1));

    A(all, seq(Nx*k, Nx*(k+1)-1)) = MatrixXd::Identity(Nx, Nx) + dt*A2 + 0.5*dt*dt*A2*A1;
    B(all, seq(Nu*k, Nu*(k+1)-1)) = dt*B2 + 0.5*dt*dt*A2*B1;
}


/**
  * Simulates the pendulum's dynamics. Used for forward step with runge-kutta integrator.
  *
  @ t, current simulation time
  @ x, current state vector (Nx, 1)
  @ u, control input (Nu, 1)
  @ xdot, state vector derivative (return value)
  @ dxdot, error-state vector jacobian [A, B] (return value)
  */
void satelliteDynamics(double t, const MatrixXd& x, const MatrixXd& u, MatrixXd& xdot, MatrixXd& dxdot) {

    // parameters (Probably should be passed in as a configuration variable)
    const double m = 1.0; // kg

    // Non-linear EOM's  (Returning xdot vector)
    xdot(0, 0) = x(1, 0);
    xdot(1, 0) = (u(0, 0) - m * g * lc * sin(x(0, 0)) - b * x(1, 0)) / I;

    // Returning concatenated matrices of linearized dynamics (jacobians)
    // dxdot = [A, B]
    dxdot(0, 0) = 0;
    dxdot(0, 1) = 1;
    dxdot(0, 2) = 0;
    dxdot(1, 0) = -m * g * lc * cos(x(0, 0)) / I;
    dxdot(1, 1) = -b / I;
    dxdot(1, 2) = 1.0 / I;

}