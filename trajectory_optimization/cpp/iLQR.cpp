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
 * Iterative Linear-Quadratic-Regulator algorithm for a satellite
 * @return, xtraj, utraj, K
 */
bool iLQR(const MatrixXd& xg,
                const MatrixXd& Qw,
                const MatrixXd& R,
                const MatrixXd& Qwf,
                const double Qqf,
                const double dt,
                const double tol,
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
        J = J + min(1 + (xg(seq(0, 3), all).transpose() * xtraj(seq(0, 3), k))(0), 1 - (xg(seq(0, 3), all).transpose() * xtraj(seq(0, 3), k))(0)) +
                (0.5 * (xtraj(seq(4, 7), k) - xg(seq(4, 7), all)).transpose() * Qw * (xtraj(seq(4, 7), k) - xg(seq(4, 7), all)) +
                0.5 * utraj(all, k).transpose() * R * utraj(all, k))(0);
    }
    J = J + Qqf * min(1 + (xg(seq(0, 3), all).transpose() * xtraj(seq(0, 3), N-1))(0), 1 - (xg(seq(0, 3), all).transpose() * xtraj(seq(0, 3), N-1))(0)) +
            (0.5 * (xtraj(seq(4, 7), N-1) - xg(seq(4, 7), all)).transpose() * Qwf * (xtraj(seq(4, 7), N-1) - xg(seq(4, 7), all)))(0); 	// Add terminal cost
    Jhist.push_back(J);

    // Initialize matrices for optimization
    MatrixXd S = MatrixXd::Zero(Nx, Nx);
    MatrixXd s = MatrixXd::Zero(Nx, 1);
    MatrixXd l = MatrixXd::Constant(Nu, N-1, 1 + tol);
    MatrixXd q(Nx, 1);
    MatrixXd r(Nu, 1);

    MatrixXd Snew(Nx, Nx); // temporary matrices for update step
    MatrixXd snew(Nx, 1);
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
                Jnew = Jnew + min(1 + (xg(seq(0, 3), all).transpose() * xnew(seq(0, 3), k))(0), 1 - (xg(seq(0, 3), all).transpose() * xnew(seq(0, 3), k))(0)) +
                    (0.5 * (xnew(seq(4, 7), k) - xg(seq(4, 7), all)).transpose() * Qw * (xnew(seq(4, 7), k) - xg(seq(4, 7), all)) +
                     0.5 * unew(all, k).transpose() * R * unew(all, k))(0);
            }
            Jnew = Jnew + Qqf * min(1 + (xg(seq(0, 3), all).transpose() * xnew(seq(0, 3), N-1))(0), 1 - (xg(seq(0, 3), all).transpose() * xnew(seq(0, 3), N-1))(0)) +
                (0.5 * (xnew(seq(4, 7), N-1) - xg(seq(4, 7), all)).transpose() * Qwf * (xnew(seq(4, 7), N-1) - xg(seq(4, 7), all)))(0); 	// Add terminal cost
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

    // parameters TODO: (Probably should be passed in as a configuration variable)
    MatrixXd J = MatrixXd::Identity(3, 3) * 0.01;  // kgm^2

    // Non-linear EOM's  (Returning xdot vector)
    

    // Returning concatenated matrices of linearized dynamics (jacobians)
    // dxdot = [A, B]
    dxdot(0, 0) = 0;
    dxdot(0, 1) = 1;
    dxdot(0, 2) = 0;
    dxdot(1, 0) = 0;
    dxdot(1, 1) = 0;
    dxdot(1, 2) = 0;
}