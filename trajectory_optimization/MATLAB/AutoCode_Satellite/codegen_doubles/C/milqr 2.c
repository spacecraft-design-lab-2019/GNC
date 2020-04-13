/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: milqr.c
 *
 * MATLAB Coder version            : 4.3
 * C/C++ source code generated on  : 19-Mar-2020 19:51:09
 */

/* Include Files */
#include "milqr.h"
#include <math.h>
#include <string.h>

/* Variable Definitions */
static const signed char iv[9] = { 1, 0, 0, 0, 1, 0, 0, 0, 1 };

static const double dv[9] = { 0.01, 0.0, 0.0, 0.0, 0.01, 0.0, 0.0, 0.0, 0.01 };

static const double dv1[9] = { 0.05, 0.0, 0.0, 0.0, 0.05, 0.0, 0.0, 0.0, 0.05 };

/* Function Declarations */
static void backwardPass(const double fx[17964], const double fu[8982], const
  double cx[3000], const double cu[1497], const double cxx[18000], const double
  cuu[4491], double lambda, const double u_lims[6], const double u[1497], double
  l[1497], double K[8982], double dV[2], boolean_T *diverge);
static void boxQPsolve(const double Quu[9], const double Qu[3], const double
  lower[3], const double upper[3], const double u0[3], double u[3], double
  *result, double Luu_data[], int Luu_size[2], boolean_T b_free[3]);
static void chol_free(const double A_data[], const int A_size[2], double L_data[],
                      int L_size[2], double *fail);
static void forwardRollout(const double x[3500], const double xg[7], const
  double u[1497], const double l[1497], const double K[8982], double alpha,
  const double u_lims[6], double xnew[3500], double unew[1497], double fx[17964],
  double fu[8982], double cx[3000], double cu[1497], double cxx[18000], double
  cuu[4491], double *cost);
static void satellite_cost(const double x[7], const double xg[7], const double
  u[3], double terminal, double *cost, double cx[6], double cu[3], double cxx[36]);

/* Function Definitions */

/*
 * Arguments    : const double fx[17964]
 *                const double fu[8982]
 *                const double cx[3000]
 *                const double cu[1497]
 *                const double cxx[18000]
 *                const double cuu[4491]
 *                double lambda
 *                const double u_lims[6]
 *                const double u[1497]
 *                double l[1497]
 *                double K[8982]
 *                double dV[2]
 *                boolean_T *diverge
 * Return Type  : void
 */
static void backwardPass(const double fx[17964], const double fu[8982], const
  double cx[3000], const double cu[1497], const double cxx[18000], const double
  cuu[4491], double lambda, const double u_lims[6], const double u[1497], double
  l[1497], double K[8982], double dV[2], boolean_T *diverge)
{
  int i;
  double Vx[6];
  int k;
  int i1;
  boolean_T exitg1;
  int b_k;
  double Vxx[36];
  double coeffs;
  int i2;
  double Qu[3];
  double dV_idx_0;
  int i3;
  int u_lims_tmp;
  double Quu[9];
  double b_Quu[9];
  double b_u_lims[3];
  double b_fu[18];
  double c_u_lims[3];
  double Qux[18];
  double b_dv[3];
  int b_u_lims_tmp;
  int c_u_lims_tmp;
  double lk[3];
  double result;
  double Luu_data[9];
  int Luu_size[2];
  boolean_T b_free[3];
  double Kk[18];
  boolean_T y;
  boolean_T exitg2;
  signed char tmp_data[3];
  int n;
  double b_cx[6];
  double Vx_tmp[6];
  double b_cxx[36];
  double b_fx[36];
  int Y_size_idx_0;
  double b_Vx_tmp[18];
  int loop_ub;
  double Y_data[18];
  int coeffs_size_idx_0;
  double c_cxx[36];
  double coeffs_data[3];
  double d;
  int b_i;
  int X_size_idx_0;
  int j;
  int LT_size_idx_0;
  int c_i;
  double LT_data[9];
  signed char b_tmp_data[3];

  /*  Perfoms the LQR backward pass to find the optimal controls */
  /*  Solves a quadratic program (QP) at each timestep for the optimal */
  /*  controls given the control limits */
  /*  Initialize matrices (for C) */
  memset(&l[0], 0, 1497U * sizeof(double));
  memset(&K[0], 0, 8982U * sizeof(double));

  /*  Change in cost */
  dV[0] = 0.0;
  dV[1] = 0.0;

  /*  Set cost-to-go Jacobain and Hessian equal to final costs */
  for (i = 0; i < 6; i++) {
    Vx[i] = cx[i + 2994];
    for (i1 = 0; i1 < 6; i1++) {
      b_k = i1 + 6 * i;
      Vxx[b_k] = cxx[b_k + 17964];
    }
  }

  *diverge = false;
  k = 0;
  exitg1 = false;
  while ((!exitg1) && (k < 499)) {
    /*  Define cost gradients */
    i = 18 * (498 - k);
    for (i1 = 0; i1 < 3; i1++) {
      coeffs = 0.0;
      for (i2 = 0; i2 < 6; i2++) {
        coeffs += fu[(i2 + 6 * i1) + i] * Vx[i2];
        dV_idx_0 = 0.0;
        for (i3 = 0; i3 < 6; i3++) {
          dV_idx_0 += fu[(i3 + 6 * i1) + i] * Vxx[i3 + 6 * i2];
        }

        b_fu[i1 + 3 * i2] = dV_idx_0;
      }

      Qu[i1] = cu[i1 + 3 * (498 - k)] + coeffs;
      for (i2 = 0; i2 < 3; i2++) {
        coeffs = 0.0;
        for (i3 = 0; i3 < 6; i3++) {
          coeffs += b_fu[i1 + 3 * i3] * fu[(i3 + 6 * i2) + 18 * (498 - k)];
        }

        b_k = i1 + 3 * i2;
        b_Quu[b_k] = cuu[b_k + 9 * (498 - k)] + coeffs;
      }
    }

    for (i = 0; i < 3; i++) {
      for (i1 = 0; i1 < 6; i1++) {
        coeffs = 0.0;
        for (i2 = 0; i2 < 6; i2++) {
          coeffs += fu[(i2 + 6 * i) + 18 * (498 - k)] * Vxx[i2 + 6 * i1];
        }

        b_fu[i + 3 * i1] = coeffs;
      }

      for (i1 = 0; i1 < 6; i1++) {
        coeffs = 0.0;
        for (i2 = 0; i2 < 6; i2++) {
          coeffs += b_fu[i + 3 * i2] * fx[(i2 + 6 * i1) + 36 * (498 - k)];
        }

        Qux[i + 3 * i1] = coeffs;
      }
    }

    /*  Regularization (for Cholesky positive definiteness) */
    /*  Solve the Quadratic program with control limits */
    for (i = 0; i < 9; i++) {
      Quu[i] = b_Quu[i] + (double)iv[i] * lambda;
    }

    u_lims_tmp = 3 * (498 - k);
    b_u_lims[0] = u_lims[0] - u[u_lims_tmp];
    c_u_lims[0] = u_lims[3] - u[u_lims_tmp];
    i = 3 * ((int)fmin(499.0, (-(double)k + 499.0) + 1.0) - 1);
    b_dv[0] = -l[i];
    b_u_lims_tmp = u_lims_tmp + 1;
    b_u_lims[1] = u_lims[1] - u[b_u_lims_tmp];
    c_u_lims[1] = u_lims[4] - u[b_u_lims_tmp];
    b_dv[1] = -l[i + 1];
    c_u_lims_tmp = u_lims_tmp + 2;
    b_u_lims[2] = u_lims[2] - u[c_u_lims_tmp];
    c_u_lims[2] = u_lims[5] - u[c_u_lims_tmp];
    b_dv[2] = -l[i + 2];
    boxQPsolve(Quu, Qu, b_u_lims, c_u_lims, b_dv, lk, &result, Luu_data,
               Luu_size, b_free);
    if (result < 1.0) {
      *diverge = true;
      exitg1 = true;
    } else {
      /*  Solve for feedback gains in non-clamped rows of u */
      /*  (using cholesky factor of Quu) */
      memset(&Kk[0], 0, 18U * sizeof(double));
      y = false;
      b_k = 0;
      exitg2 = false;
      while ((!exitg2) && (b_k < 3)) {
        if (b_free[b_k]) {
          y = true;
          exitg2 = true;
        } else {
          b_k++;
        }
      }

      if (y) {
        b_k = 0;
        if (b_free[0]) {
          tmp_data[0] = 1;
          b_k = 1;
        }

        if (b_free[1]) {
          tmp_data[b_k] = 2;
          b_k++;
        }

        if (b_free[2]) {
          tmp_data[b_k] = 3;
        }

        /*  Solves the linear system AX = B for the unknown X using the Cholesky */
        /*  decomposition L of the matrix A. Where LL' = A */
        /*  X can be a vector or a matrix of size n x m */
        /*  Solution is found in O(nm) time using back-substitution */
        /*  This implementation only works for lower triangular factorisations */
        n = Luu_size[0];

        /*  Solve LY = B for Y */
        /* ======================= */
        Y_size_idx_0 = Luu_size[0];
        loop_ub = Luu_size[0] * 6;
        if (0 <= loop_ub - 1) {
          memset(&Y_data[0], 0, loop_ub * sizeof(double));
        }

        coeffs_size_idx_0 = Luu_size[0];
        loop_ub = Luu_size[0];
        if (0 <= loop_ub - 1) {
          memset(&coeffs_data[0], 0, loop_ub * sizeof(double));
        }

        i = Luu_size[0];
        for (b_i = 0; b_i < i; b_i++) {
          if (b_i + 1 != 1) {
            for (i1 = 0; i1 < b_i; i1++) {
              coeffs_data[i1] = Luu_data[b_i + Luu_size[0] * i1];
            }
          }

          for (j = 0; j < 6; j++) {
            if (((signed char)coeffs_size_idx_0 == 1) || (Y_size_idx_0 == 1)) {
              coeffs = 0.0;
              for (i1 = 0; i1 < coeffs_size_idx_0; i1++) {
                coeffs += coeffs_data[i1] * Y_data[i1 + Y_size_idx_0 * j];
              }
            } else {
              coeffs = 0.0;
              for (i1 = 0; i1 < coeffs_size_idx_0; i1++) {
                coeffs += coeffs_data[i1] * Y_data[i1 + Y_size_idx_0 * j];
              }
            }

            Y_data[b_i + Y_size_idx_0 * j] = (Qux[(tmp_data[b_i] + 3 * j) - 1] -
              coeffs) / Luu_data[b_i + Luu_size[0] * b_i];
          }
        }

        /*  Solve L'X = Y for X */
        /* ======================= */
        X_size_idx_0 = Luu_size[0];
        loop_ub = Luu_size[0] * 6;
        if (0 <= loop_ub - 1) {
          memset(&b_fu[0], 0, loop_ub * sizeof(double));
        }

        LT_size_idx_0 = Luu_size[1];
        loop_ub = Luu_size[0];
        coeffs_size_idx_0 = Luu_size[0];
        for (i = 0; i < loop_ub; i++) {
          b_k = Luu_size[1];
          for (i1 = 0; i1 < b_k; i1++) {
            LT_data[i1 + LT_size_idx_0 * i] = Luu_data[i + Luu_size[0] * i1];
          }

          coeffs_data[i] = 0.0;
        }

        i = (int)(((-1.0 - (double)Luu_size[0]) + 1.0) / -1.0);
        for (b_i = 0; b_i < i; b_i++) {
          c_i = (n - b_i) - 1;
          if (c_i + 1 != n) {
            if (c_i + 2 > n) {
              i1 = 0;
              i2 = 0;
              i3 = 0;
            } else {
              i1 = c_i + 1;
              i2 = n;
              i3 = c_i + 1;
            }

            loop_ub = i2 - i1;
            for (i2 = 0; i2 < loop_ub; i2++) {
              coeffs_data[i3 + i2] = LT_data[c_i + LT_size_idx_0 * (i1 + i2)];
            }
          }

          for (j = 0; j < 6; j++) {
            if (((signed char)coeffs_size_idx_0 == 1) || (X_size_idx_0 == 1)) {
              coeffs = 0.0;
              for (i1 = 0; i1 < coeffs_size_idx_0; i1++) {
                coeffs += coeffs_data[i1] * b_fu[i1 + X_size_idx_0 * j];
              }
            } else {
              coeffs = 0.0;
              for (i1 = 0; i1 < coeffs_size_idx_0; i1++) {
                coeffs += coeffs_data[i1] * b_fu[i1 + X_size_idx_0 * j];
              }
            }

            b_fu[c_i + X_size_idx_0 * j] = (Y_data[c_i + Y_size_idx_0 * j] -
              coeffs) / LT_data[c_i + LT_size_idx_0 * c_i];
          }
        }

        b_k = 0;
        if (b_free[0]) {
          b_tmp_data[0] = 1;
          b_k = 1;
        }

        if (b_free[1]) {
          b_tmp_data[b_k] = 2;
          b_k++;
        }

        if (b_free[2]) {
          b_tmp_data[b_k] = 3;
        }

        for (i = 0; i < 6; i++) {
          for (i1 = 0; i1 < X_size_idx_0; i1++) {
            Kk[(b_tmp_data[i1] + 3 * i) - 1] = -b_fu[i1 + X_size_idx_0 * i];
          }
        }
      }

      /*  Update Cost to Go */
      coeffs = 0.0;
      for (i = 0; i < 3; i++) {
        coeffs += ((0.5 * lk[0] * b_Quu[3 * i] + 0.5 * lk[1] * b_Quu[3 * i + 1])
                   + 0.5 * lk[2] * b_Quu[3 * i + 2]) * lk[i];
      }

      dV_idx_0 = dV[0] + ((lk[0] * Qu[0] + lk[1] * Qu[1]) + lk[2] * Qu[2]);
      coeffs += dV[1];
      dV[0] = dV_idx_0;
      dV[1] = coeffs;
      for (i = 0; i < 6; i++) {
        coeffs = Kk[3 * i];
        i1 = 3 * i + 1;
        i2 = 3 * i + 2;
        for (i3 = 0; i3 < 3; i3++) {
          b_Vx_tmp[i + 6 * i3] = (coeffs * b_Quu[3 * i3] + Kk[i1] * b_Quu[3 * i3
            + 1]) + Kk[i2] * b_Quu[3 * i3 + 2];
        }

        dV_idx_0 = 0.0;
        for (i3 = 0; i3 < 6; i3++) {
          dV_idx_0 += fx[(i3 + 6 * i) + 36 * (498 - k)] * Vx[i3];
        }

        b_cx[i] = ((cx[i + 6 * (498 - k)] + dV_idx_0) + ((b_Vx_tmp[i] * lk[0] +
          b_Vx_tmp[i + 6] * lk[1]) + b_Vx_tmp[i + 12] * lk[2])) + ((coeffs * Qu
          [0] + Kk[i1] * Qu[1]) + Kk[i2] * Qu[2]);
        Vx_tmp[i] = (Qux[3 * i] * lk[0] + Qux[i1] * lk[1]) + Qux[i2] * lk[2];
      }

      for (i = 0; i < 6; i++) {
        Vx[i] = b_cx[i] + Vx_tmp[i];
        for (i1 = 0; i1 < 6; i1++) {
          coeffs = 0.0;
          for (i2 = 0; i2 < 6; i2++) {
            coeffs += fx[(i2 + 6 * i) + 36 * (498 - k)] * Vxx[i2 + 6 * i1];
          }

          b_fx[i + 6 * i1] = coeffs;
        }

        for (i1 = 0; i1 < 6; i1++) {
          coeffs = 0.0;
          for (i2 = 0; i2 < 6; i2++) {
            coeffs += b_fx[i + 6 * i2] * fx[(i2 + 6 * i1) + 36 * (498 - k)];
          }

          b_k = i + 6 * i1;
          c_cxx[b_k] = cxx[b_k + 36 * (498 - k)] + coeffs;
        }

        coeffs = b_Vx_tmp[i + 6];
        dV_idx_0 = b_Vx_tmp[i + 12];
        i1 = 3 * i + 1;
        i2 = 3 * i + 2;
        for (i3 = 0; i3 < 6; i3++) {
          d = Kk[3 * i3];
          loop_ub = 3 * i3 + 1;
          Y_size_idx_0 = 3 * i3 + 2;
          b_k = i + 6 * i3;
          b_cxx[b_k] = (c_cxx[b_k] + ((b_Vx_tmp[i] * d + coeffs * Kk[loop_ub]) +
            dV_idx_0 * Kk[Y_size_idx_0])) + ((Kk[3 * i] * Qux[3 * i3] + Kk[i1] *
            Qux[loop_ub]) + Kk[i2] * Qux[Y_size_idx_0]);
          b_fx[b_k] = (Qux[3 * i] * d + Qux[i1] * Kk[loop_ub]) + Qux[i2] *
            Kk[Y_size_idx_0];
        }
      }

      for (i = 0; i < 36; i++) {
        Vxx[i] = b_cxx[i] + b_fx[i];
      }

      for (i = 0; i < 6; i++) {
        for (i1 = 0; i1 < 6; i1++) {
          b_k = i1 + 6 * i;
          b_fx[b_k] = 0.5 * (Vxx[b_k] + Vxx[i + 6 * i1]);
        }
      }

      memcpy(&Vxx[0], &b_fx[0], 36U * sizeof(double));

      /*  Make sure Hessian is symmetric */
      /*  Update Control Vectors */
      l[u_lims_tmp] = -lk[0];
      l[b_u_lims_tmp] = -lk[1];
      l[c_u_lims_tmp] = -lk[2];
      for (i = 0; i < 6; i++) {
        b_k = 3 * i + 18 * (498 - k);
        K[b_k] = -Kk[3 * i];
        K[b_k + 1] = -Kk[3 * i + 1];
        K[b_k + 2] = -Kk[3 * i + 2];
      }

      k++;
    }
  }
}

/*
 * Arguments    : const double Quu[9]
 *                const double Qu[3]
 *                const double lower[3]
 *                const double upper[3]
 *                const double u0[3]
 *                double u[3]
 *                double *result
 *                double Luu_data[]
 *                int Luu_size[2]
 *                boolean_T b_free[3]
 * Return Type  : void
 */
static void boxQPsolve(const double Quu[9], const double Qu[3], const double
  lower[3], const double upper[3], const double u0[3], double u[3], double
  *result, double Luu_data[], int Luu_size[2], boolean_T b_free[3])
{
  boolean_T clamped[3];
  double oldvalue;
  double scale;
  double z1[3];
  double absxk;
  double t;
  int i;
  double value;
  int iter;
  int b_iter;
  boolean_T exitg1;
  int k;
  boolean_T y;
  double grad[3];
  boolean_T prev_clamped[3];
  boolean_T exitg2;
  boolean_T factorize;
  boolean_T guard1 = false;
  int trueCount;
  signed char tmp_data[3];
  signed char b_tmp_data[3];
  int Quu_size[2];
  double Quu_data[9];
  double indef;
  int i1;
  double gnorm;
  double grad_clamped[3];
  double deltaX[3];
  signed char c_tmp_data[3];
  int n;
  int coeffs_size_idx_0;
  double Y_data[3];
  double coeffs_data[3];
  int X_size_idx_0;
  int LT_size_idx_0;
  int b_i;
  int i2;
  int i3;
  double LT_data[9];
  double sdotg;
  double step;
  double uc_idx_0;
  double uc_idx_1;
  double uc_idx_2;
  double vc;

  /*  Finds the optimal control with limits to minimize a quadratic cost */
  /*  Minimize 0.5*u'*Quu*u + u'*Qu  s.t. lower <= u <= upper */
  /*  */
  /*   inputs: */
  /*      Quu       - positive definite matrix              (m * m) */
  /*      Qu        - bias vector                           (m) */
  /*      lower     - lower bounds                          (m) */
  /*      upper     - upper bounds                          (m) */
  /*      u0        - initial control input for warm-start  (m) */
  /*  */
  /*   outputs: */
  /*      u         - solution                   (m) */
  /*      result    - result type (roughly, higher is better, see below) */
  /*      Luu       - cholesky factor            (m * m) */
  /*      free      - set of free dimensions     (m) */
  /*  Initialize arrays */
  clamped[0] = false;
  b_free[0] = true;
  clamped[1] = false;
  b_free[1] = true;
  clamped[2] = false;
  b_free[2] = true;
  Luu_size[0] = 3;
  Luu_size[1] = 3;
  memset(&Luu_data[0], 0, 9U * sizeof(double));

  /*  Placeholder to return if Luu not assigned */
  /*  Initialize scalars */
  oldvalue = 0.0;
  *result = 0.0;

  /*  options */
  /*  maximum number of iterations */
  /*  minimum norm of non-fixed gradient */
  /*  minimum relative improvement */
  /*  factor for decreasing stepsize */
  /*  minimal stepsize for linesearch */
  /*  Armijo parameter (fraction of linear improvement required) */
  /*  initial state */
  /*  Returns array x with all values clamped between lower and upper */
  /*  initial objective value */
  scale = fmax(lower[0], fmin(upper[0], u0[0]));
  z1[0] = scale;
  u[0] = scale;
  absxk = scale * Qu[0];
  scale = fmax(lower[1], fmin(upper[1], u0[1]));
  z1[1] = scale;
  u[1] = scale;
  absxk += scale * Qu[1];
  scale = fmax(lower[2], fmin(upper[2], u0[2]));
  z1[2] = scale;
  u[2] = scale;
  absxk += scale * Qu[2];
  t = 0.0;
  for (i = 0; i < 3; i++) {
    t += ((0.5 * z1[0] * Quu[3 * i] + 0.5 * z1[1] * Quu[3 * i + 1]) + 0.5 *
          scale * Quu[3 * i + 2]) * z1[i];
  }

  value = absxk + t;

  /*  main loop */
  iter = 1;
  b_iter = 1;
  exitg1 = false;
  while ((!exitg1) && (b_iter - 1 < 100)) {
    iter = b_iter;
    if (*result != 0.0) {
      exitg1 = true;
    } else {
      /*  check relative improvement */
      if ((b_iter > 1) && (oldvalue - value < 1.0E-8 * fabs(oldvalue))) {
        *result = 3.0;
        exitg1 = true;
      } else {
        oldvalue = value;

        /*  get gradient */
        /*  find clamped dimensions */
        for (k = 0; k < 3; k++) {
          scale = Qu[k] + ((Quu[k] * u[0] + Quu[k + 3] * u[1]) + Quu[k + 6] * u
                           [2]);
          grad[k] = scale;
          prev_clamped[k] = clamped[k];
          y = false;
          clamped[k] = false;
          if ((u[k] == lower[k]) && (scale > 0.0)) {
            y = true;
            clamped[k] = true;
          }

          if ((u[k] == upper[k]) && (scale < 0.0)) {
            y = true;
            clamped[k] = true;
          }

          b_free[k] = !y;
        }

        /*  check for all clamped */
        y = true;
        k = 0;
        exitg2 = false;
        while ((!exitg2) && (k < 3)) {
          if (!clamped[k]) {
            y = false;
            exitg2 = true;
          } else {
            k++;
          }
        }

        if (y) {
          *result = 5.0;
          exitg1 = true;
        } else {
          /*  Cholesky factorize if clamped controls have changed */
          if (b_iter == 1) {
            factorize = true;
          } else {
            factorize = false;
            k = 0;
            exitg2 = false;
            while ((!exitg2) && (k < 3)) {
              if (prev_clamped[k] != clamped[k]) {
                factorize = true;
                exitg2 = true;
              } else {
                k++;
              }
            }
          }

          /*  Cholesky (check for non PD) */
          guard1 = false;
          if (factorize) {
            trueCount = 0;
            if (b_free[0]) {
              trueCount = 1;
            }

            if (b_free[1]) {
              trueCount++;
            }

            if (b_free[2]) {
              trueCount++;
            }

            k = 0;
            if (b_free[0]) {
              tmp_data[0] = 1;
              k = 1;
            }

            if (b_free[1]) {
              tmp_data[k] = 2;
              k++;
            }

            if (b_free[2]) {
              tmp_data[k] = 3;
            }

            Quu_size[0] = trueCount;
            Quu_size[1] = trueCount;
            for (i = 0; i < trueCount; i++) {
              for (i1 = 0; i1 < trueCount; i1++) {
                Quu_data[i1 + trueCount * i] = Quu[(tmp_data[i1] + 3 *
                  (tmp_data[i] - 1)) - 1];
              }
            }

            chol_free(Quu_data, Quu_size, Luu_data, Luu_size, &indef);
            if (indef != 0.0) {
              *result = -1.0;
              exitg1 = true;
            } else {
              guard1 = true;
            }
          } else {
            guard1 = true;
          }

          if (guard1) {
            /*  check gradient norm */
            trueCount = 0;
            if (b_free[0]) {
              trueCount = 1;
            }

            if (b_free[1]) {
              trueCount++;
            }

            if (b_free[2]) {
              trueCount++;
            }

            k = 0;
            if (b_free[0]) {
              b_tmp_data[0] = 1;
              k = 1;
            }

            if (b_free[1]) {
              b_tmp_data[k] = 2;
              k++;
            }

            if (b_free[2]) {
              b_tmp_data[k] = 3;
            }

            if (trueCount == 0) {
              gnorm = 0.0;
            } else {
              gnorm = 0.0;
              if (trueCount == 1) {
                gnorm = fabs(grad[b_tmp_data[0] - 1]);
              } else {
                scale = 3.3121686421112381E-170;
                for (k = 0; k < trueCount; k++) {
                  absxk = fabs(grad[b_tmp_data[k] - 1]);
                  if (absxk > scale) {
                    t = scale / absxk;
                    gnorm = gnorm * t * t + 1.0;
                    scale = absxk;
                  } else {
                    t = absxk / scale;
                    gnorm += t * t;
                  }
                }

                gnorm = scale * sqrt(gnorm);
              }
            }

            if (gnorm < 1.0E-8) {
              *result = 4.0;
              exitg1 = true;
            } else {
              /*  get search direction */
              scale = u[0] * (double)clamped[0];
              absxk = u[1] * (double)clamped[1];
              t = u[2] * (double)clamped[2];
              for (k = 0; k < 3; k++) {
                grad_clamped[k] = Qu[k] + ((Quu[k] * scale + Quu[k + 3] * absxk)
                  + Quu[k + 6] * t);
                deltaX[k] = 0.0;
              }

              k = 0;
              if (b_free[0]) {
                c_tmp_data[0] = 1;
                k = 1;
              }

              if (b_free[1]) {
                c_tmp_data[k] = 2;
                k++;
              }

              if (b_free[2]) {
                c_tmp_data[k] = 3;
              }

              /*  Solves the linear system AX = B for the unknown X using the Cholesky */
              /*  decomposition L of the matrix A. Where LL' = A */
              /*  X can be a vector or a matrix of size n x m */
              /*  Solution is found in O(nm) time using back-substitution */
              /*  This implementation only works for lower triangular factorisations */
              n = Luu_size[0];

              /*  Solve LY = B for Y */
              /* ======================= */
              trueCount = Luu_size[0];
              coeffs_size_idx_0 = Luu_size[0];
              if (0 <= trueCount - 1) {
                memset(&Y_data[0], 0, trueCount * sizeof(double));
                memset(&coeffs_data[0], 0, trueCount * sizeof(double));
              }

              i = Luu_size[0];

              /*  Solve L'X = Y for X */
              /* ======================= */
              X_size_idx_0 = Luu_size[0];
              LT_size_idx_0 = Luu_size[1];
              for (k = 0; k < i; k++) {
                if (k + 1 != 1) {
                  for (i1 = 0; i1 < k; i1++) {
                    coeffs_data[i1] = Luu_data[k + Luu_size[0] * i1];
                  }
                }

                if (((signed char)coeffs_size_idx_0 == 1) || (Luu_size[0] == 1))
                {
                  scale = 0.0;
                  for (i1 = 0; i1 < coeffs_size_idx_0; i1++) {
                    scale += coeffs_data[i1] * Y_data[i1];
                  }
                } else {
                  scale = 0.0;
                  for (i1 = 0; i1 < coeffs_size_idx_0; i1++) {
                    scale += coeffs_data[i1] * Y_data[i1];
                  }
                }

                Y_data[k] = (grad_clamped[c_tmp_data[k] - 1] - scale) /
                  Luu_data[k + Luu_size[0] * k];
                z1[k] = 0.0;
                trueCount = Luu_size[1];
                for (i1 = 0; i1 < trueCount; i1++) {
                  LT_data[i1 + LT_size_idx_0 * k] = Luu_data[k + Luu_size[0] *
                    i1];
                }
              }

              coeffs_size_idx_0 = Luu_size[0];
              trueCount = Luu_size[0];
              if (0 <= trueCount - 1) {
                memset(&coeffs_data[0], 0, trueCount * sizeof(double));
              }

              i = (int)(((-1.0 - (double)Luu_size[0]) + 1.0) / -1.0);
              for (k = 0; k < i; k++) {
                b_i = (n - k) - 1;
                if (b_i + 1 != n) {
                  if (b_i + 2 > n) {
                    i1 = 0;
                    i2 = 0;
                    i3 = 0;
                  } else {
                    i1 = b_i + 1;
                    i2 = n;
                    i3 = b_i + 1;
                  }

                  trueCount = i2 - i1;
                  for (i2 = 0; i2 < trueCount; i2++) {
                    coeffs_data[i3 + i2] = LT_data[b_i + LT_size_idx_0 * (i1 +
                      i2)];
                  }
                }

                if (((signed char)coeffs_size_idx_0 == 1) || (X_size_idx_0 == 1))
                {
                  scale = 0.0;
                  for (i1 = 0; i1 < coeffs_size_idx_0; i1++) {
                    scale += coeffs_data[i1] * z1[i1];
                  }
                } else {
                  scale = 0.0;
                  for (i1 = 0; i1 < coeffs_size_idx_0; i1++) {
                    scale += coeffs_data[i1] * z1[i1];
                  }
                }

                z1[b_i] = (Y_data[b_i] - scale) / LT_data[b_i + LT_size_idx_0 *
                  b_i];
              }

              for (i = 0; i < X_size_idx_0; i++) {
                z1[i] = -z1[i];
              }

              k = 0;

              /*  cholesky solver */
              /*  check for descent direction */
              if (b_free[0]) {
                deltaX[0] = z1[0] - u[0];
                k = 1;
              }

              if (b_free[1]) {
                deltaX[1] = z1[k] - u[1];
                k++;
              }

              if (b_free[2]) {
                deltaX[2] = z1[k] - u[2];
              }

              sdotg = (deltaX[0] * grad[0] + deltaX[1] * grad[1]) + deltaX[2] *
                grad[2];
              if (sdotg >= 0.0) {
                /*  (should not happen) */
                exitg1 = true;
              } else {
                /*  Armijo linesearch */
                step = 1.0;

                /*  Returns array x with all values clamped between lower and upper */
                scale = fmax(lower[0], fmin(upper[0], u[0] + deltaX[0]));
                z1[0] = scale;
                uc_idx_0 = scale;
                absxk = scale * Qu[0];
                scale = fmax(lower[1], fmin(upper[1], u[1] + deltaX[1]));
                z1[1] = scale;
                uc_idx_1 = scale;
                absxk += scale * Qu[1];
                scale = fmax(lower[2], fmin(upper[2], u[2] + deltaX[2]));
                z1[2] = scale;
                uc_idx_2 = scale;
                absxk += scale * Qu[2];
                t = 0.0;
                for (i = 0; i < 3; i++) {
                  t += ((0.5 * z1[0] * Quu[3 * i] + 0.5 * z1[1] * Quu[3 * i + 1])
                        + 0.5 * scale * Quu[3 * i + 2]) * z1[i];
                }

                vc = absxk + t;
                exitg2 = false;
                while ((!exitg2) && ((vc - value) / (step * sdotg) < 0.1)) {
                  step *= 0.6;

                  /*  Returns array x with all values clamped between lower and upper */
                  scale = fmax(lower[0], fmin(upper[0], u[0] + step * deltaX[0]));
                  z1[0] = scale;
                  uc_idx_0 = scale;
                  absxk = scale * Qu[0];
                  scale = fmax(lower[1], fmin(upper[1], u[1] + step * deltaX[1]));
                  z1[1] = scale;
                  uc_idx_1 = scale;
                  absxk += scale * Qu[1];
                  scale = fmax(lower[2], fmin(upper[2], u[2] + step * deltaX[2]));
                  z1[2] = scale;
                  uc_idx_2 = scale;
                  absxk += scale * Qu[2];
                  t = 0.0;
                  for (i = 0; i < 3; i++) {
                    t += ((0.5 * z1[0] * Quu[3 * i] + 0.5 * z1[1] * Quu[3 * i +
                           1]) + 0.5 * scale * Quu[3 * i + 2]) * z1[i];
                  }

                  vc = absxk + t;
                  if (step < 1.0E-20) {
                    *result = 2.0;
                    exitg2 = true;
                  }
                }

                /*  accept candidate */
                u[0] = uc_idx_0;
                u[1] = uc_idx_1;
                u[2] = uc_idx_2;
                value = vc;
                b_iter++;
              }
            }
          }
        }
      }
    }
  }

  if (iter >= 100) {
    *result = 1.0;
  }

  /*  Results */
  /*  =========================== */
  /*  -1: Hessian is not positive definite */
  /*   0: No descent direction found          (SHOULD NOT OCCUR) */
  /*   1: Maximum main iterations exceeded        */
  /*   2: Maximum line-search iterations exceeded   */
  /*   3: Improvement smaller than tolerance      */
  /*   4: Gradient norm smaller than tolerance     */
  /*   5: All dimensions are clamped  */
}

/*
 * Arguments    : const double A_data[]
 *                const int A_size[2]
 *                double L_data[]
 *                int L_size[2]
 *                double *fail
 * Return Type  : void
 */
static void chol_free(const double A_data[], const int A_size[2], double L_data[],
                      int L_size[2], double *fail)
{
  int A_size_idx_0;
  int jmax;
  int idxAjj;
  double b_A_data[9];
  int n;
  int info;
  int b_info;
  int j;
  boolean_T exitg1;
  int idxA1j;
  double ssq;
  int ix;
  int iy;
  int i;
  int nmj;
  int i1;
  int ia0;
  double c_A_data[9];
  double c;

  /*  Wrapper for MATLAB chol for use with auto coder */
  /*  Inputs: */
  /* =========== */
  /*  A - positive semi-definite matrix */
  A_size_idx_0 = A_size[0];
  jmax = A_size[1];
  idxAjj = A_size[0] * A_size[1];
  if (0 <= idxAjj - 1) {
    memcpy(&b_A_data[0], &A_data[0], idxAjj * sizeof(double));
  }

  n = A_size[1];
  info = 0;
  if (A_size[1] != 0) {
    b_info = -1;
    j = 0;
    exitg1 = false;
    while ((!exitg1) && (j <= n - 1)) {
      idxA1j = j * n;
      idxAjj = idxA1j + j;
      ssq = 0.0;
      if (j >= 1) {
        ix = idxA1j;
        iy = idxA1j;
        for (jmax = 0; jmax < j; jmax++) {
          ssq += b_A_data[ix] * b_A_data[iy];
          ix++;
          iy++;
        }
      }

      ssq = b_A_data[idxAjj] - ssq;
      if (ssq > 0.0) {
        ssq = sqrt(ssq);
        b_A_data[idxAjj] = ssq;
        if (j + 1 < n) {
          nmj = (n - j) - 2;
          ia0 = (idxA1j + n) + 1;
          idxAjj += n;
          if ((j != 0) && (nmj + 1 != 0)) {
            iy = idxAjj;
            i = ia0 + n * nmj;
            for (jmax = ia0; n < 0 ? jmax >= i : jmax <= i; jmax += n) {
              ix = idxA1j;
              c = 0.0;
              i1 = (jmax + j) - 1;
              for (info = jmax; info <= i1; info++) {
                c += b_A_data[info - 1] * b_A_data[ix];
                ix++;
              }

              b_A_data[iy] += -c;
              iy += n;
            }
          }

          ssq = 1.0 / ssq;
          i = (idxAjj + n * nmj) + 1;
          for (jmax = idxAjj + 1; n < 0 ? jmax >= i : jmax <= i; jmax += n) {
            b_A_data[jmax - 1] *= ssq;
          }
        }

        j++;
      } else {
        b_A_data[idxAjj] = ssq;
        b_info = j;
        exitg1 = true;
      }
    }

    info = b_info + 1;
    if (b_info + 1 == 0) {
      jmax = A_size[1];
    } else {
      jmax = b_info;
    }

    for (j = 0; j < jmax; j++) {
      i = j + 2;
      for (idxAjj = i; idxAjj <= jmax; idxAjj++) {
        b_A_data[(idxAjj + A_size_idx_0 * j) - 1] = 0.0;
      }
    }

    if (1 > jmax) {
      idxAjj = 0;
      jmax = 0;
    } else {
      idxAjj = jmax;
    }

    for (i = 0; i < jmax; i++) {
      for (i1 = 0; i1 < idxAjj; i1++) {
        c_A_data[i1 + idxAjj * i] = b_A_data[i1 + A_size_idx_0 * i];
      }
    }

    A_size_idx_0 = idxAjj;
    idxAjj *= jmax;
    if (0 <= idxAjj - 1) {
      memcpy(&b_A_data[0], &c_A_data[0], idxAjj * sizeof(double));
    }
  }

  *fail = info;
  L_size[0] = A_size_idx_0;
  L_size[1] = jmax;
  idxAjj = A_size_idx_0 * jmax;
  if (0 <= idxAjj - 1) {
    memcpy(&L_data[0], &b_A_data[0], idxAjj * sizeof(double));
  }
}

/*
 * Arguments    : const double x[3500]
 *                const double xg[7]
 *                const double u[1497]
 *                const double l[1497]
 *                const double K[8982]
 *                double alpha
 *                const double u_lims[6]
 *                double xnew[3500]
 *                double unew[1497]
 *                double fx[17964]
 *                double fu[8982]
 *                double cx[3000]
 *                double cu[1497]
 *                double cxx[18000]
 *                double cuu[4491]
 *                double *cost
 * Return Type  : void
 */
static void forwardRollout(const double x[3500], const double xg[7], const
  double u[1497], const double l[1497], const double K[8982], double alpha,
  const double u_lims[6], double xnew[3500], double unew[1497], double fx[17964],
  double fu[8982], double cx[3000], double cu[1497], double cxx[18000], double
  cuu[4491], double *cost)
{
  int i;
  int k;
  double z1[3];
  int dx_tmp;
  double dx[6];
  int b_dx_tmp;
  double c;
  double dq[3];
  int c_dx_tmp;
  int q_tmp;
  int b_q_tmp;
  int c_q_tmp;
  double q[16];
  double q_error;
  double xnew_tmp;
  double d;
  double b_q_error[4];
  double d1;
  double d2;
  double d3;
  int b_k;
  double qdot_tmp[9];
  double wdot_tmp[9];
  double b_wdot_tmp[9];
  double c_wdot_tmp[9];
  int i1;
  int wdot_tmp_tmp;
  double dxdot1[70];
  double a[9];
  double b_q[4];
  double b_dv[12];
  short i2;
  static const short b_a[9] = { -200, 0, 0, 0, -200, 0, 0, 0, -200 };

  short i3;
  int dxdot1_tmp;
  double b_x[7];
  double c_a[3];
  static const signed char d_a[9] = { 100, 0, 0, 0, 100, 0, 0, 0, 100 };

  static const signed char b_iv[21] = { 0, 0, 0, 0, 100, 0, 0, 0, 0, 0, 0, 0,
    100, 0, 0, 0, 0, 0, 0, 0, 100 };

  double dxdot2[70];
  double x1[7];
  double E1[42];
  static const signed char iv1[18] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1,
    0, 0, 0, 1 };

  signed char b_I[49];
  double fx_tmp[49];
  double d4;
  double c_I[49];
  double b_xnew[42];
  double b_E1[42];
  double b_dv1[21];

  /*  Uses an rk method to roll out a trajectory */
  /*  Returns the new trajectory, cost and the derivatives along the trajectory */
  /*  If feed-forward and feed-back controls l, K are non-zero */
  /*  then the unew returned is the new control sequeunce. Otherwise unew = u */
  /*  Sizes */
  /*  Initialize outputs (Don't overwrite x or u) */
  memset(&xnew[0], 0, 3500U * sizeof(double));
  memset(&unew[0], 0, 1497U * sizeof(double));
  memset(&cx[0], 0, 3000U * sizeof(double));
  memset(&cxx[0], 0, 18000U * sizeof(double));
  *cost = 0.0;
  for (i = 0; i < 7; i++) {
    xnew[i] = x[i];
  }

  for (k = 0; k < 499; k++) {
    /*  Find the state error vector dx */
    dx_tmp = 7 * k + 4;
    dx[3] = xnew[dx_tmp] - x[dx_tmp];
    b_dx_tmp = 7 * k + 5;
    dx[4] = xnew[b_dx_tmp] - x[b_dx_tmp];
    c_dx_tmp = 7 * k + 6;
    dx[5] = xnew[c_dx_tmp] - x[c_dx_tmp];

    /*  Calculate error between qnew and q_nom */
    /*  Defined as conj(q_nom)*qnew */
    q_tmp = 7 * k + 1;
    b_q_tmp = 7 * k + 2;
    c_q_tmp = 7 * k + 3;

    /*  conjugate */
    /*  Forms the "left matrix" for quaternion multiplication */
    /*  where q*q1 = L(q)q1 */
    /*  Input */
    /*  q - the quaternion to build the left matrix from */
    /*  Output */
    /*  L - the left matrix */
    /*  L = [s, -v'; */
    /*       v, sI+skew(v)] */
    /* --------------------------------------------------- */
    q[0] = x[7 * k];
    q[4] = x[q_tmp];
    q[8] = x[b_q_tmp];
    q[12] = x[c_q_tmp];
    q[1] = -x[q_tmp];
    q[5] = x[7 * k];
    q[9] = x[c_q_tmp];
    q[13] = -x[b_q_tmp];
    q[2] = -x[b_q_tmp];
    q[6] = -x[c_q_tmp];
    q[10] = x[7 * k];
    q[14] = x[q_tmp];
    q[3] = -x[c_q_tmp];
    q[7] = x[b_q_tmp];
    q[11] = -x[q_tmp];
    q[15] = x[7 * k];
    q_error = 0.0;
    xnew_tmp = xnew[7 * k];
    for (i = 0; i < 4; i++) {
      d = ((q[i] * xnew_tmp + q[i + 4] * xnew[q_tmp]) + q[i + 8] * xnew[b_q_tmp])
        + q[i + 12] * xnew[c_q_tmp];
      b_q_error[i] = d;
      q_error += d * d;
    }

    q_error = sqrt(q_error);
    d = b_q_error[0] / q_error;
    d1 = b_q_error[1] / q_error;
    d2 = b_q_error[2] / q_error;
    d3 = b_q_error[3] / q_error;

    /*  re-normalize */
    /*  inverse Cayley Map */
    dx[0] = d1 / d;
    dx[1] = d2 / d;
    dx[2] = d3 / d;

    /*  Find the new control and ensure it is within the limits */
    /*  Step the dynamics forward */
    /*  Steps the dynamics forward using a 2nd order rk-method */
    /*  Returns: new state, discrete time Jacobians */
    /*  Step dynamics */
    /* Explicit midpoint step from x_k to x_{k+1} */
    /*  Calculates the continuous time state derivative and Jacobians */
    /*  kgm^2 */
    /*  Angular velocity */
    /*  Quaternion components */
    /*  Non-linear dynamics */
    for (b_k = 0; b_k < 3; b_k++) {
      d = 0.0;
      for (i = 0; i < 6; i++) {
        d += K[(b_k + 3 * i) + 18 * k] * dx[i];
      }

      i = b_k + 3 * k;
      d = fmin(u_lims[b_k + 3], fmax(u_lims[b_k], (u[i] - alpha * l[i]) - d));
      unew[i] = d;
      dq[b_k] = -xnew[(b_k + 7 * k) + 1];
    }

    /*  Returns skew symmetric - cross product matrix of a 3x1 vector */
    qdot_tmp[0] = 0.0;
    qdot_tmp[3] = -xnew[c_q_tmp];
    qdot_tmp[6] = xnew[b_q_tmp];
    qdot_tmp[1] = xnew[c_q_tmp];
    qdot_tmp[4] = 0.0;
    qdot_tmp[7] = -xnew[q_tmp];
    qdot_tmp[2] = -xnew[b_q_tmp];
    qdot_tmp[5] = xnew[q_tmp];
    qdot_tmp[8] = 0.0;
    for (i = 0; i < 9; i++) {
      qdot_tmp[i] += xnew_tmp * (double)iv[i];
    }

    /*  Returns skew symmetric - cross product matrix of a 3x1 vector */
    wdot_tmp[0] = 0.0;
    wdot_tmp[3] = -xnew[c_dx_tmp];
    wdot_tmp[6] = xnew[b_dx_tmp];
    wdot_tmp[1] = xnew[c_dx_tmp];
    wdot_tmp[4] = 0.0;
    wdot_tmp[7] = -xnew[dx_tmp];
    wdot_tmp[2] = -xnew[b_dx_tmp];
    wdot_tmp[5] = xnew[dx_tmp];
    wdot_tmp[8] = 0.0;

    /*  Jacobians  */
    for (i = 0; i < 3; i++) {
      d = wdot_tmp[i + 3];
      d1 = wdot_tmp[i + 6];
      d2 = 0.0;
      for (i1 = 0; i1 < 3; i1++) {
        wdot_tmp_tmp = i + 3 * i1;
        c_wdot_tmp[wdot_tmp_tmp] = (wdot_tmp[i] * dv[3 * i1] + d * dv[3 * i1 + 1])
          + d1 * dv[3 * i1 + 2];
        d2 += dv[wdot_tmp_tmp] * xnew[(i1 + 7 * k) + 4];
      }

      z1[i] = d2;
    }

    /*  Returns skew symmetric - cross product matrix of a 3x1 vector */
    /* !!!!! This changes slightly with magnetorquers !!!! */
    b_wdot_tmp[0] = c_wdot_tmp[0];
    b_wdot_tmp[3] = c_wdot_tmp[3] - (-z1[2]);
    b_wdot_tmp[6] = c_wdot_tmp[6] - z1[1];
    b_wdot_tmp[1] = c_wdot_tmp[1] - z1[2];
    b_wdot_tmp[4] = c_wdot_tmp[4];
    b_wdot_tmp[7] = c_wdot_tmp[7] - (-z1[0]);
    b_wdot_tmp[2] = c_wdot_tmp[2] - (-z1[1]);
    b_wdot_tmp[5] = c_wdot_tmp[5] - z1[0];
    b_wdot_tmp[8] = c_wdot_tmp[8];
    dxdot1[0] = 0.0;
    for (i = 0; i < 3; i++) {
      d = xnew[(i + 7 * k) + 4];
      b_k = 7 * (i + 1);
      dxdot1[b_k] = 0.5 * -d;
      wdot_tmp_tmp = 7 * (i + 4);
      dxdot1[wdot_tmp_tmp] = 0.5 * dq[i];
      dxdot1[i + 1] = 0.5 * d;
      i2 = b_a[i + 3];
      i3 = b_a[i + 6];
      for (i1 = 0; i1 < 3; i1++) {
        a[i + 3 * i1] = ((double)b_a[i] * b_wdot_tmp[3 * i1] + (double)i2 *
                         b_wdot_tmp[3 * i1 + 1]) + (double)i3 * b_wdot_tmp[3 *
          i1 + 2];
        dxdot1_tmp = i1 + 3 * i;
        dxdot1[(i1 + b_k) + 1] = 0.5 * -wdot_tmp[dxdot1_tmp];
        dxdot1[(i1 + wdot_tmp_tmp) + 1] = 0.5 * qdot_tmp[dxdot1_tmp];
      }
    }

    for (i = 0; i < 4; i++) {
      dxdot1[7 * i + 4] = 0.0;
      dxdot1[7 * i + 5] = 0.0;
      dxdot1[7 * i + 6] = 0.0;
    }

    for (i = 0; i < 3; i++) {
      b_k = 7 * (i + 4);
      dxdot1[b_k + 4] = 0.5 * a[3 * i];
      wdot_tmp_tmp = 3 * i + 1;
      dxdot1[b_k + 5] = 0.5 * a[wdot_tmp_tmp];
      dxdot1_tmp = 3 * i + 2;
      dxdot1[b_k + 6] = 0.5 * a[dxdot1_tmp];
      for (i1 = 0; i1 < 7; i1++) {
        dxdot1[i1 + 7 * (i + 7)] = b_iv[i1 + 7 * i];
      }

      i1 = i << 2;
      b_dv[i1] = 0.5 * dq[i];
      b_dv[i1 + 1] = 0.5 * qdot_tmp[3 * i];
      b_dv[i1 + 2] = 0.5 * qdot_tmp[wdot_tmp_tmp];
      b_dv[i1 + 3] = 0.5 * qdot_tmp[dxdot1_tmp];
    }

    for (i = 0; i < 4; i++) {
      b_q[i] = (b_dv[i] * xnew[dx_tmp] + b_dv[i + 4] * xnew[b_dx_tmp]) + b_dv[i
        + 8] * xnew[c_dx_tmp];
    }

    for (i = 0; i < 3; i++) {
      z1[i] = unew[i + 3 * k] - ((c_wdot_tmp[i] * xnew[dx_tmp] + c_wdot_tmp[i +
        3] * xnew[b_dx_tmp]) + c_wdot_tmp[i + 6] * xnew[c_dx_tmp]);
    }

    for (i = 0; i < 3; i++) {
      c_a[i] = ((double)d_a[i] * z1[0] + (double)d_a[i + 3] * z1[1]) + (double)
        d_a[i + 6] * z1[2];
    }

    b_x[0] = xnew_tmp + 0.015 * b_q[0];
    d = xnew[q_tmp] + 0.015 * b_q[1];
    b_x[1] = d;
    d1 = xnew[b_q_tmp] + 0.015 * b_q[2];
    b_x[2] = d1;
    d2 = xnew[c_q_tmp] + 0.015 * b_q[3];
    b_x[3] = d2;

    /*  Calculates the continuous time state derivative and Jacobians */
    /*  kgm^2 */
    /*  Angular velocity */
    /*  Quaternion components */
    /*  Non-linear dynamics */
    b_x[4] = xnew[dx_tmp] + 0.015 * c_a[0];
    dq[0] = -d;
    b_x[5] = xnew[b_dx_tmp] + 0.015 * c_a[1];
    dq[1] = -d1;
    b_x[6] = xnew[c_dx_tmp] + 0.015 * c_a[2];
    dq[2] = -d2;

    /*  Returns skew symmetric - cross product matrix of a 3x1 vector */
    qdot_tmp[0] = 0.0;
    qdot_tmp[3] = -d2;
    qdot_tmp[6] = d1;
    qdot_tmp[1] = d2;
    qdot_tmp[4] = 0.0;
    qdot_tmp[7] = -d;
    qdot_tmp[2] = -d1;
    qdot_tmp[5] = d;
    qdot_tmp[8] = 0.0;
    for (i = 0; i < 9; i++) {
      qdot_tmp[i] += b_x[0] * (double)iv[i];
    }

    /*  Returns skew symmetric - cross product matrix of a 3x1 vector */
    wdot_tmp[0] = 0.0;
    wdot_tmp[3] = -b_x[6];
    wdot_tmp[6] = b_x[5];
    wdot_tmp[1] = b_x[6];
    wdot_tmp[4] = 0.0;
    wdot_tmp[7] = -b_x[4];
    wdot_tmp[2] = -b_x[5];
    wdot_tmp[5] = b_x[4];
    wdot_tmp[8] = 0.0;

    /*  Jacobians  */
    for (i = 0; i < 3; i++) {
      d = wdot_tmp[i + 3];
      d1 = wdot_tmp[i + 6];
      d2 = 0.0;
      for (i1 = 0; i1 < 3; i1++) {
        wdot_tmp_tmp = i + 3 * i1;
        c_wdot_tmp[wdot_tmp_tmp] = (wdot_tmp[i] * dv[3 * i1] + d * dv[3 * i1 + 1])
          + d1 * dv[3 * i1 + 2];
        d2 += dv[wdot_tmp_tmp] * b_x[i1 + 4];
      }

      z1[i] = d2;
    }

    /*  Returns skew symmetric - cross product matrix of a 3x1 vector */
    /* !!!!! This changes slightly with magnetorquers !!!! */
    b_wdot_tmp[0] = c_wdot_tmp[0];
    b_wdot_tmp[3] = c_wdot_tmp[3] - (-z1[2]);
    b_wdot_tmp[6] = c_wdot_tmp[6] - z1[1];
    b_wdot_tmp[1] = c_wdot_tmp[1] - z1[2];
    b_wdot_tmp[4] = c_wdot_tmp[4];
    b_wdot_tmp[7] = c_wdot_tmp[7] - (-z1[0]);
    b_wdot_tmp[2] = c_wdot_tmp[2] - (-z1[1]);
    b_wdot_tmp[5] = c_wdot_tmp[5] - z1[0];
    b_wdot_tmp[8] = c_wdot_tmp[8];
    dxdot2[0] = 0.0;
    for (i = 0; i < 3; i++) {
      q_error = b_x[i + 4];
      b_k = 7 * (i + 1);
      dxdot2[b_k] = 0.5 * -q_error;
      wdot_tmp_tmp = 7 * (i + 4);
      dxdot2[wdot_tmp_tmp] = 0.5 * dq[i];
      dxdot2[i + 1] = 0.5 * q_error;
      i2 = b_a[i + 3];
      i3 = b_a[i + 6];
      for (i1 = 0; i1 < 3; i1++) {
        a[i + 3 * i1] = ((double)b_a[i] * b_wdot_tmp[3 * i1] + (double)i2 *
                         b_wdot_tmp[3 * i1 + 1]) + (double)i3 * b_wdot_tmp[3 *
          i1 + 2];
        dxdot1_tmp = i1 + 3 * i;
        dxdot2[(i1 + b_k) + 1] = 0.5 * -wdot_tmp[dxdot1_tmp];
        dxdot2[(i1 + wdot_tmp_tmp) + 1] = 0.5 * qdot_tmp[dxdot1_tmp];
      }
    }

    for (i = 0; i < 4; i++) {
      dxdot2[7 * i + 4] = 0.0;
      dxdot2[7 * i + 5] = 0.0;
      dxdot2[7 * i + 6] = 0.0;
    }

    for (i = 0; i < 3; i++) {
      b_k = 7 * (i + 4);
      dxdot2[b_k + 4] = 0.5 * a[3 * i];
      wdot_tmp_tmp = 3 * i + 1;
      dxdot2[b_k + 5] = 0.5 * a[wdot_tmp_tmp];
      dxdot1_tmp = 3 * i + 2;
      dxdot2[b_k + 6] = 0.5 * a[dxdot1_tmp];
      for (i1 = 0; i1 < 7; i1++) {
        dxdot2[i1 + 7 * (i + 7)] = b_iv[i1 + 7 * i];
      }

      i1 = i << 2;
      b_dv[i1] = 0.5 * dq[i];
      b_dv[i1 + 1] = 0.5 * qdot_tmp[3 * i];
      b_dv[i1 + 2] = 0.5 * qdot_tmp[wdot_tmp_tmp];
      b_dv[i1 + 3] = 0.5 * qdot_tmp[dxdot1_tmp];
    }

    d = b_x[4];
    d1 = b_x[5];
    d2 = b_x[6];
    for (i = 0; i < 4; i++) {
      b_q[i] = (b_dv[i] * d + b_dv[i + 4] * d1) + b_dv[i + 8] * d2;
    }

    for (i = 0; i < 3; i++) {
      z1[i] = unew[i + 3 * k] - ((c_wdot_tmp[i] * d + c_wdot_tmp[i + 3] * d1) +
        c_wdot_tmp[i + 6] * d2);
    }

    for (i = 0; i < 3; i++) {
      c_a[i] = ((double)d_a[i] * z1[0] + (double)d_a[i + 3] * z1[1]) + (double)
        d_a[i + 6] * z1[2];
    }

    d = xnew_tmp + 0.03 * b_q[0];
    x1[0] = d;
    d1 = xnew[q_tmp] + 0.03 * b_q[1];
    x1[1] = d1;
    d2 = xnew[b_q_tmp] + 0.03 * b_q[2];
    x1[2] = d2;
    d3 = xnew[c_q_tmp] + 0.03 * b_q[3];
    x1[3] = d3;
    x1[4] = xnew[dx_tmp] + 0.03 * c_a[0];
    x1[5] = xnew[b_dx_tmp] + 0.03 * c_a[1];
    x1[6] = xnew[c_dx_tmp] + 0.03 * c_a[2];

    /*  Normalize the quaternion */
    q_error = sqrt(((d * d + d1 * d1) + d2 * d2) + d3 * d3);

    /*  Continuous time Jacobians */
    /*  Discrete time Jacobians */
    /*  First form the state error Jacobians E(x_k) & E(x_k+1) */
    /*  Returns skew symmetric - cross product matrix of a 3x1 vector */
    /*  Returns skew symmetric - cross product matrix of a 3x1 vector */
    qdot_tmp[0] = 0.0;
    qdot_tmp[3] = -d3;
    qdot_tmp[6] = d2;
    qdot_tmp[1] = d3;
    qdot_tmp[4] = 0.0;
    qdot_tmp[7] = -d1;
    qdot_tmp[2] = -d2;
    qdot_tmp[5] = d1;
    qdot_tmp[8] = 0.0;
    for (i = 0; i < 3; i++) {
      E1[7 * i] = -x1[i + 1];
      b_k = 7 * (i + 3);
      E1[b_k] = 0.0;
      E1[7 * i + 1] = d * (double)iv[3 * i] + qdot_tmp[3 * i];
      E1[b_k + 1] = 0.0;
      wdot_tmp_tmp = 3 * i + 1;
      E1[7 * i + 2] = d * (double)iv[wdot_tmp_tmp] + qdot_tmp[wdot_tmp_tmp];
      E1[b_k + 2] = 0.0;
      wdot_tmp_tmp = 3 * i + 2;
      E1[7 * i + 3] = d * (double)iv[wdot_tmp_tmp] + qdot_tmp[wdot_tmp_tmp];
      E1[b_k + 3] = 0.0;
    }

    for (i = 0; i < 6; i++) {
      E1[7 * i + 4] = iv1[3 * i];
      E1[7 * i + 5] = iv1[3 * i + 1];
      E1[7 * i + 6] = iv1[3 * i + 2];
    }

    for (i = 0; i < 7; i++) {
      for (i1 = 0; i1 < 7; i1++) {
        b_k = i1 + 7 * i;
        fx_tmp[b_k] = 0.00045 * dxdot2[b_k];
      }
    }

    for (i = 0; i < 49; i++) {
      b_I[i] = 0;
    }

    for (b_k = 0; b_k < 7; b_k++) {
      b_I[b_k + 7 * b_k] = 1;
    }

    for (i = 0; i < 7; i++) {
      for (i1 = 0; i1 < 7; i1++) {
        d4 = 0.0;
        for (b_k = 0; b_k < 7; b_k++) {
          d4 += fx_tmp[i + 7 * b_k] * dxdot1[b_k + 7 * i1];
        }

        b_k = i + 7 * i1;
        c_I[b_k] = ((double)b_I[b_k] + 0.03 * dxdot2[b_k]) + d4;
      }
    }

    qdot_tmp[0] = 0.0;
    qdot_tmp[3] = -xnew[c_q_tmp];
    qdot_tmp[6] = xnew[b_q_tmp];
    qdot_tmp[1] = xnew[c_q_tmp];
    qdot_tmp[4] = 0.0;
    qdot_tmp[7] = -xnew[q_tmp];
    qdot_tmp[2] = -xnew[b_q_tmp];
    qdot_tmp[5] = xnew[q_tmp];
    qdot_tmp[8] = 0.0;
    for (i = 0; i < 6; i++) {
      for (i1 = 0; i1 < 7; i1++) {
        d4 = 0.0;
        for (b_k = 0; b_k < 7; b_k++) {
          d4 += E1[b_k + 7 * i] * c_I[b_k + 7 * i1];
        }

        b_E1[i + 6 * i1] = d4;
      }
    }

    for (i = 0; i < 3; i++) {
      b_xnew[7 * i] = -xnew[(i + 7 * k) + 1];
      wdot_tmp_tmp = 7 * (i + 3);
      b_xnew[wdot_tmp_tmp] = 0.0;
      b_xnew[7 * i + 1] = xnew_tmp * (double)iv[3 * i] + qdot_tmp[3 * i];
      b_xnew[wdot_tmp_tmp + 1] = 0.0;
      b_k = 3 * i + 1;
      b_xnew[7 * i + 2] = xnew_tmp * (double)iv[b_k] + qdot_tmp[b_k];
      b_xnew[wdot_tmp_tmp + 2] = 0.0;
      b_k = 3 * i + 2;
      b_xnew[7 * i + 3] = xnew_tmp * (double)iv[b_k] + qdot_tmp[b_k];
      b_xnew[wdot_tmp_tmp + 3] = 0.0;
    }

    for (i = 0; i < 6; i++) {
      b_xnew[7 * i + 4] = iv1[3 * i];
      b_xnew[7 * i + 5] = iv1[3 * i + 1];
      b_xnew[7 * i + 6] = iv1[3 * i + 2];
    }

    for (i = 0; i < 6; i++) {
      for (i1 = 0; i1 < 6; i1++) {
        d4 = 0.0;
        for (b_k = 0; b_k < 7; b_k++) {
          d4 += b_E1[i + 6 * b_k] * b_xnew[b_k + 7 * i1];
        }

        fx[(i + 6 * i1) + 36 * k] = d4;
      }
    }

    for (i = 0; i < 7; i++) {
      for (i1 = 0; i1 < 3; i1++) {
        d4 = 0.0;
        for (b_k = 0; b_k < 7; b_k++) {
          d4 += fx_tmp[i + 7 * b_k] * dxdot1[b_k + 7 * (i1 + 7)];
        }

        b_dv1[i + 7 * i1] = 0.03 * dxdot2[i + 7 * (i1 + 7)] + d4;
      }
    }

    for (i = 0; i < 6; i++) {
      for (i1 = 0; i1 < 3; i1++) {
        d4 = 0.0;
        for (b_k = 0; b_k < 7; b_k++) {
          d4 += E1[b_k + 7 * i] * b_dv1[b_k + 7 * i1];
        }

        fu[(i + 6 * i1) + 18 * k] = d4;
      }
    }

    wdot_tmp_tmp = 7 * (k + 1);
    xnew[wdot_tmp_tmp] = d / q_error;
    xnew[wdot_tmp_tmp + 1] = d1 / q_error;
    xnew[wdot_tmp_tmp + 2] = d2 / q_error;
    xnew[wdot_tmp_tmp + 3] = d3 / q_error;

    /*  Calculate the cost */
    for (i = 0; i < 3; i++) {
      xnew[(i + wdot_tmp_tmp) + 4] = x1[i + 4];
      b_k = 3 * i + 9 * k;
      cuu[b_k] = dv1[3 * i];
      cuu[b_k + 1] = dv1[3 * i + 1];
      cuu[b_k + 2] = dv1[3 * i + 2];
    }

    satellite_cost(*(double (*)[7])&xnew[7 * k], xg, *(double (*)[3])&unew[3 * k],
                   0.0, &c, *(double (*)[6])&cx[6 * k], *(double (*)[3])&cu[3 *
                   k], *(double (*)[36])&cxx[36 * k]);
    *cost += c;
  }

  /*  Final cost */
  z1[0] = 0.0;
  z1[1] = 0.0;
  z1[2] = 0.0;
  satellite_cost(*(double (*)[7])&xnew[3493], xg, z1, 1.0, &c, *(double (*)[6])&
                 cx[2994], dq, *(double (*)[36])&cxx[17964]);
  *cost += c;
}

/*
 * Arguments    : const double x[7]
 *                const double xg[7]
 *                const double u[3]
 *                double terminal
 *                double *cost
 *                double cx[6]
 *                double cu[3]
 *                double cxx[36]
 * Return Type  : void
 */
static void satellite_cost(const double x[7], const double xg[7], const double
  u[3], double terminal, double *cost, double cx[6], double cu[3], double cxx[36])
{
  double Qw[9];
  int i;
  double xg_tmp;
  double quat_cost;
  int b_sign;
  double b[3];
  double d;
  double b_dv[3];
  int cxx_tmp;
  int sign_tmp;
  double b_dv1[9];
  double c_sign[12];
  int b_cxx_tmp;

  /*  Calculates the cost contribution of a given state and control  */
  /*  Also calculates the 2nd order expansion of the cost function */
  /*  Inputs */
  /* ===================================== */
  /*  x        - [quaternion; omega] (7x1) */
  /*  u        - (3x1) (passed in as zeros for final time-step/terminal cost) */
  /*  terminal - int 0 or 1 */
  /* --------------------------------------------------- */
  /*  control hessian */
  if (terminal != 0.0) {
    for (i = 0; i < 9; i++) {
      Qw[i] = iv[i];
    }

    /*  terminal angular velocity hessian  */
  } else {
    memcpy(&Qw[0], &dv[0], 9U * sizeof(double));

    /*  angular velocity hessian */
  }

  /*  Find the linear quaternion-error cost (qcost = min(1+qg'q, 1-qg'q)) */
  /*  Also record the sign for use in the backward pass */
  xg_tmp = ((xg[0] * x[0] + xg[1] * x[1]) + xg[2] * x[2]) + xg[3] * x[3];
  if (xg_tmp + 1.0 < 1.0 - xg_tmp) {
    quat_cost = xg_tmp + 1.0;
    b_sign = 1;
  } else {
    quat_cost = 1.0 - xg_tmp;
    b_sign = -1;
  }

  b[0] = x[4] - xg[4];
  b[1] = x[5] - xg[5];
  b[2] = x[6] - xg[6];
  d = 0.0;
  for (i = 0; i < 3; i++) {
    cxx_tmp = 3 * i + 1;
    sign_tmp = 3 * i + 2;
    d += ((0.5 * b[0] * Qw[3 * i] + 0.5 * b[1] * Qw[cxx_tmp]) + 0.5 * b[2] *
          Qw[sign_tmp]) * b[i];
    b_dv[i] = (0.5 * u[0] * dv1[3 * i] + 0.5 * u[1] * dv1[cxx_tmp]) + 0.5 * u[2]
      * dv1[sign_tmp];
  }

  *cost = (quat_cost + d) + ((b_dv[0] * u[0] + b_dv[1] * u[1]) + b_dv[2] * u[2]);

  /*  State cost Hessian */
  /*  State cost Jacobian */
  /*  Returns skew symmetric - cross product matrix of a 3x1 vector */
  b_dv1[0] = 0.0;
  b_dv1[3] = -x[3];
  b_dv1[6] = x[2];
  b_dv1[1] = x[3];
  b_dv1[4] = 0.0;
  b_dv1[7] = -x[1];
  b_dv1[2] = -x[2];
  b_dv1[5] = x[1];
  b_dv1[8] = 0.0;
  for (i = 0; i < 3; i++) {
    c_sign[i] = (double)b_sign * -x[i + 1];
    cxx[6 * i] = (double)(-b_sign * iv[3 * i]) * xg_tmp;
    cxx_tmp = 6 * (i + 3);
    cxx[cxx_tmp] = 0.0;
    cxx[6 * i + 3] = 0.0;
    cxx[cxx_tmp + 3] = Qw[3 * i];
    sign_tmp = 3 * (i + 1);
    c_sign[sign_tmp] = (double)b_sign * (x[0] * (double)iv[i] + b_dv1[i]);
    b_cxx_tmp = 3 * i + 1;
    cxx[6 * i + 1] = (double)(-b_sign * iv[b_cxx_tmp]) * xg_tmp;
    cxx[cxx_tmp + 1] = 0.0;
    cxx[6 * i + 4] = 0.0;
    cxx[cxx_tmp + 4] = Qw[b_cxx_tmp];
    c_sign[sign_tmp + 1] = (double)b_sign * (x[0] * (double)iv[i + 3] + b_dv1[i
      + 3]);
    b_cxx_tmp = 3 * i + 2;
    cxx[6 * i + 2] = (double)(-b_sign * iv[b_cxx_tmp]) * xg_tmp;
    cxx[cxx_tmp + 2] = 0.0;
    cxx[6 * i + 5] = 0.0;
    cxx[cxx_tmp + 5] = Qw[b_cxx_tmp];
    c_sign[sign_tmp + 2] = (double)b_sign * (x[0] * (double)iv[i + 6] + b_dv1[i
      + 6]);
  }

  /*  Control cost Hessian & Jacobian */
  for (i = 0; i < 3; i++) {
    cx[i] = ((c_sign[i] * xg[0] + c_sign[i + 3] * xg[1]) + c_sign[i + 6] * xg[2])
      + c_sign[i + 9] * xg[3];
    cx[i + 3] = (Qw[i] * b[0] + Qw[i + 3] * b[1]) + Qw[i + 6] * b[2];
    cu[i] = (dv1[i] * u[0] + dv1[i + 3] * u[1]) + dv1[i + 6] * u[2];
  }
}

/*
 * Arguments    : const double x0[3500]
 *                const double xg[7]
 *                const double u0[1497]
 *                const double u_lims[6]
 *                double x[3500]
 *                double u[1497]
 *                double K[8982]
 *                boolean_T *result
 * Return Type  : void
 */
void milqr(const double x0[3500], const double xg[7], const double u0[1497],
           const double u_lims[6], double x[3500], double u[1497], double K[8982],
           boolean_T *result)
{
  int k;
  double Alphas[11];
  double lambda;
  double dlambda;
  static double x_n[3500];
  static double u_n[1497];
  static double fx_n[17964];
  static double fu_n[8982];
  static double cx_n[3000];
  static double cu_n[1497];
  static double cxx_n[18000];
  static double cuu_n[4491];
  double cost_n;
  static double l[1497];
  double dV[2];
  static double y[1497];
  static double b_dv[8982];
  static double fx[17964];
  static double fu[8982];
  static double cx[3000];
  static double cu[1497];
  static double cxx[18000];
  static double cuu[4491];
  double cost;
  int iter;
  int b_iter;
  boolean_T exitg1;
  boolean_T backPassDone;
  int exitg2;
  boolean_T diverge;
  double b_y;
  double maxval[499];
  double d;
  boolean_T fwdPassDone;
  boolean_T exitg3;
  double expectedChange;
  double z;
  double dcost;

  /*  Solves finite horizon optimal control problem using a */
  /*  multiplicative iterative linear quadratic regulator */
  /*  Note that this solver will not work without control limits */
  /*  Inputs */
  /*  =========================================== */
  /*  x0 - The intial trajectory (n, N) */
  /*  */
  /*  xg - The goal state (n, 1) */
  /*  */
  /*  u0 - The initial control sequeunce (m, N-1) */
  /*  */
  /*  u_lims - The control limits (m, 2) (lower, upper) */
  /*  Outputs */
  /*  =========================================== */
  /*  x - Final nominal trajectory (n, N) */
  /*  */
  /*  u - Final open-loop controls (m, N-1) */
  /*  */
  /*  K - Feedback control gains (n-1, m, N-1)  */
  /*  */
  /*  result - Indicates convergence (boolean) */
  /*  Options (pass in as array) */
  /*  Timestep (Should be highest we can get away with) */
  /*  maximum iterations */
  /*  cost reduction exit tolerance */
  /*  gradient exit criterion */
  /*  lambda criterion for gradient exit */
  /*  minimum accepted cost reduction ratio */
  /*  maximum regularization parameter */
  /*  set lambda = 0 below this value */
  /*  amount to scale dlambda by */
  /*  CONSTANTS */
  for (k = 0; k < 11; k++) {
    Alphas[k] = pow(10.0, -0.3 * (double)k);
  }

  /*  line search param */
  lambda = 1.0;
  dlambda = 1.0;

  /*  error state size (3 param. error representation for attitude) */
  /*  Init matrices for update (otherwise MATLAB coder throws an error) */
  memset(&x_n[0], 0, 3500U * sizeof(double));
  memset(&u_n[0], 0, 1497U * sizeof(double));
  memset(&fx_n[0], 0, 17964U * sizeof(double));
  memset(&fu_n[0], 0, 8982U * sizeof(double));
  memset(&cx_n[0], 0, 3000U * sizeof(double));
  memset(&cu_n[0], 0, 1497U * sizeof(double));
  memset(&cxx_n[0], 0, 18000U * sizeof(double));
  memset(&cuu_n[0], 0, 4491U * sizeof(double));
  cost_n = 0.0;

  /*  Initial Forward rollout */
  memset(&l[0], 0, 1497U * sizeof(double));
  memset(&K[0], 0, 8982U * sizeof(double));
  dV[0] = 0.0;
  dV[1] = 0.0;
  memset(&y[0], 0, 1497U * sizeof(double));
  memset(&b_dv[0], 0, 8982U * sizeof(double));
  forwardRollout(x0, xg, u0, y, b_dv, 0.0, u_lims, x, u, fx, fu, cx, cu, cxx,
                 cuu, &cost);

  /*  Convergence check params */
  /*  Expected cost change */
  /*  Ratio of cost change to expected cost change */
  *result = false;
  iter = 1;
  b_iter = 0;
  exitg1 = false;
  while ((!exitg1) && (b_iter < 500)) {
    iter = b_iter + 1;

    /*  Backward Pass */
    /* ======================================= */
    backPassDone = false;
    do {
      exitg2 = 0;
      if (!backPassDone) {
        backwardPass(fx, fu, cx, cu, cxx, cuu, lambda, u_lims, u, l, K, dV,
                     &diverge);
        if (diverge) {
          /*  Increase regularization parameter (lambda) */
          dlambda = fmax(1.6 * dlambda, 1.6);
          lambda = fmax(lambda * dlambda, 1.0E-6);
          if (lambda > 1.0E+10) {
            exitg2 = 1;
          }
        } else {
          backPassDone = true;
        }
      } else {
        exitg2 = 1;
      }
    } while (exitg2 == 0);

    /*  Check gradient of control, defined as l/u */
    /*  Terminate if sufficiently small (success) */
    for (k = 0; k < 1497; k++) {
      y[k] = fabs(l[k]) / (fabs(u[k]) + 1.0);
    }

    for (k = 0; k < 499; k++) {
      b_y = y[3 * k];
      d = y[3 * k + 1];
      if (b_y < d) {
        b_y = d;
      }

      d = y[3 * k + 2];
      if (b_y < d) {
        b_y = d;
      }

      maxval[k] = b_y;
    }

    b_y = maxval[0];
    for (k = 0; k < 498; k++) {
      b_y += maxval[k + 1];
    }

    /*  Avg of max grad at each time step */
    if ((b_y / 499.0 < 0.0001) && (lambda < 1.0E-5)) {
      *result = true;
      exitg1 = true;
    } else {
      /*  Forward Line-Search */
      /* =========================================== */
      fwdPassDone = false;
      if (backPassDone) {
        k = 0;
        exitg3 = false;
        while ((!exitg3) && (k < 11)) {
          forwardRollout(x, xg, u, l, K, Alphas[k], u_lims, x_n, u_n, fx_n, fu_n,
                         cx_n, cu_n, cxx_n, cuu_n, &cost_n);
          expectedChange = -Alphas[k] * (dV[0] + Alphas[k] * dV[1]);
          if (expectedChange > 0.0) {
            z = (cost - cost_n) / expectedChange;
          } else {
            z = cost - cost_n;
            if (z < 0.0) {
              z = -1.0;
            } else {
              if (z > 0.0) {
                z = 1.0;
              }
            }
          }

          if (z > 0.0) {
            fwdPassDone = true;
            exitg3 = true;
          } else {
            k++;
          }
        }
      }

      /*  Parameter Updates */
      /* ============================================= */
      if (fwdPassDone) {
        /*  Decrease Lambda */
        dlambda = fmin(dlambda / 1.6, 0.625);
        lambda = lambda * dlambda * (double)(lambda > 1.0E-6);

        /*  set = 0 if lambda too small */
        dcost = cost - cost_n;

        /*  Update trajectory and controls */
        memcpy(&x[0], &x_n[0], 3500U * sizeof(double));
        memcpy(&u[0], &u_n[0], 1497U * sizeof(double));
        memcpy(&fx[0], &fx_n[0], 17964U * sizeof(double));
        memcpy(&fu[0], &fu_n[0], 8982U * sizeof(double));
        memcpy(&cx[0], &cx_n[0], 3000U * sizeof(double));
        memcpy(&cu[0], &cu_n[0], 1497U * sizeof(double));
        memcpy(&cxx[0], &cxx_n[0], 18000U * sizeof(double));
        memcpy(&cuu[0], &cuu_n[0], 4491U * sizeof(double));
        cost = cost_n;

        /*  Terminate if cost change sufficiently small */
        if (dcost < 1.0E-7) {
          *result = true;
          exitg1 = true;
        } else {
          b_iter++;
        }
      } else {
        /*  No cost reduction (based on z-value) */
        /*  Increase lambda */
        dlambda = fmax(1.6 * dlambda, 1.6);
        lambda = fmax(lambda * dlambda, 1.0E-6);
        if (lambda > 1.0E+10) {
          /*  Lambda too large - solver diverged */
          exitg1 = true;
        } else {
          b_iter++;
        }
      }
    }
  }

  if (iter == 500) {
    /*  Ddin't converge completely */
    *result = false;
  }
}

/*
 * File trailer for milqr.c
 *
 * [EOF]
 */
