/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: milqr_efficient.c
 *
 * MATLAB Coder version            : 4.1
 * C/C++ source code generated on  : 22-Apr-2020 16:18:57
 */

/* Include Files */
#include <math.h>
#include <string.h>
#include "milqr_efficient.h"
#include <stdio.h>

/* Variable Definitions */
static const signed char iv0[9] = { 1, 0, 0, 0, 1, 0, 0, 0, 1 };

static const double dv0[9] = { 0.001, 0.0, 0.0, 0.0, 0.001, 0.0, 0.0, 0.0, 0.001
};

static const signed char iv1[18] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0,
  0, 0, 1 };

/* Function Declarations */
static bool all(const bool x[3]);
static bool any(const bool x[3]);
static void b_abs(const double x[297], double y[297]);
static void b_sign(double *x);
static void backwardPass(double lambda, const double u_lims[6], const double u
  [297], const double x[700], const double xg[7], double dt, const double B_ECI
  [297], double l[297], double K[1782], double dV[2], bool *diverged);
static void boxQPsolve(const double Quu[9], const double Qu[3], const double
  lower_lim[3], const double upper_lim[3], const double u0[3], double u[3],
  double *result, double Luu_data[], int Luu_size[2], bool b_free[3]);
static void chol_free(const double A_data[], const int A_size[2], double L_data[],
                      int L_size[2], double *fail);
static void chol_solve(const double L_data[], const int L_size[2], const double
  B_data[], const int B_size[2], double X_data[], int X_size[2]);
static void clamp(const double x[3], const double lower[3], const double upper[3],
                  double clampedVals[3]);
static void eye(double b_I[49]);
static void forwardRollout(const double x[700], const double xg[7], const double
  u[297], const double l[297], const double K[1782], double alpha, const double
  u_lims[6], double dt, const double B_ECI[297], double xnew[700], double unew
  [297], double *cost);
static double mean(const double x[99]);
static void power(double y[11]);
static void qmult(const double q1[4], const double q2[4], double q_out[4]);
static void satellite_cost_derivatives(const double x[7], const double xg[7],
  const double u[3], double cx[6], double cu[3], double cxx[36]);
static double satellite_cost_efficient(const double x[7], const double xg[7],
  const double u[3], double terminal);
static void satellite_derivatives(const double x0[7], const double u0[3], double
  dt, const double B_ECI[3], double fx[36], double fu[18]);
static void satellite_dynamics(const double x[7], const double u[3], const
  double B_ECI[3], double xdot[7], double dxdot[70]);
static void updateLambda(double *lambda, double direction);

/* Function Definitions */

/*
 * Arguments    : const bool x[3]
 * Return Type  : bool
 */
static bool all(const bool x[3])
{
  bool y;
  int k;
  bool exitg1;
  y = true;
  k = 0;
  exitg1 = false;
  while ((!exitg1) && (k < 3)) {
    if (!x[k]) {
      y = false;
      exitg1 = true;
    } else {
      k++;
    }
  }

  return y;
}

/*
 * Arguments    : const bool x[3]
 * Return Type  : bool
 */
static bool any(const bool x[3])
{
  bool y;
  int k;
  bool exitg1;
  y = false;
  k = 0;
  exitg1 = false;
  while ((!exitg1) && (k < 3)) {
    if (x[k]) {
      y = true;
      exitg1 = true;
    } else {
      k++;
    }
  }

  return y;
}

/*
 * Arguments    : const double x[297]
 *                double y[297]
 * Return Type  : void
 */
static void b_abs(const double x[297], double y[297])
{
  int k;
  for (k = 0; k < 297; k++) {
    y[k] = fabs(x[k]);
  }
}

/*
 * Arguments    : double *x
 * Return Type  : void
 */
static void b_sign(double *x)
{
  if (*x < 0.0) {
    *x = -1.0;
  } else {
    if (*x > 0.0) {
      *x = 1.0;
    }
  }
}

/*
 * Arguments    : double lambda
 *                const double u_lims[6]
 *                const double u[297]
 *                const double x[700]
 *                const double xg[7]
 *                double dt
 *                const double B_ECI[297]
 *                double l[297]
 *                double K[1782]
 *                double dV[2]
 *                bool *diverged
 * Return Type  : void
 */
static void backwardPass(double lambda, const double u_lims[6], const double u
  [297], const double x[700], const double xg[7], double dt, const double B_ECI
  [297], double l[297], double K[1782], double dV[2], bool *diverged)
{
  double dv9[3];
  double Vx[6];
  double unusedU0[3];
  double Vxx[36];
  int k;
  bool exitg1;
  double fx[36];
  double fu[18];
  double cx[6];
  double cu[3];
  double cxx[36];
  int i5;
  int i6;
  double lk_idx_1;
  double Qu_tmp[18];
  double Qx_tmp[36];
  double Qu[3];
  int i7;
  int u_lims_tmp_tmp;
  double Quu[9];
  double b_Quu[9];
  double Quu_tmp[18];
  double b_u_lims[3];
  double d4;
  int b_Quu_tmp;
  int b_u_lims_tmp_tmp;
  double Qux[18];
  int c_u_lims_tmp_tmp;
  double lk[3];
  double result;
  double Luu_data[9];
  int Luu_size[2];
  bool b_free[3];
  double Kk[18];
  int trueCount;
  double Qux_data[18];
  signed char tmp_data[3];
  double b_cx[6];
  double c_Quu_tmp[6];
  int Qux_size[2];
  int tmp_size[2];
  double b_Qx_tmp[36];
  double b_cxx[36];
  signed char b_tmp_data[3];

  /*  Perfoms the LQR backward pass to find the optimal controls */
  /*  Solves a quadratic program (QP) at each timestep for the optimal */
  /*  controls given the control limits */
  /*  insert section where I evaluate the cost and dynamics derivatives */
  /*  function [cx,cu,cxx,cuu,cxu] = cost_derivatives(x,u) */
  /*  function [fx,fu] = state_derivatives(x,u) */
  /*  Initialize matrices (for C code, not needed in MATLAB) */
  memset(&l[0], 0, 297U * sizeof(double));
  memset(&K[0], 0, 1782U * sizeof(double));

  /*  Change in cost */
  dV[0] = 0.0;
  dV[1] = 0.0;

  /*  Set cost-to-go Jacobian and Hessian equal to final costs */
  dv9[0] = 0.0;
  dv9[1] = 0.0;
  dv9[2] = 0.0;
  satellite_cost_derivatives(*(double (*)[7])&x[693], xg, dv9, Vx, unusedU0, Vxx);
  *diverged = false;
  k = 0;
  exitg1 = false;
  while ((!exitg1) && (k < 99)) {
    /*  calculate dynamics derivatives */
    satellite_derivatives(*(double (*)[7])&x[7 * (98 - k)], *(double (*)[3])&u[3
                          * (98 - k)], dt, *(double (*)[3])&B_ECI[3 * (98 - k)],
                          fx, fu);

    /*  Define cost gradients and hessians */
    /*  convert cost and dynamics derivatives into functions */
    satellite_cost_derivatives(*(double (*)[7])&x[7 * (98 - k)], xg, *(double (*)
      [3])&u[3 * (98 - k)], cx, cu, cxx);
    for (i5 = 0; i5 < 6; i5++) {
      for (i6 = 0; i6 < 6; i6++) {
        Qx_tmp[i6 + 6 * i5] = fx[i5 + 6 * i6];
      }

      Qu_tmp[3 * i5] = fu[i5];
      Qu_tmp[1 + 3 * i5] = fu[i5 + 6];
      Qu_tmp[2 + 3 * i5] = fu[i5 + 12];
    }

    for (i5 = 0; i5 < 3; i5++) {
      lk_idx_1 = 0.0;
      for (i6 = 0; i6 < 6; i6++) {
        i7 = i5 + 3 * i6;
        lk_idx_1 += Qu_tmp[i7] * Vx[i6];
        Quu_tmp[i7] = 0.0;
        d4 = 0.0;
        for (b_Quu_tmp = 0; b_Quu_tmp < 6; b_Quu_tmp++) {
          d4 += Qu_tmp[i5 + 3 * b_Quu_tmp] * Vxx[b_Quu_tmp + 6 * i6];
        }

        Quu_tmp[i7] = d4;
      }

      Qu[i5] = cu[i5] + lk_idx_1;
      for (i6 = 0; i6 < 3; i6++) {
        lk_idx_1 = 0.0;
        for (i7 = 0; i7 < 6; i7++) {
          lk_idx_1 += Quu_tmp[i5 + 3 * i7] * fu[i7 + 6 * i6];
        }

        b_Quu_tmp = i5 + 3 * i6;
        b_Quu[b_Quu_tmp] = dv0[b_Quu_tmp] + lk_idx_1;
      }

      for (i6 = 0; i6 < 6; i6++) {
        b_Quu_tmp = i5 + 3 * i6;
        Qux[b_Quu_tmp] = 0.0;
        lk_idx_1 = 0.0;
        for (i7 = 0; i7 < 6; i7++) {
          lk_idx_1 += Quu_tmp[i5 + 3 * i7] * fx[i7 + 6 * i6];
        }

        Qux[b_Quu_tmp] = lk_idx_1;
      }
    }

    /*  Regularization (for Cholesky positive definiteness) */
    /*  Solve the Quadratic program with control limits */
    for (i5 = 0; i5 < 9; i5++) {
      Quu[i5] = b_Quu[i5] + (double)iv0[i5] * lambda;
    }

    u_lims_tmp_tmp = 3 * (98 - k);
    unusedU0[0] = u_lims[0] - u[u_lims_tmp_tmp];
    b_u_lims[0] = u_lims[3] - u[u_lims_tmp_tmp];
    i5 = 3 * ((int)fmin(99.0, (99.0 + -(double)k) + 1.0) - 1);
    dv9[0] = -l[i5];
    b_u_lims_tmp_tmp = 1 + u_lims_tmp_tmp;
    unusedU0[1] = u_lims[1] - u[b_u_lims_tmp_tmp];
    b_u_lims[1] = u_lims[4] - u[b_u_lims_tmp_tmp];
    dv9[1] = -l[1 + i5];
    c_u_lims_tmp_tmp = 2 + u_lims_tmp_tmp;
    unusedU0[2] = u_lims[2] - u[c_u_lims_tmp_tmp];
    b_u_lims[2] = u_lims[5] - u[c_u_lims_tmp_tmp];
    dv9[2] = -l[2 + i5];
    boxQPsolve(Quu, Qu, unusedU0, b_u_lims, dv9, lk, &result, Luu_data, Luu_size,
               b_free);
    if (result < 2.0) {
      *diverged = true;

      /*  fprintf('\nDiverged with lambda = %f\n',lambda); */
      exitg1 = true;
    } else {
      /*  Solve for feedback gains in non-clamped rows of u */
      /*  (using cholesky factor of Quu) */
      memset(&Kk[0], 0, 18U * sizeof(double));
      if (any(b_free)) {
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

        b_Quu_tmp = 0;
        if (b_free[0]) {
          tmp_data[0] = 1;
          b_Quu_tmp = 1;
        }

        if (b_free[1]) {
          tmp_data[b_Quu_tmp] = 2;
          b_Quu_tmp++;
        }

        if (b_free[2]) {
          tmp_data[b_Quu_tmp] = 3;
        }

        Qux_size[0] = trueCount;
        Qux_size[1] = 6;
        for (i5 = 0; i5 < 6; i5++) {
          for (i6 = 0; i6 < trueCount; i6++) {
            Qux_data[i6 + trueCount * i5] = Qux[(tmp_data[i6] + 3 * i5) - 1];
          }
        }

        chol_solve(Luu_data, Luu_size, Qux_data, Qux_size, Qu_tmp, tmp_size);
        b_Quu_tmp = tmp_size[0] * tmp_size[1] - 1;
        for (i5 = 0; i5 <= b_Quu_tmp; i5++) {
          Qu_tmp[i5] = -Qu_tmp[i5];
        }

        b_Quu_tmp = 0;
        if (b_free[0]) {
          b_tmp_data[0] = 1;
          b_Quu_tmp = 1;
        }

        if (b_free[1]) {
          b_tmp_data[b_Quu_tmp] = 2;
          b_Quu_tmp++;
        }

        if (b_free[2]) {
          b_tmp_data[b_Quu_tmp] = 3;
        }

        b_Quu_tmp = tmp_size[0];
        for (i5 = 0; i5 < 6; i5++) {
          for (i6 = 0; i6 < b_Quu_tmp; i6++) {
            Kk[(b_tmp_data[i6] + 3 * i5) - 1] = Qu_tmp[i6 + tmp_size[0] * i5];
          }
        }
      }

      /*  Update Cost to Go Jacobian and Hessian */
      for (i5 = 0; i5 < 3; i5++) {
        for (i6 = 0; i6 < 6; i6++) {
          Qux_data[i6 + 6 * i5] = Kk[i5 + 3 * i6];
        }
      }

      memcpy(&Qu_tmp[0], &Qux_data[0], 18U * sizeof(double));
      for (i5 = 0; i5 < 6; i5++) {
        for (i6 = 0; i6 < 3; i6++) {
          b_Quu_tmp = i5 + 6 * i6;
          Qux_data[b_Quu_tmp] = 0.0;
          Qux_data[b_Quu_tmp] = (Qu_tmp[i5] * b_Quu[3 * i6] + Qu_tmp[i5 + 6] *
            b_Quu[1 + 3 * i6]) + Qu_tmp[i5 + 12] * b_Quu[2 + 3 * i6];
        }
      }

      memcpy(&Quu_tmp[0], &Qux_data[0], 18U * sizeof(double));
      for (i5 = 0; i5 < 3; i5++) {
        for (i6 = 0; i6 < 6; i6++) {
          Qux_data[i6 + 6 * i5] = Qux[i5 + 3 * i6];
        }
      }

      for (i5 = 0; i5 < 6; i5++) {
        lk_idx_1 = 0.0;
        for (i6 = 0; i6 < 6; i6++) {
          lk_idx_1 += Qx_tmp[i5 + 6 * i6] * Vx[i6];
        }

        b_cx[i5] = cx[i5] + lk_idx_1;
        c_Quu_tmp[i5] = 0.0;
        c_Quu_tmp[i5] = (Quu_tmp[i5] * lk[0] + Quu_tmp[i5 + 6] * lk[1]) +
          Quu_tmp[i5 + 12] * lk[2];
      }

      for (i5 = 0; i5 < 6; i5++) {
        cx[i5] = (b_cx[i5] + c_Quu_tmp[i5]) + ((Qu_tmp[i5] * Qu[0] + Qu_tmp[i5 +
          6] * Qu[1]) + Qu_tmp[i5 + 12] * Qu[2]);
      }

      for (i5 = 0; i5 < 6; i5++) {
        Vx[i5] = cx[i5] + ((Qux_data[i5] * lk[0] + Qux_data[i5 + 6] * lk[1]) +
                           Qux_data[i5 + 12] * lk[2]);
        for (i6 = 0; i6 < 6; i6++) {
          b_Quu_tmp = i5 + 6 * i6;
          b_Qx_tmp[b_Quu_tmp] = 0.0;
          lk_idx_1 = 0.0;
          for (i7 = 0; i7 < 6; i7++) {
            lk_idx_1 += Qx_tmp[i5 + 6 * i7] * Vxx[i7 + 6 * i6];
          }

          b_Qx_tmp[b_Quu_tmp] = lk_idx_1;
        }

        for (i6 = 0; i6 < 6; i6++) {
          lk_idx_1 = 0.0;
          for (i7 = 0; i7 < 6; i7++) {
            lk_idx_1 += b_Qx_tmp[i5 + 6 * i7] * fx[i7 + 6 * i6];
          }

          trueCount = i5 + 6 * i6;
          b_cxx[trueCount] = cxx[trueCount] + lk_idx_1;
        }
      }

      for (i5 = 0; i5 < 6; i5++) {
        for (i6 = 0; i6 < 6; i6++) {
          i7 = 1 + 3 * i6;
          b_Quu_tmp = 2 + 3 * i6;
          trueCount = i5 + 6 * i6;
          cxx[trueCount] = (b_cxx[trueCount] + ((Quu_tmp[i5] * Kk[3 * i6] +
            Quu_tmp[i5 + 6] * Kk[i7]) + Quu_tmp[i5 + 12] * Kk[b_Quu_tmp])) +
            ((Qu_tmp[i5] * Qux[3 * i6] + Qu_tmp[i5 + 6] * Qux[i7]) + Qu_tmp[i5 +
             12] * Qux[b_Quu_tmp]);
        }
      }

      for (i5 = 0; i5 < 6; i5++) {
        for (i6 = 0; i6 < 6; i6++) {
          b_Quu_tmp = i5 + 6 * i6;
          Qx_tmp[b_Quu_tmp] = 0.0;
          Qx_tmp[b_Quu_tmp] = (Qux_data[i5] * Kk[3 * i6] + Qux_data[i5 + 6] *
                               Kk[1 + 3 * i6]) + Qux_data[i5 + 12] * Kk[2 + 3 *
            i6];
        }
      }

      for (i5 = 0; i5 < 36; i5++) {
        Vxx[i5] = cxx[i5] + Qx_tmp[i5];
      }

      for (i5 = 0; i5 < 6; i5++) {
        for (i6 = 0; i6 < 6; i6++) {
          b_Quu_tmp = i6 + 6 * i5;
          Qx_tmp[b_Quu_tmp] = 0.5 * (Vxx[b_Quu_tmp] + Vxx[i5 + 6 * i6]);
        }
      }

      memcpy(&Vxx[0], &Qx_tmp[0], 36U * sizeof(double));

      /*  Ensure Hessian is symmetric */
      /*  Record control cost change to check convergence */
      lk_idx_1 = 0.0;
      for (i5 = 0; i5 < 3; i5++) {
        lk_idx_1 += ((0.5 * lk[0] * b_Quu[3 * i5] + 0.5 * lk[1] * b_Quu[1 + 3 *
                      i5]) + 0.5 * lk[2] * b_Quu[2 + 3 * i5]) * lk[i5];
      }

      dV[0] += (lk[0] * Qu[0] + lk[1] * Qu[1]) + lk[2] * Qu[2];
      dV[1] += lk_idx_1;

      /*  Update Control Vectors */
      l[u_lims_tmp_tmp] = -lk[0];
      l[b_u_lims_tmp_tmp] = -lk[1];
      l[c_u_lims_tmp_tmp] = -lk[2];
      for (i5 = 0; i5 < 6; i5++) {
        b_Quu_tmp = 3 * i5 + 18 * (98 - k);
        K[b_Quu_tmp] = -Kk[3 * i5];
        K[b_Quu_tmp + 1] = -Kk[1 + 3 * i5];
        K[b_Quu_tmp + 2] = -Kk[2 + 3 * i5];
      }

      k++;
    }
  }
}

/*
 * Arguments    : const double Quu[9]
 *                const double Qu[3]
 *                const double lower_lim[3]
 *                const double upper_lim[3]
 *                const double u0[3]
 *                double u[3]
 *                double *result
 *                double Luu_data[]
 *                int Luu_size[2]
 *                bool b_free[3]
 * Return Type  : void
 */
static void boxQPsolve(const double Quu[9], const double Qu[3], const double
  lower_lim[3], const double upper_lim[3], const double u0[3], double u[3],
  double *result, double Luu_data[], int Luu_size[2], bool b_free[3])
{
  bool clamped[3];
  double old_cost;
  double scale;
  int i12;
  double cost;
  int iter;
  int b_iter;
  bool exitg1;
  int i;
  double grad[3];
  bool prev_clamped[3];
  bool b0;
  bool factorize;
  bool b_clamped;
  int k;
  bool exitg2;
  bool guard1 = false;
  int trueCount;
  signed char tmp_data[3];
  signed char b_tmp_data[3];
  int Quu_size[2];
  double Quu_data[9];
  double indef;
  int i13;
  double grad_norm;
  double absxk;
  double t;
  double grad_clamped[3];
  signed char c_tmp_data[3];
  int n;
  double Y_data[3];
  double coeffs_data[3];
  int a_size_idx_1;
  int loop_ub;
  int b_loop_ub;
  double X_data[3];
  double a_data[3];
  int LT_size_idx_0;
  int c_loop_ub;
  double b_data[3];
  double LT_data[9];
  int d_loop_ub;
  int e_loop_ub;
  int b_i;
  double delta_u;
  double x[3];
  double b_delta_u[3];
  double expected_change;
  double step;
  double u_c_idx_0;
  double u_c_idx_1;
  double u_c_idx_2;
  double cost_c;

  /*  Copyright (c) 2020 Robotic Exploration Lab */
  /*   */
  /*  Permission is hereby granted, free of charge, to any person obtaining a copy */
  /*  of this software and associated documentation files (the "Software"), to deal */
  /*  in the Software without restriction, including without limitation the rights */
  /*  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell */
  /*  copies of the Software, and to permit persons to whom the Software is */
  /*  furnished to do so, subject to the following conditions: */
  /*   */
  /*  The above copyright notice and this permission notice shall be included in all */
  /*  copies or substantial portions of the Software. */
  /*   */
  /*  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR */
  /*  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, */
  /*  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE */
  /*  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER */
  /*  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, */
  /*  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE */
  /*  SOFTWARE. */
  /*  Finds the optimal control within limits to minimize a quadratic cost */
  /*  Minimizes 0.5*u'*Quu*u + u'*Qu  s.t. lower_lim <= u <= upper_lim */
  /*  */
  /*  Inputs: */
  /*  ========================================== */
  /*  Quu       - control cost Hessian (positive definite)   (m, m) */
  /*  Qu        - control cost Jacobian                      (m) */
  /*  lower     - control lower limit                        (m) */
  /*  upper     - control upper limit                        (m) */
  /*  u0        - initial control input for warm-start       (m) */
  /*  */
  /*  Outputs: */
  /*  ===================================== */
  /*  u         - optimal feed-forward control                (m) */
  /*  result    - gives exit criterion (explained below)          */
  /*  Luu       - cholesky factor                             (m, m) */
  /*  free      - set of free dimensions                      (m) */
  /*  Results */
  /*  =========================== */
  /*   0: No descent direction found */
  /*   1: Hessian is not positive definite */
  /*   2: Maximum main iterations exceeded        */
  /*   3: Maximum line-search iterations exceeded  */
  /*  */
  /*   4: Cost reduction smaller than tolerance      */
  /*   5: Gradient smaller than tolerance     */
  /*   6: All controls are clamped  */
  /*  Initialize arrays */
  /*  Indicies of clamped controls */
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
  old_cost = 0.0;
  *result = 0.0;

  /*  Solver Options */
  /*  max iterations */
  /*  min norm of non-clamped gradient */
  /*  min relative improvement */
  /*  factor for decreasing stepsize */
  /*  min stepsize for linesearch */
  /*  Armijo tolerance (fraction of linear improvement required) */
  /*  Initial controls */
  clamp(u0, lower_lim, upper_lim, u);

  /*  Initial cost value */
  scale = 0.0;
  for (i12 = 0; i12 < 3; i12++) {
    scale += ((0.5 * u[0] * Quu[3 * i12] + 0.5 * u[1] * Quu[1 + 3 * i12]) + 0.5 *
              u[2] * Quu[2 + 3 * i12]) * u[i12];
  }

  cost = ((u[0] * Qu[0] + u[1] * Qu[1]) + u[2] * Qu[2]) + scale;

  /*  Start optimisation */
  iter = 1;
  b_iter = 1;
  exitg1 = false;
  while ((!exitg1) && (b_iter - 1 < 100)) {
    iter = b_iter;
    if (*result != 0.0) {
      exitg1 = true;
    } else {
      /*  Check relative cost change for convergence */
      if ((b_iter > 1) && (old_cost - cost < 1.0E-8 * fabs(old_cost))) {
        *result = 4.0;
        exitg1 = true;
      } else {
        old_cost = cost;

        /*  Gradient of cost function */
        /*  Find clamped controls */
        for (i = 0; i < 3; i++) {
          scale = Qu[i] + ((Quu[i] * u[0] + Quu[i + 3] * u[1]) + Quu[i + 6] * u
                           [2]);
          grad[i] = scale;
          prev_clamped[i] = clamped[i];
          b0 = false;
          b_clamped = false;
          if ((u[i] == lower_lim[i]) && (scale > 0.0)) {
            b0 = true;
            b_clamped = true;
          }

          if ((u[i] == upper_lim[i]) && (scale < 0.0)) {
            b0 = true;
            b_clamped = true;
          }

          b_free[i] = !b0;
          clamped[i] = b_clamped;
        }

        /*  Check if all controls clamped */
        if (all(clamped)) {
          *result = 6.0;
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

          /*  Cholesky (check for non-PD) */
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
            for (i12 = 0; i12 < trueCount; i12++) {
              for (i13 = 0; i13 < trueCount; i13++) {
                Quu_data[i13 + trueCount * i12] = Quu[(tmp_data[i13] + 3 *
                  (tmp_data[i12] - 1)) - 1];
              }
            }

            chol_free(Quu_data, Quu_size, Luu_data, Luu_size, &indef);
            if (indef != 0.0) {
              *result = 1.0;
              exitg1 = true;
            } else {
              guard1 = true;
            }
          } else {
            guard1 = true;
          }

          if (guard1) {
            /*  check gradient-norm of free controls */
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
              grad_norm = 0.0;
            } else {
              grad_norm = 0.0;
              if (trueCount == 1) {
                grad_norm = fabs(grad[b_tmp_data[0] - 1]);
              } else {
                scale = 3.3121686421112381E-170;
                for (k = 0; k < trueCount; k++) {
                  absxk = fabs(grad[b_tmp_data[k] - 1]);
                  if (absxk > scale) {
                    t = scale / absxk;
                    grad_norm = 1.0 + grad_norm * t * t;
                    scale = absxk;
                  } else {
                    t = absxk / scale;
                    grad_norm += t * t;
                  }
                }

                grad_norm = scale * sqrt(grad_norm);
              }
            }

            if (grad_norm < 1.0E-8) {
              *result = 5.0;
              exitg1 = true;
            } else {
              /*  get search direction */
              scale = u[0] * (double)clamped[0];
              absxk = u[1] * (double)clamped[1];
              t = u[2] * (double)clamped[2];
              k = 0;
              for (i = 0; i < 3; i++) {
                grad_clamped[i] = Qu[i] + ((Quu[i] * scale + Quu[i + 3] * absxk)
                  + Quu[i + 6] * t);
                if (b_free[i]) {
                  c_tmp_data[k] = (signed char)(i + 1);
                  k++;
                }
              }

              /*  Copyright (c) 2020 Robotic Exploration Lab */
              /*   */
              /*  Permission is hereby granted, free of charge, to any person obtaining a copy */
              /*  of this software and associated documentation files (the "Software"), to deal */
              /*  in the Software without restriction, including without limitation the rights */
              /*  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell */
              /*  copies of the Software, and to permit persons to whom the Software is */
              /*  furnished to do so, subject to the following conditions: */
              /*   */
              /*  The above copyright notice and this permission notice shall be included in all */
              /*  copies or substantial portions of the Software. */
              /*   */
              /*  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR */
              /*  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, */
              /*  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE */
              /*  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER */
              /*  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, */
              /*  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE */
              /*  SOFTWARE. */
              /*  Solves the linear system AX = B for the unknown X using the Cholesky */
              /*  decomposition L of the matrix A. Where LL' = A */
              /*  X can be a vector or a matrix of size n x m */
              /*  Solution is found in O(nm) time using back-substitution */
              /*  This implementation only works for lower triangular factorisations */
              n = Luu_size[0];

              /*  Check sizes match and L is lower-triangular */
              /*  Solve LY = B for Y */
              /* ======================= */
              if (0 <= Luu_size[0] - 1) {
                memset(&Y_data[0], 0, (unsigned int)(Luu_size[0] * (int)sizeof
                        (double)));
              }

              if (0 <= Luu_size[0] - 1) {
                memset(&coeffs_data[0], 0, (unsigned int)(Luu_size[0] * (int)
                        sizeof(double)));
              }

              i12 = Luu_size[0];
              if (0 <= Luu_size[0] - 1) {
                a_size_idx_1 = Luu_size[0];
                loop_ub = Luu_size[0];
                b_loop_ub = Luu_size[0];
              }

              for (i = 0; i < i12; i++) {
                if (1 + i != 1) {
                  for (i13 = 0; i13 < i; i13++) {
                    coeffs_data[i13] = Luu_data[i + Luu_size[0] * i13];
                  }
                }

                if (0 <= loop_ub - 1) {
                  memcpy(&a_data[0], &coeffs_data[0], (unsigned int)(loop_ub *
                          (int)sizeof(double)));
                }

                if (0 <= b_loop_ub - 1) {
                  memcpy(&b_data[0], &Y_data[0], (unsigned int)(b_loop_ub * (int)
                          sizeof(double)));
                }

                if ((a_size_idx_1 == 1) || (Luu_size[0] == 1)) {
                  scale = 0.0;
                  for (i13 = 0; i13 < a_size_idx_1; i13++) {
                    scale += a_data[i13] * b_data[i13];
                  }
                } else {
                  scale = 0.0;
                  for (i13 = 0; i13 < a_size_idx_1; i13++) {
                    scale += a_data[i13] * b_data[i13];
                  }
                }

                Y_data[i] = (grad_clamped[c_tmp_data[i] - 1] - scale) /
                  Luu_data[i + Luu_size[0] * i];
              }

              /*  Solve L'X = Y for X */
              /* ======================= */
              if (0 <= Luu_size[0] - 1) {
                memset(&X_data[0], 0, (unsigned int)(Luu_size[0] * (int)sizeof
                        (double)));
              }

              LT_size_idx_0 = Luu_size[1];
              c_loop_ub = Luu_size[0];
              for (i12 = 0; i12 < c_loop_ub; i12++) {
                k = Luu_size[1];
                for (i13 = 0; i13 < k; i13++) {
                  LT_data[i13 + LT_size_idx_0 * i12] = Luu_data[i12 + Luu_size[0]
                    * i13];
                }
              }

              if (0 <= Luu_size[0] - 1) {
                memset(&coeffs_data[0], 0, (unsigned int)(Luu_size[0] * (int)
                        sizeof(double)));
              }

              i12 = (int)((1.0 + (-1.0 - (double)Luu_size[0])) / -1.0);
              if (0 <= i12 - 1) {
                a_size_idx_1 = Luu_size[0];
                d_loop_ub = Luu_size[0];
                e_loop_ub = Luu_size[0];
              }

              for (i = 0; i < i12; i++) {
                b_i = (n - i) - 1;
                if (b_i + 1 != n) {
                  if (b_i + 2 > n) {
                    i13 = -1;
                    k = 0;
                    trueCount = -1;
                  } else {
                    i13 = b_i;
                    k = n;
                    trueCount = b_i;
                  }

                  c_loop_ub = k - i13;
                  for (k = 0; k <= c_loop_ub - 2; k++) {
                    coeffs_data[(trueCount + k) + 1] = LT_data[b_i +
                      LT_size_idx_0 * ((i13 + k) + 1)];
                  }
                }

                if (0 <= d_loop_ub - 1) {
                  memcpy(&a_data[0], &coeffs_data[0], (unsigned int)(d_loop_ub *
                          (int)sizeof(double)));
                }

                if (0 <= e_loop_ub - 1) {
                  memcpy(&b_data[0], &X_data[0], (unsigned int)(e_loop_ub * (int)
                          sizeof(double)));
                }

                if ((a_size_idx_1 == 1) || (Luu_size[0] == 1)) {
                  scale = 0.0;
                  for (i13 = 0; i13 < a_size_idx_1; i13++) {
                    scale += a_data[i13] * b_data[i13];
                  }
                } else {
                  scale = 0.0;
                  for (i13 = 0; i13 < a_size_idx_1; i13++) {
                    scale += a_data[i13] * b_data[i13];
                  }
                }

                X_data[b_i] = (Y_data[b_i] - scale) / LT_data[b_i +
                  LT_size_idx_0 * b_i];
              }

              c_loop_ub = Luu_size[0];
              for (i12 = 0; i12 < c_loop_ub; i12++) {
                X_data[i12] = -X_data[i12];
              }

              k = 0;

              /*  cholesky solver */
              /*  check projected change in cost is a reduction */
              delta_u = 0.0;
              if (b_free[0]) {
                delta_u = X_data[0] - u[0];
                k = 1;
              }

              x[0] = delta_u * grad[0];
              b_delta_u[0] = delta_u;
              delta_u = 0.0;
              if (b_free[1]) {
                delta_u = X_data[k] - u[1];
                k++;
              }

              x[1] = delta_u * grad[1];
              b_delta_u[1] = delta_u;
              delta_u = 0.0;
              if (b_free[2]) {
                delta_u = X_data[k] - u[2];
              }

              expected_change = (x[0] + x[1]) + delta_u * grad[2];
              if (expected_change >= 0.0) {
                /*  (should not happen) */
                exitg1 = true;
              } else {
                /*  Armijo linesearch */
                step = 1.0;

                /*  Returns array x with all values clamped between lower and upper */
                scale = fmax(lower_lim[0], fmin(upper_lim[0], u[0] + b_delta_u[0]));
                x[0] = scale;
                u_c_idx_0 = scale;
                absxk = scale * Qu[0];
                scale = fmax(lower_lim[1], fmin(upper_lim[1], u[1] + b_delta_u[1]));
                x[1] = scale;
                u_c_idx_1 = scale;
                absxk += scale * Qu[1];
                scale = fmax(lower_lim[2], fmin(upper_lim[2], u[2] + delta_u));
                x[2] = scale;
                u_c_idx_2 = scale;
                absxk += scale * Qu[2];
                t = 0.0;
                for (i12 = 0; i12 < 3; i12++) {
                  t += ((0.5 * x[0] * Quu[3 * i12] + 0.5 * x[1] * Quu[1 + 3 *
                         i12]) + 0.5 * scale * Quu[2 + 3 * i12]) * x[i12];
                }

                cost_c = absxk + t;
                exitg2 = false;
                while ((!exitg2) && ((cost_c - cost) / (step * expected_change) <
                                     0.1)) {
                  step *= 0.6;

                  /*  Returns array x with all values clamped between lower and upper */
                  scale = fmax(lower_lim[0], fmin(upper_lim[0], u[0] + step *
                    b_delta_u[0]));
                  x[0] = scale;
                  u_c_idx_0 = scale;
                  absxk = scale * Qu[0];
                  scale = fmax(lower_lim[1], fmin(upper_lim[1], u[1] + step *
                    b_delta_u[1]));
                  x[1] = scale;
                  u_c_idx_1 = scale;
                  absxk += scale * Qu[1];
                  scale = fmax(lower_lim[2], fmin(upper_lim[2], u[2] + step *
                    delta_u));
                  x[2] = scale;
                  u_c_idx_2 = scale;
                  absxk += scale * Qu[2];
                  t = 0.0;
                  for (i12 = 0; i12 < 3; i12++) {
                    t += ((0.5 * x[0] * Quu[3 * i12] + 0.5 * x[1] * Quu[1 + 3 *
                           i12]) + 0.5 * scale * Quu[2 + 3 * i12]) * x[i12];
                  }

                  cost_c = absxk + t;
                  if (step < 1.0E-20) {
                    *result = 3.0;
                    exitg2 = true;
                  }
                }

                /*  Update the control */
                u[0] = u_c_idx_0;
                u[1] = u_c_idx_1;
                u[2] = u_c_idx_2;
                cost = cost_c;
                b_iter++;
              }
            }
          }
        }
      }
    }
  }

  if (iter >= 100) {
    *result = 2.0;
  }
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
  int nmj;
  int k;
  double b_A_data[9];
  int n;
  int info;
  int b_n;
  int j;
  bool exitg1;
  int jj;
  double ajj;
  int ix;
  int iy;
  int i14;
  int i15;
  double c_A_data[9];
  int jp1j;
  int iac;
  double c;
  int ia;

  /*  Wrapper for MATLAB chol for use with auto coder */
  /*  Inputs: */
  /* =========== */
  /*  A - positive semi-definite matrix */
  A_size_idx_0 = A_size[0];
  nmj = A_size[1];
  k = A_size[0] * A_size[1];
  if (0 <= k - 1) {
    memcpy(&b_A_data[0], &A_data[0], (unsigned int)(k * (int)sizeof(double)));
  }

  n = A_size[1];
  info = 0;
  if (A_size[1] != 0) {
    b_n = A_size[0];
    info = 0;
    if (A_size[0] != 0) {
      j = 0;
      exitg1 = false;
      while ((!exitg1) && (j <= b_n - 1)) {
        jj = j + j * n;
        ajj = 0.0;
        if (j >= 1) {
          ix = j;
          iy = j;
          for (k = 0; k < j; k++) {
            ajj += b_A_data[ix] * b_A_data[iy];
            ix += n;
            iy += n;
          }
        }

        ajj = b_A_data[jj] - ajj;
        if (ajj > 0.0) {
          ajj = sqrt(ajj);
          b_A_data[jj] = ajj;
          if (j + 1 < b_n) {
            nmj = (b_n - j) - 1;
            k = j + 2;
            jp1j = jj + 2;
            if ((nmj == 0) || (j == 0)) {
            } else {
              ix = j;
              i14 = (j + n * (j - 1)) + 2;
              for (iac = k; n < 0 ? iac >= i14 : iac <= i14; iac += n) {
                c = -b_A_data[ix];
                iy = jj + 1;
                i15 = (iac + nmj) - 1;
                for (ia = iac; ia <= i15; ia++) {
                  b_A_data[iy] += b_A_data[ia - 1] * c;
                  iy++;
                }

                ix += n;
              }
            }

            ajj = 1.0 / ajj;
            i14 = jj + nmj;
            for (k = jp1j; k <= i14 + 1; k++) {
              b_A_data[k - 1] *= ajj;
            }
          }

          j++;
        } else {
          b_A_data[jj] = ajj;
          info = j + 1;
          exitg1 = true;
        }
      }
    }

    if (info == 0) {
      nmj = A_size[1];
    } else {
      nmj = info - 1;
    }

    for (j = 2; j <= nmj; j++) {
      for (k = 0; k <= j - 2; k++) {
        b_A_data[k + A_size_idx_0 * (j - 1)] = 0.0;
      }
    }

    if (1 > nmj) {
      k = 0;
      nmj = 0;
    } else {
      k = nmj;
    }

    for (i14 = 0; i14 < nmj; i14++) {
      for (i15 = 0; i15 < k; i15++) {
        c_A_data[i15 + k * i14] = b_A_data[i15 + A_size_idx_0 * i14];
      }
    }

    A_size_idx_0 = k;
    k *= nmj;
    if (0 <= k - 1) {
      memcpy(&b_A_data[0], &c_A_data[0], (unsigned int)(k * (int)sizeof(double)));
    }
  }

  *fail = info;
  L_size[0] = A_size_idx_0;
  L_size[1] = nmj;
  k = A_size_idx_0 * nmj;
  if (0 <= k - 1) {
    memcpy(&L_data[0], &b_A_data[0], (unsigned int)(k * (int)sizeof(double)));
  }
}

/*
 * Arguments    : const double L_data[]
 *                const int L_size[2]
 *                const double B_data[]
 *                const int B_size[2]
 *                double X_data[]
 *                int X_size[2]
 * Return Type  : void
 */
static void chol_solve(const double L_data[], const int L_size[2], const double
  B_data[], const int B_size[2], double X_data[], int X_size[2])
{
  int n;
  int Y_size_idx_0;
  int loop_ub;
  double Y_data[18];
  double coeffs_data[3];
  int i16;
  int a_size_idx_1;
  int i;
  int b_loop_ub;
  int c_loop_ub;
  int i17;
  double a_data[3];
  int j;
  int LT_size_idx_0;
  double b_data[3];
  double a;
  double LT_data[9];
  int d_loop_ub;
  int e_loop_ub;
  int b_i;

  /*  Copyright (c) 2020 Robotic Exploration Lab */
  /*   */
  /*  Permission is hereby granted, free of charge, to any person obtaining a copy */
  /*  of this software and associated documentation files (the "Software"), to deal */
  /*  in the Software without restriction, including without limitation the rights */
  /*  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell */
  /*  copies of the Software, and to permit persons to whom the Software is */
  /*  furnished to do so, subject to the following conditions: */
  /*   */
  /*  The above copyright notice and this permission notice shall be included in all */
  /*  copies or substantial portions of the Software. */
  /*   */
  /*  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR */
  /*  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, */
  /*  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE */
  /*  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER */
  /*  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, */
  /*  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE */
  /*  SOFTWARE. */
  /*  Solves the linear system AX = B for the unknown X using the Cholesky */
  /*  decomposition L of the matrix A. Where LL' = A */
  /*  X can be a vector or a matrix of size n x m */
  /*  Solution is found in O(nm) time using back-substitution */
  /*  This implementation only works for lower triangular factorisations */
  n = L_size[0];

  /*  Check sizes match and L is lower-triangular */
  /*  Solve LY = B for Y */
  /* ======================= */
  Y_size_idx_0 = L_size[0];
  loop_ub = L_size[0] * 6;
  if (0 <= loop_ub - 1) {
    memset(&Y_data[0], 0, (unsigned int)(loop_ub * (int)sizeof(double)));
  }

  if (0 <= L_size[0] - 1) {
    memset(&coeffs_data[0], 0, (unsigned int)(L_size[0] * (int)sizeof(double)));
  }

  i16 = L_size[0];
  if (0 <= L_size[0] - 1) {
    a_size_idx_1 = L_size[0];
    b_loop_ub = L_size[0];
    c_loop_ub = L_size[0];
  }

  for (i = 0; i < i16; i++) {
    if (1 + i != 1) {
      for (i17 = 0; i17 < i; i17++) {
        coeffs_data[i17] = L_data[i + L_size[0] * i17];
      }
    }

    if (0 <= b_loop_ub - 1) {
      memcpy(&a_data[0], &coeffs_data[0], (unsigned int)(b_loop_ub * (int)sizeof
              (double)));
    }

    for (j = 0; j < 6; j++) {
      for (i17 = 0; i17 < c_loop_ub; i17++) {
        b_data[i17] = Y_data[i17 + Y_size_idx_0 * j];
      }

      if ((a_size_idx_1 == 1) || (Y_size_idx_0 == 1)) {
        a = 0.0;
        for (i17 = 0; i17 < a_size_idx_1; i17++) {
          a += a_data[i17] * b_data[i17];
        }
      } else {
        a = 0.0;
        for (i17 = 0; i17 < a_size_idx_1; i17++) {
          a += a_data[i17] * b_data[i17];
        }
      }

      Y_data[i + Y_size_idx_0 * j] = (B_data[i + B_size[0] * j] - a) / L_data[i
        + L_size[0] * i];
    }
  }

  /*  Solve L'X = Y for X */
  /* ======================= */
  X_size[0] = L_size[0];
  X_size[1] = 6;
  loop_ub = L_size[0] * 6;
  if (0 <= loop_ub - 1) {
    memset(&X_data[0], 0, (unsigned int)(loop_ub * (int)sizeof(double)));
  }

  LT_size_idx_0 = L_size[1];
  loop_ub = L_size[0];
  for (i16 = 0; i16 < loop_ub; i16++) {
    b_loop_ub = L_size[1];
    for (i17 = 0; i17 < b_loop_ub; i17++) {
      LT_data[i17 + LT_size_idx_0 * i16] = L_data[i16 + L_size[0] * i17];
    }
  }

  if (0 <= L_size[0] - 1) {
    memset(&coeffs_data[0], 0, (unsigned int)(L_size[0] * (int)sizeof(double)));
  }

  i16 = (int)((1.0 + (-1.0 - (double)L_size[0])) / -1.0);
  if (0 <= i16 - 1) {
    a_size_idx_1 = L_size[0];
    d_loop_ub = L_size[0];
    e_loop_ub = L_size[0];
  }

  for (i = 0; i < i16; i++) {
    b_i = (n - i) - 1;
    if (b_i + 1 != n) {
      if (b_i + 2 > n) {
        i17 = -1;
        c_loop_ub = 0;
        b_loop_ub = -1;
      } else {
        i17 = b_i;
        c_loop_ub = n;
        b_loop_ub = b_i;
      }

      loop_ub = c_loop_ub - i17;
      for (c_loop_ub = 0; c_loop_ub <= loop_ub - 2; c_loop_ub++) {
        coeffs_data[(b_loop_ub + c_loop_ub) + 1] = LT_data[b_i + LT_size_idx_0 *
          ((i17 + c_loop_ub) + 1)];
      }
    }

    if (0 <= d_loop_ub - 1) {
      memcpy(&a_data[0], &coeffs_data[0], (unsigned int)(d_loop_ub * (int)sizeof
              (double)));
    }

    for (j = 0; j < 6; j++) {
      for (i17 = 0; i17 < e_loop_ub; i17++) {
        b_data[i17] = X_data[i17 + X_size[0] * j];
      }

      if ((a_size_idx_1 == 1) || (X_size[0] == 1)) {
        a = 0.0;
        for (i17 = 0; i17 < a_size_idx_1; i17++) {
          a += a_data[i17] * b_data[i17];
        }
      } else {
        a = 0.0;
        for (i17 = 0; i17 < a_size_idx_1; i17++) {
          a += a_data[i17] * b_data[i17];
        }
      }

      X_data[b_i + X_size[0] * j] = (Y_data[b_i + Y_size_idx_0 * j] - a) /
        LT_data[b_i + LT_size_idx_0 * b_i];
    }
  }
}

/*
 * Arguments    : const double x[3]
 *                const double lower[3]
 *                const double upper[3]
 *                double clampedVals[3]
 * Return Type  : void
 */
static void clamp(const double x[3], const double lower[3], const double upper[3],
                  double clampedVals[3])
{
  /*  Returns array x with all values clamped between lower and upper */
  clampedVals[0] = fmax(lower[0], fmin(upper[0], x[0]));
  clampedVals[1] = fmax(lower[1], fmin(upper[1], x[1]));
  clampedVals[2] = fmax(lower[2], fmin(upper[2], x[2]));
}

/*
 * Arguments    : double b_I[49]
 * Return Type  : void
 */
static void eye(double b_I[49])
{
  int k;
  memset(&b_I[0], 0, 49U * sizeof(double));
  for (k = 0; k < 7; k++) {
    b_I[k + 7 * k] = 1.0;
  }
}

/*
 * Arguments    : const double x[700]
 *                const double xg[7]
 *                const double u[297]
 *                const double l[297]
 *                const double K[1782]
 *                double alpha
 *                const double u_lims[6]
 *                double dt
 *                const double B_ECI[297]
 *                double xnew[700]
 *                double unew[297]
 *                double *cost
 * Return Type  : void
 */
static void forwardRollout(const double x[700], const double xg[7], const double
  u[297], const double l[297], const double K[1782], double alpha, const double
  u_lims[6], double dt, const double B_ECI[297], double xnew[700], double unew
  [297], double *cost)
{
  int i1;
  double a;
  int k;
  double dv3[3];
  int dx_tmp;
  double dx[6];
  double q_idx_0;
  int b_k;
  int q_idx_3_tmp;
  double q[16];
  double q_error[4];
  double d1;
  double b[7];
  double dxdot1[70];
  double b_xnew[7];
  double dxdot2[70];
  double x1[7];

  /*  get rid of everything after unew except cost */
  /*  Uses an rk method to roll out a trajectory */
  /*  Returns the new trajectory, cost and the derivatives along the trajectory */
  /*  If feed-forward and feed-back controls l, K are non-zero */
  /*  then the unew returned is the new control sequeunce. Otherwise unew = u */
  /*  Sizes */
  /*  Initialize outputs (Don't overwrite x or u) */
  memset(&xnew[0], 0, 700U * sizeof(double));
  memset(&unew[0], 0, 297U * sizeof(double));
  *cost = 0.0;
  for (i1 = 0; i1 < 7; i1++) {
    xnew[i1] = x[i1];
  }

  a = 0.5 * dt;
  for (k = 0; k < 99; k++) {
    /*  Find the state error vector dx */
    dx_tmp = 7 * k + 4;
    dx[3] = xnew[dx_tmp] - x[dx_tmp];
    dx_tmp = 7 * k + 5;
    dx[4] = xnew[dx_tmp] - x[dx_tmp];
    dx_tmp = 7 * k + 6;
    dx[5] = xnew[dx_tmp] - x[dx_tmp];

    /*  Calculate error between qk and q_nom */
    /*  Defined as conj(q_nom)*qnew */
    /*  Returns error as Rodrigues parameters (3x1) */
    q_idx_0 = x[7 * k];
    dx_tmp = 1 + 7 * k;
    b_k = 2 + 7 * k;
    q_idx_3_tmp = 3 + 7 * k;

    /*  conjugate */
    /*  Copyright (c) 2020 Robotic Exploration Lab */
    /*   */
    /*  Permission is hereby granted, free of charge, to any person obtaining a copy */
    /*  of this software and associated documentation files (the "Software"), to deal */
    /*  in the Software without restriction, including without limitation the rights */
    /*  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell */
    /*  copies of the Software, and to permit persons to whom the Software is */
    /*  furnished to do so, subject to the following conditions: */
    /*   */
    /*  The above copyright notice and this permission notice shall be included in all */
    /*  copies or substantial portions of the Software. */
    /*   */
    /*  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR */
    /*  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, */
    /*  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE */
    /*  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER */
    /*  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, */
    /*  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE */
    /*  SOFTWARE. */
    /*  Forms the "left matrix" for quaternion multiplication */
    /*  where q*q1 = L(q)q1 */
    /*  Input */
    /*  q - the quaternion to build the left matrix from */
    /*  Output */
    /*  L - the left matrix */
    /*  L = [s, -v'; */
    /*       v, sI+skew(v)] */
    /* --------------------------------------------------- */
    q[0] = q_idx_0;
    q[4] = x[dx_tmp];
    q[8] = x[b_k];
    q[12] = x[q_idx_3_tmp];
    q[1] = -x[dx_tmp];
    q[5] = q_idx_0;
    q[9] = x[q_idx_3_tmp];
    q[13] = -x[b_k];
    q[2] = -x[b_k];
    q[6] = -x[q_idx_3_tmp];
    q[10] = q_idx_0;
    q[14] = x[dx_tmp];
    q[3] = -x[q_idx_3_tmp];
    q[7] = x[b_k];
    q[11] = -x[dx_tmp];
    q[15] = q_idx_0;
    q_idx_0 = 0.0;
    for (i1 = 0; i1 < 4; i1++) {
      q_error[i1] = 0.0;
      d1 = ((q[i1] * xnew[7 * k] + q[i1 + 4] * xnew[dx_tmp]) + q[i1 + 8] *
            xnew[b_k]) + q[i1 + 12] * xnew[q_idx_3_tmp];
      q_error[i1] = d1;
      q_idx_0 += d1 * d1;
    }

    q_idx_0 = sqrt(q_idx_0);
    q_error[0] /= q_idx_0;
    q_error[1] /= q_idx_0;
    q_error[2] /= q_idx_0;
    q_error[3] /= q_idx_0;

    /*  re-normalize */
    /*  inverse Cayley Map */
    dx[0] = q_error[1] / q_error[0];
    dx[1] = q_error[2] / q_error[0];
    dx[2] = q_error[3] / q_error[0];

    /*  Find the new control and ensure it is within the limits */
    for (b_k = 0; b_k < 3; b_k++) {
      d1 = 0.0;
      for (i1 = 0; i1 < 6; i1++) {
        d1 += K[(b_k + 3 * i1) + 18 * k] * dx[i1];
      }

      dx_tmp = b_k + 3 * k;
      unew[dx_tmp] = fmin(u_lims[3 + b_k], fmax(u_lims[b_k], (u[dx_tmp] - alpha *
        l[dx_tmp]) - d1));
    }

    /*  Step the dynamics forward */
    /*  recreate dynamics function to separate into dynamics step and  */
    /*  derivative function to get derivatives */
    /*  Copyright (c) 2020 Robotic Exploration Lab */
    /*   */
    /*  Permission is hereby granted, free of charge, to any person obtaining a copy */
    /*  of this software and associated documentation files (the "Software"), to deal */
    /*  in the Software without restriction, including without limitation the rights */
    /*  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell */
    /*  copies of the Software, and to permit persons to whom the Software is */
    /*  furnished to do so, subject to the following conditions: */
    /*   */
    /*  The above copyright notice and this permission notice shall be included in all */
    /*  copies or substantial portions of the Software. */
    /*   */
    /*  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR */
    /*  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, */
    /*  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE */
    /*  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER */
    /*  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, */
    /*  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE */
    /*  SOFTWARE. */
    /*  Steps the dynamics forward using a 2nd order rk-method */
    /*  Returns: new state, discrete time Jacobians */
    /*  Step dynamics */
    /*  Explicit midpoint step from x_k to x_{k+1} */
    satellite_dynamics(*(double (*)[7])&xnew[7 * k], *(double (*)[3])&unew[3 * k],
                       *(double (*)[3])&B_ECI[3 * k], b, dxdot1);
    for (i1 = 0; i1 < 7; i1++) {
      b_xnew[i1] = xnew[i1 + 7 * k] + a * b[i1];
    }

    satellite_dynamics(b_xnew, *(double (*)[3])&unew[3 * k], *(double (*)[3])&
                       B_ECI[3 * k], b, dxdot2);
    for (i1 = 0; i1 < 7; i1++) {
      x1[i1] = xnew[i1 + 7 * k] + dt * b[i1];
    }

    /*  Re-normalize the quaternion */
    q_idx_0 = sqrt(((x1[0] * x1[0] + x1[1] * x1[1]) + x1[2] * x1[2]) + x1[3] *
                   x1[3]);
    dx_tmp = 7 * (k + 1);
    xnew[dx_tmp] = x1[0] / q_idx_0;
    xnew[1 + dx_tmp] = x1[1] / q_idx_0;
    xnew[2 + dx_tmp] = x1[2] / q_idx_0;
    xnew[3 + dx_tmp] = x1[3] / q_idx_0;
    xnew[dx_tmp + 4] = x1[4];
    xnew[dx_tmp + 5] = x1[5];
    xnew[dx_tmp + 6] = x1[6];

    /*  Continuous time Jacobians */
    /*  Discrete time Jacobians */
    /*  First form the state error Jacobians E(x_k) & E(x_k+1) */
    /*  Calculate the cost */
    /*  make second function that's satellite cost derivatives */
    *cost += satellite_cost_efficient(*(double (*)[7])&xnew[7 * k], xg, *(double
      (*)[3])&unew[3 * k], 0.0);
  }

  /*  Final cost */
  dv3[0] = 0.0;
  dv3[1] = 0.0;
  dv3[2] = 0.0;
  *cost += satellite_cost_efficient(*(double (*)[7])&xnew[693], xg, dv3, 1.0);
}

/*
 * Arguments    : const double x[99]
 * Return Type  : double
 */
static double mean(const double x[99])
{
  double y;
  int k;
  y = x[0];
  for (k = 0; k < 98; k++) {
    y += x[k + 1];
  }

  y /= 99.0;
  return y;
}

/*
 * Arguments    : double y[11]
 * Return Type  : void
 */
static void power(double y[11])
{
  int k;
  for (k = 0; k < 11; k++) {
    y[k] = pow(10.0, -0.3 * (double)k);
  }
}

/*
 * Arguments    : const double q1[4]
 *                const double q2[4]
 *                double q_out[4]
 * Return Type  : void
 */
static void qmult(const double q1[4], const double q2[4], double q_out[4])
{
  /*  this function multiplies quaternions */
  q_out[0] = q1[0] * q2[0] - ((q1[1] * q2[1] + q1[2] * q2[2]) + q1[3] * q2[3]);
  q_out[1] = (q1[0] * q2[1] + q2[0] * q1[1]) + (q1[2] * q2[3] - q1[3] * q2[2]);
  q_out[2] = (q1[0] * q2[2] + q2[0] * q1[2]) + (q1[3] * q2[1] - q1[1] * q2[3]);
  q_out[3] = (q1[0] * q2[3] + q2[0] * q1[3]) + (q1[1] * q2[2] - q1[2] * q2[1]);
}

/*
 * Arguments    : const double x[7]
 *                const double xg[7]
 *                const double u[3]
 *                double cx[6]
 *                double cu[3]
 *                double cxx[36]
 * Return Type  : void
 */
static void satellite_cost_derivatives(const double x[7], const double xg[7],
  const double u[3], double cx[6], double cu[3], double cxx[36])
{
  double xg_tmp;
  double b_xg;
  int c_sign;
  int a;
  int i8;
  int cxx_tmp;
  double dv10[9];
  double b_a[12];
  double b_x[3];

  /*  Copyright (c) 2020 Robotic Exploration Lab */
  /*   */
  /*  Permission is hereby granted, free of charge, to any person obtaining a copy */
  /*  of this software and associated documentation files (the "Software"), to deal */
  /*  in the Software without restriction, including without limitation the rights */
  /*  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell */
  /*  copies of the Software, and to permit persons to whom the Software is */
  /*  furnished to do so, subject to the following conditions: */
  /*   */
  /*  The above copyright notice and this permission notice shall be included in all */
  /*  copies or substantial portions of the Software. */
  /*   */
  /*  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR */
  /*  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, */
  /*  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE */
  /*  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER */
  /*  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, */
  /*  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE */
  /*  SOFTWARE. */
  /*  Calculates the cost contribution of a given state and control */
  /*  Utilizes a standard quadratic cost function for the angular velocity */
  /*  and a geodesic cost for the attitude quaternion */
  /*  Also returns the 2nd order expansion of the cost function */
  /*  Inputs */
  /* ===================================== */
  /*  x        - [quaternion; omega] (7x1) */
  /*  u        - (3x1) (passed in as zeros for final time-step/terminal cost) */
  /*  terminal - integer 0 or 1 */
  /* --------------------------------------------------- */
  /*  control hessian */
  /*  terminal angular velocity cost hessian */
  /*  terminal geodesic cost weight */
  /*  Finds the geodesic quaternion-error cost */
  /*  quat_cost = min(1+qd'q, 1-qd'q)) where qd is the desired attitude quaternion  */
  /*  Also records the sign (+ or -) which minimizes quat_cost */
  /*  this is used when calculating the Jacobain and Hessian */
  xg_tmp = xg[0] * x[0];
  b_xg = ((xg_tmp + xg[1] * x[1]) + xg[2] * x[2]) + xg[3] * x[3];
  if (1.0 + b_xg < 1.0 - b_xg) {
    c_sign = 1;
  } else {
    c_sign = -1;
  }

  /*  State cost Hessian */
  a = -c_sign * 10;
  b_xg = ((xg_tmp + xg[1] * x[1]) + xg[2] * x[2]) + xg[3] * x[3];
  for (i8 = 0; i8 < 3; i8++) {
    cxx[6 * i8] = (double)(a * iv0[3 * i8]) * b_xg;
    cxx_tmp = 6 * (i8 + 3);
    cxx[cxx_tmp] = 0.0;
    cxx[1 + 6 * i8] = (double)(a * iv0[1 + 3 * i8]) * b_xg;
    cxx[1 + cxx_tmp] = 0.0;
    cxx[2 + 6 * i8] = (double)(a * iv0[2 + 3 * i8]) * b_xg;
    cxx[2 + cxx_tmp] = 0.0;
  }

  for (i8 = 0; i8 < 6; i8++) {
    cxx[6 * i8 + 3] = iv1[3 * i8];
    cxx[6 * i8 + 4] = iv1[1 + 3 * i8];
    cxx[6 * i8 + 5] = iv1[2 + 3 * i8];
  }

  /*  State cost Jacobian */
  /*  Copyright (c) 2020 Robotic Exploration Lab */
  /*   */
  /*  Permission is hereby granted, free of charge, to any person obtaining a copy */
  /*  of this software and associated documentation files (the "Software"), to deal */
  /*  in the Software without restriction, including without limitation the rights */
  /*  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell */
  /*  copies of the Software, and to permit persons to whom the Software is */
  /*  furnished to do so, subject to the following conditions: */
  /*   */
  /*  The above copyright notice and this permission notice shall be included in all */
  /*  copies or substantial portions of the Software. */
  /*   */
  /*  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR */
  /*  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, */
  /*  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE */
  /*  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER */
  /*  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, */
  /*  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE */
  /*  SOFTWARE. */
  /*  Returns skew symmetric - cross product matrix of a 3x1 vector */
  a = c_sign * 10;
  dv10[0] = 0.0;
  dv10[3] = -x[3];
  dv10[6] = x[2];
  dv10[1] = x[3];
  dv10[4] = 0.0;
  dv10[7] = -x[1];
  dv10[2] = -x[2];
  dv10[5] = x[1];
  dv10[8] = 0.0;
  for (i8 = 0; i8 < 3; i8++) {
    b_a[i8] = (double)a * -x[1 + i8];
    cxx_tmp = 3 * (i8 + 1);
    b_a[cxx_tmp] = (double)a * (x[0] * (double)iv0[i8] + dv10[i8]);
    b_a[1 + cxx_tmp] = (double)a * (x[0] * (double)iv0[i8 + 3] + dv10[i8 + 3]);
    b_a[2 + cxx_tmp] = (double)a * (x[0] * (double)iv0[i8 + 6] + dv10[i8 + 6]);
    b_x[i8] = x[4 + i8] - xg[4 + i8];
  }

  /*  Control cost Hessian & Jacobian */
  for (i8 = 0; i8 < 3; i8++) {
    cx[i8] = ((b_a[i8] * xg[0] + b_a[i8 + 3] * xg[1]) + b_a[i8 + 6] * xg[2]) +
      b_a[i8 + 9] * xg[3];
    cx[i8 + 3] = ((double)iv0[i8] * b_x[0] + (double)iv0[i8 + 3] * b_x[1]) +
      (double)iv0[i8 + 6] * b_x[2];
    cu[i8] = 0.0;
    cu[i8] = (dv0[i8] * u[0] + dv0[i8 + 3] * u[1]) + dv0[i8 + 6] * u[2];
  }
}

/*
 * Arguments    : const double x[7]
 *                const double xg[7]
 *                const double u[3]
 *                double terminal
 * Return Type  : double
 */
static double satellite_cost_efficient(const double x[7], const double xg[7],
  const double u[3], double terminal)
{
  double cost;
  double b[9];
  static const double dv8[9] = { 0.1, 0.0, 0.0, 0.0, 0.1, 0.0, 0.0, 0.0, 0.1 };

  int i4;
  int w;
  double b_xg;
  double cost_tmp[3];
  double d2;
  double d3;

  /*  Copyright (c) 2020 Robotic Exploration Lab */
  /*   */
  /*  Permission is hereby granted, free of charge, to any person obtaining a copy */
  /*  of this software and associated documentation files (the "Software"), to deal */
  /*  in the Software without restriction, including without limitation the rights */
  /*  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell */
  /*  copies of the Software, and to permit persons to whom the Software is */
  /*  furnished to do so, subject to the following conditions: */
  /*   */
  /*  The above copyright notice and this permission notice shall be included in all */
  /*  copies or substantial portions of the Software. */
  /*   */
  /*  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR */
  /*  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, */
  /*  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE */
  /*  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER */
  /*  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, */
  /*  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE */
  /*  SOFTWARE. */
  /*  Calculates the cost contribution of a given state and control */
  /*  Utilizes a standard quadratic cost function for the angular velocity */
  /*  and a geodesic cost for the attitude quaternion */
  /*  Also returns the 2nd order expansion of the cost function */
  /*  Inputs */
  /* ===================================== */
  /*  x        - [quaternion; omega] (7x1) */
  /*  u        - (3x1) (passed in as zeros for final time-step/terminal cost) */
  /*  terminal - integer 0 or 1 */
  /* --------------------------------------------------- */
  /*  control hessian */
  if (terminal != 0.0) {
    for (i4 = 0; i4 < 9; i4++) {
      b[i4] = iv0[i4];
    }

    /*  terminal angular velocity cost hessian */
    w = 10;

    /*  terminal geodesic cost weight */
  } else {
    memcpy(&b[0], &dv8[0], 9U * sizeof(double));

    /*  cumulative angular velocity cost hessian */
    w = 1;

    /*  cumulative geodesic cost weighting */
  }

  /*  Finds the geodesic quaternion-error cost */
  /*  quat_cost = min(1+qd'q, 1-qd'q)) where qd is the desired attitude quaternion  */
  /*  Also records the sign (+ or -) which minimizes quat_cost */
  /*  this is used when calculating the Jacobain and Hessian */
  b_xg = ((xg[0] * x[0] + xg[1] * x[1]) + xg[2] * x[2]) + xg[3] * x[3];
  cost_tmp[0] = x[4] - xg[4];
  cost_tmp[1] = x[5] - xg[5];
  cost_tmp[2] = x[6] - xg[6];
  d2 = 0.0;
  for (i4 = 0; i4 < 3; i4++) {
    d2 += ((0.5 * cost_tmp[0] * b[3 * i4] + 0.5 * cost_tmp[1] * b[1 + 3 * i4]) +
           0.5 * cost_tmp[2] * b[2 + 3 * i4]) * cost_tmp[i4];
  }

  d3 = 0.0;
  for (i4 = 0; i4 < 3; i4++) {
    d3 += ((0.5 * u[0] * dv0[3 * i4] + 0.5 * u[1] * dv0[1 + 3 * i4]) + 0.5 * u[2]
           * dv0[2 + 3 * i4]) * u[i4];
  }

  cost = ((double)w * fmin(1.0 + b_xg, 1.0 - b_xg) + d2) + d3;

  /*  State cost Hessian */
  /*  State cost Jacobian */
  /*  Control cost Hessian & Jacobian */
  return cost;
}

/*
 * Arguments    : const double x0[7]
 *                const double u0[3]
 *                double dt
 *                const double B_ECI[3]
 *                double fx[36]
 *                double fu[18]
 * Return Type  : void
 */
static void satellite_derivatives(const double x0[7], const double u0[3], double
  dt, const double B_ECI[3], double fx[36], double fu[18])
{
  double b[7];
  double dxdot1[70];
  double a;
  int i;
  double b_x0[7];
  double dxdot2[70];
  double dv11[9];
  double x1[7];
  int i9;
  int fx_tmp_tmp;
  double fx_tmp[49];
  double dv12[49];
  double b_fx_tmp[42];
  int i10;
  int i11;
  double dv13[49];
  double c_x0[42];
  double c_fx_tmp[42];
  double b_dt[21];

  /*  Copyright (c) 2020 Robotic Exploration Lab */
  /*   */
  /*  Permission is hereby granted, free of charge, to any person obtaining a copy */
  /*  of this software and associated documentation files (the "Software"), to deal */
  /*  in the Software without restriction, including without limitation the rights */
  /*  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell */
  /*  copies of the Software, and to permit persons to whom the Software is */
  /*  furnished to do so, subject to the following conditions: */
  /*   */
  /*  The above copyright notice and this permission notice shall be included in all */
  /*  copies or substantial portions of the Software. */
  /*   */
  /*  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR */
  /*  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, */
  /*  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE */
  /*  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER */
  /*  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, */
  /*  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE */
  /*  SOFTWARE. */
  /*  Steps the dynamics forward using a 2nd order rk-method */
  /*  Returns: new state, discrete time Jacobians */
  /*  Step dynamics */
  /*  Explicit midpoint step from x_k to x_{k+1} */
  /*  update to rk4 */
  satellite_dynamics(x0, u0, B_ECI, b, dxdot1);
  a = 0.5 * dt;
  for (i = 0; i < 7; i++) {
    b_x0[i] = x0[i] + a * b[i];
  }

  satellite_dynamics(b_x0, u0, B_ECI, b, dxdot2);

  /*  Re-normalize the quaternion */
  /*  Continuous time Jacobians */
  /*  Discrete time Jacobians */
  /*  First form the state error Jacobians E(x_k) & E(x_k+1) */
  /*  Copyright (c) 2020 Robotic Exploration Lab */
  /*   */
  /*  Permission is hereby granted, free of charge, to any person obtaining a copy */
  /*  of this software and associated documentation files (the "Software"), to deal */
  /*  in the Software without restriction, including without limitation the rights */
  /*  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell */
  /*  copies of the Software, and to permit persons to whom the Software is */
  /*  furnished to do so, subject to the following conditions: */
  /*   */
  /*  The above copyright notice and this permission notice shall be included in all */
  /*  copies or substantial portions of the Software. */
  /*   */
  /*  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR */
  /*  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, */
  /*  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE */
  /*  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER */
  /*  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, */
  /*  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE */
  /*  SOFTWARE. */
  /*  Returns skew symmetric - cross product matrix of a 3x1 vector */
  /*  Copyright (c) 2020 Robotic Exploration Lab */
  /*   */
  /*  Permission is hereby granted, free of charge, to any person obtaining a copy */
  /*  of this software and associated documentation files (the "Software"), to deal */
  /*  in the Software without restriction, including without limitation the rights */
  /*  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell */
  /*  copies of the Software, and to permit persons to whom the Software is */
  /*  furnished to do so, subject to the following conditions: */
  /*   */
  /*  The above copyright notice and this permission notice shall be included in all */
  /*  copies or substantial portions of the Software. */
  /*   */
  /*  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR */
  /*  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, */
  /*  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE */
  /*  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER */
  /*  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, */
  /*  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE */
  /*  SOFTWARE. */
  /*  Returns skew symmetric - cross product matrix of a 3x1 vector */
  a = 0.5 * dt * dt;
  for (i = 0; i < 7; i++) {
    x1[i] = x0[i] + dt * b[i];
    for (i9 = 0; i9 < 7; i9++) {
      fx_tmp_tmp = i9 + 7 * i;
      fx_tmp[fx_tmp_tmp] = a * dxdot2[fx_tmp_tmp];
    }
  }

  dv11[0] = 0.0;
  dv11[3] = -x1[3];
  dv11[6] = x1[2];
  dv11[1] = x1[3];
  dv11[4] = 0.0;
  dv11[7] = -x1[1];
  dv11[2] = -x1[2];
  dv11[5] = x1[1];
  dv11[8] = 0.0;
  for (i9 = 0; i9 < 3; i9++) {
    b_fx_tmp[i9] = -x1[1 + i9];
    b_fx_tmp[i9 + 3] = 0.0;
    fx_tmp_tmp = 6 * (i9 + 1);
    b_fx_tmp[fx_tmp_tmp] = x1[0] * (double)iv0[i9] + dv11[i9];
    b_fx_tmp[fx_tmp_tmp + 3] = 0.0;
    b_fx_tmp[1 + fx_tmp_tmp] = x1[0] * (double)iv0[i9 + 3] + dv11[i9 + 3];
    b_fx_tmp[fx_tmp_tmp + 4] = 0.0;
    b_fx_tmp[2 + fx_tmp_tmp] = x1[0] * (double)iv0[i9 + 6] + dv11[i9 + 6];
    b_fx_tmp[fx_tmp_tmp + 5] = 0.0;
    for (i10 = 0; i10 < 6; i10++) {
      b_fx_tmp[i10 + 6 * (i9 + 4)] = iv1[i9 + 3 * i10];
    }
  }

  eye(dv12);
  for (i9 = 0; i9 < 7; i9++) {
    for (i10 = 0; i10 < 7; i10++) {
      a = 0.0;
      for (i11 = 0; i11 < 7; i11++) {
        a += fx_tmp[i9 + 7 * i11] * dxdot1[i11 + 7 * i10];
      }

      i11 = i9 + 7 * i10;
      dv13[i11] = (dv12[i11] + dt * dxdot2[i11]) + a;
    }
  }

  dv11[0] = 0.0;
  dv11[3] = -x0[3];
  dv11[6] = x0[2];
  dv11[1] = x0[3];
  dv11[4] = 0.0;
  dv11[7] = -x0[1];
  dv11[2] = -x0[2];
  dv11[5] = x0[1];
  dv11[8] = 0.0;
  for (i9 = 0; i9 < 6; i9++) {
    for (i10 = 0; i10 < 7; i10++) {
      fx_tmp_tmp = i9 + 6 * i10;
      c_fx_tmp[fx_tmp_tmp] = 0.0;
      a = 0.0;
      for (i11 = 0; i11 < 7; i11++) {
        a += b_fx_tmp[i9 + 6 * i11] * dv13[i11 + 7 * i10];
      }

      c_fx_tmp[fx_tmp_tmp] = a;
    }
  }

  for (i9 = 0; i9 < 3; i9++) {
    c_x0[7 * i9] = -x0[1 + i9];
    i = 7 * (i9 + 3);
    c_x0[i] = 0.0;
    c_x0[7 * i9 + 1] = x0[0] * (double)iv0[3 * i9] + dv11[3 * i9];
    c_x0[i + 1] = 0.0;
    fx_tmp_tmp = 1 + 3 * i9;
    c_x0[7 * i9 + 2] = x0[0] * (double)iv0[fx_tmp_tmp] + dv11[fx_tmp_tmp];
    c_x0[i + 2] = 0.0;
    fx_tmp_tmp = 2 + 3 * i9;
    c_x0[7 * i9 + 3] = x0[0] * (double)iv0[fx_tmp_tmp] + dv11[fx_tmp_tmp];
    c_x0[i + 3] = 0.0;
  }

  for (i9 = 0; i9 < 6; i9++) {
    c_x0[7 * i9 + 4] = iv1[3 * i9];
    c_x0[7 * i9 + 5] = iv1[1 + 3 * i9];
    c_x0[7 * i9 + 6] = iv1[2 + 3 * i9];
  }

  for (i9 = 0; i9 < 6; i9++) {
    for (i10 = 0; i10 < 6; i10++) {
      i = i9 + 6 * i10;
      fx[i] = 0.0;
      a = 0.0;
      for (i11 = 0; i11 < 7; i11++) {
        a += c_fx_tmp[i9 + 6 * i11] * c_x0[i11 + 7 * i10];
      }

      fx[i] = a;
    }
  }

  for (i9 = 0; i9 < 7; i9++) {
    for (i10 = 0; i10 < 3; i10++) {
      a = 0.0;
      for (i11 = 0; i11 < 7; i11++) {
        a += fx_tmp[i9 + 7 * i11] * dxdot1[i11 + 7 * (7 + i10)];
      }

      b_dt[i9 + 7 * i10] = dt * dxdot2[i9 + 7 * (7 + i10)] + a;
    }
  }

  for (i9 = 0; i9 < 6; i9++) {
    for (i10 = 0; i10 < 3; i10++) {
      i = i9 + 6 * i10;
      fu[i] = 0.0;
      a = 0.0;
      for (i11 = 0; i11 < 7; i11++) {
        a += b_fx_tmp[i9 + 6 * i11] * b_dt[i11 + 7 * i10];
      }

      fu[i] = a;
    }
  }
}

/*
 * Arguments    : const double x[7]
 *                const double u[3]
 *                const double B_ECI[3]
 *                double xdot[7]
 *                double dxdot[70]
 * Return Type  : void
 */
static void satellite_dynamics(const double x[7], const double u[3], const
  double B_ECI[3], double xdot[7], double dxdot[70])
{
  double dv4[4];
  double dv5[4];
  double b_x[4];
  double output[4];
  double qdot_tmp[9];
  int i2;
  double wdot_tmp[9];
  double b_output[3];
  int i3;
  double dv6[12];
  double b_wdot_tmp[3];
  int a_tmp;
  static const signed char a[9] = { 100, 0, 0, 0, 100, 0, 0, 0, 100 };

  int wdot_tmp_tmp;
  double c_wdot_tmp[9];
  double dxdot_tmp;
  static const double b[9] = { 0.01, 0.0, 0.0, 0.0, 0.01, 0.0, 0.0, 0.0, 0.01 };

  double dv7[9];
  double b_a[9];
  static const short c_a[9] = { -200, 0, 0, 0, -200, 0, 0, 0, -200 };

  int b_dxdot_tmp;
  int c_dxdot_tmp;
  static const signed char d_a[9] = { -100, 0, 0, 0, -100, 0, 0, 0, -100 };

  /*  Calculates the continuous time state derivative and Jacobians */
  /*  Inputs */
  /* ======================= */
  /*  x - the current state */
  /*  u - the current control (magentic moment vector) */
  /*  B - Earth magnetic field vector in ECI co-ordinates (3x1) */
  /*  TODO: ideally, we'd like to pass in a B_B (body frame magnetic field */
  /*  vector) */
  /*  kgm^2 */
  /*  Angular velocity */
  /*  Quaternion */
  /*  Magnetic field section  */
  /*  given B_ECI (ECI magnetic field at the time step) */
  /*  this function multiplie s a vector by a quaternion rotation */
  dv4[0] = 0.0;
  dv4[1] = B_ECI[0];
  dv4[2] = B_ECI[1];
  dv4[3] = B_ECI[2];
  qmult(dv4, *(double (*)[4])&x[0], dv5);
  b_x[0] = x[0];
  b_x[1] = -x[1];
  b_x[2] = -x[2];
  b_x[3] = -x[3];
  qmult(b_x, dv5, output);
  dv4[0] = 0.0;
  dv4[1] = B_ECI[0];
  dv4[2] = B_ECI[1];
  dv4[3] = B_ECI[2];
  qmult(dv4, *(double (*)[4])&x[0], dv5);
  b_x[0] = x[0];
  b_x[1] = -x[1];
  b_x[2] = -x[2];
  b_x[3] = -x[3];
  qmult(b_x, dv5, dv4);

  /*  Non-linear dynamics */
  /*  Copyright (c) 2020 Robotic Exploration Lab */
  /*   */
  /*  Permission is hereby granted, free of charge, to any person obtaining a copy */
  /*  of this software and associated documentation files (the "Software"), to deal */
  /*  in the Software without restriction, including without limitation the rights */
  /*  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell */
  /*  copies of the Software, and to permit persons to whom the Software is */
  /*  furnished to do so, subject to the following conditions: */
  /*   */
  /*  The above copyright notice and this permission notice shall be included in all */
  /*  copies or substantial portions of the Software. */
  /*   */
  /*  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR */
  /*  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, */
  /*  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE */
  /*  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER */
  /*  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, */
  /*  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE */
  /*  SOFTWARE. */
  /*  Returns skew symmetric - cross product matrix of a 3x1 vector */
  qdot_tmp[0] = 0.0;
  qdot_tmp[3] = -x[3];
  qdot_tmp[6] = x[2];
  qdot_tmp[1] = x[3];
  qdot_tmp[4] = 0.0;
  qdot_tmp[7] = -x[1];
  qdot_tmp[2] = -x[2];
  qdot_tmp[5] = x[1];
  qdot_tmp[8] = 0.0;
  for (i2 = 0; i2 < 9; i2++) {
    qdot_tmp[i2] += x[0] * (double)iv0[i2];
  }

  /*  wdot = Jinv*(u - skew_mat(w)*J*w); */
  /*  Copyright (c) 2020 Robotic Exploration Lab */
  /*   */
  /*  Permission is hereby granted, free of charge, to any person obtaining a copy */
  /*  of this software and associated documentation files (the "Software"), to deal */
  /*  in the Software without restriction, including without limitation the rights */
  /*  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell */
  /*  copies of the Software, and to permit persons to whom the Software is */
  /*  furnished to do so, subject to the following conditions: */
  /*   */
  /*  The above copyright notice and this permission notice shall be included in all */
  /*  copies or substantial portions of the Software. */
  /*   */
  /*  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR */
  /*  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, */
  /*  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE */
  /*  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER */
  /*  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, */
  /*  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE */
  /*  SOFTWARE. */
  /*  Returns skew symmetric - cross product matrix of a 3x1 vector */
  wdot_tmp[0] = 0.0;
  wdot_tmp[3] = -x[6];
  wdot_tmp[6] = x[5];
  wdot_tmp[1] = x[6];
  wdot_tmp[4] = 0.0;
  wdot_tmp[7] = -x[4];
  wdot_tmp[2] = -x[5];
  wdot_tmp[5] = x[4];
  wdot_tmp[8] = 0.0;

  /*  with magnetorquers */
  b_output[0] = -(output[2] * u[2] - output[3] * u[1]);
  b_output[1] = -(output[3] * u[0] - output[1] * u[2]);
  b_output[2] = -(output[1] * u[1] - output[2] * u[0]);
  for (i2 = 0; i2 < 3; i2++) {
    i3 = i2 << 2;
    dv6[i3] = 0.5 * -x[1 + i2];
    b_wdot_tmp[i2] = 0.0;
    for (a_tmp = 0; a_tmp < 3; a_tmp++) {
      wdot_tmp_tmp = i2 + 3 * a_tmp;
      c_wdot_tmp[wdot_tmp_tmp] = 0.0;
      dxdot_tmp = (wdot_tmp[i2] * b[3 * a_tmp] + wdot_tmp[i2 + 3] * b[1 + 3 *
                   a_tmp]) + wdot_tmp[i2 + 6] * b[2 + 3 * a_tmp];
      c_wdot_tmp[wdot_tmp_tmp] = dxdot_tmp;
      dv6[(a_tmp + i3) + 1] = 0.5 * qdot_tmp[a_tmp + 3 * i2];
      b_wdot_tmp[i2] += dxdot_tmp * x[4 + a_tmp];
    }

    b_output[i2] -= b_wdot_tmp[i2];
  }

  for (i2 = 0; i2 < 4; i2++) {
    dv5[i2] = 0.0;
    dv5[i2] = (dv6[i2] * x[4] + dv6[i2 + 4] * x[5]) + dv6[i2 + 8] * x[6];
  }

  for (i2 = 0; i2 < 3; i2++) {
    b_wdot_tmp[i2] = 0.0;
    b_wdot_tmp[i2] = ((double)a[i2] * b_output[0] + (double)a[i2 + 3] *
                      b_output[1]) + (double)a[i2 + 6] * b_output[2];
  }

  xdot[0] = dv5[0];
  xdot[1] = dv5[1];
  xdot[2] = dv5[2];
  xdot[3] = dv5[3];

  /*  Jacobians  */
  for (i2 = 0; i2 < 3; i2++) {
    xdot[i2 + 4] = b_wdot_tmp[i2];
    b_wdot_tmp[i2] = (b[i2] * x[4] + b[i2 + 3] * x[5]) + b[i2 + 6] * x[6];
  }

  /*  Copyright (c) 2020 Robotic Exploration Lab */
  /*   */
  /*  Permission is hereby granted, free of charge, to any person obtaining a copy */
  /*  of this software and associated documentation files (the "Software"), to deal */
  /*  in the Software without restriction, including without limitation the rights */
  /*  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell */
  /*  copies of the Software, and to permit persons to whom the Software is */
  /*  furnished to do so, subject to the following conditions: */
  /*   */
  /*  The above copyright notice and this permission notice shall be included in all */
  /*  copies or substantial portions of the Software. */
  /*   */
  /*  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR */
  /*  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, */
  /*  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE */
  /*  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER */
  /*  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, */
  /*  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE */
  /*  SOFTWARE. */
  /*  Returns skew symmetric - cross product matrix of a 3x1 vector */
  /*  B = [zeros(4,3); */
  /*       Jinv]; */
  /*  with magnetorquers: */
  /*  Copyright (c) 2020 Robotic Exploration Lab */
  /*   */
  /*  Permission is hereby granted, free of charge, to any person obtaining a copy */
  /*  of this software and associated documentation files (the "Software"), to deal */
  /*  in the Software without restriction, including without limitation the rights */
  /*  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell */
  /*  copies of the Software, and to permit persons to whom the Software is */
  /*  furnished to do so, subject to the following conditions: */
  /*   */
  /*  The above copyright notice and this permission notice shall be included in all */
  /*  copies or substantial portions of the Software. */
  /*   */
  /*  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR */
  /*  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, */
  /*  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE */
  /*  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER */
  /*  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, */
  /*  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE */
  /*  SOFTWARE. */
  /*  Returns skew symmetric - cross product matrix of a 3x1 vector */
  dv7[0] = 0.0;
  dv7[3] = -b_wdot_tmp[2];
  dv7[6] = b_wdot_tmp[1];
  dv7[1] = b_wdot_tmp[2];
  dv7[4] = 0.0;
  dv7[7] = -b_wdot_tmp[0];
  dv7[2] = -b_wdot_tmp[1];
  dv7[5] = b_wdot_tmp[0];
  dv7[8] = 0.0;
  for (i2 = 0; i2 < 9; i2++) {
    c_wdot_tmp[i2] -= dv7[i2];
  }

  for (i2 = 0; i2 < 3; i2++) {
    for (i3 = 0; i3 < 3; i3++) {
      a_tmp = i2 + 3 * i3;
      b_a[a_tmp] = 0.0;
      b_a[a_tmp] = ((double)c_a[i2] * c_wdot_tmp[3 * i3] + (double)c_a[i2 + 3] *
                    c_wdot_tmp[1 + 3 * i3]) + (double)c_a[i2 + 6] * c_wdot_tmp[2
        + 3 * i3];
    }
  }

  dv7[0] = 0.0;
  dv7[3] = -dv4[3];
  dv7[6] = dv4[2];
  dv7[1] = dv4[3];
  dv7[4] = 0.0;
  dv7[7] = -dv4[1];
  dv7[2] = -dv4[2];
  dv7[5] = dv4[1];
  dv7[8] = 0.0;
  dxdot[0] = 0.0;
  for (i2 = 0; i2 < 3; i2++) {
    dxdot_tmp = x[4 + i2];
    b_dxdot_tmp = 7 * (i2 + 1);
    dxdot[b_dxdot_tmp] = 0.5 * -dxdot_tmp;
    c_dxdot_tmp = 7 * (i2 + 4);
    dxdot[c_dxdot_tmp] = 0.5 * -x[1 + i2];
    dxdot[i2 + 1] = 0.5 * dxdot_tmp;
    for (i3 = 0; i3 < 3; i3++) {
      wdot_tmp_tmp = i2 + 3 * i3;
      c_wdot_tmp[wdot_tmp_tmp] = 0.0;
      c_wdot_tmp[wdot_tmp_tmp] = ((double)d_a[i2] * dv7[3 * i3] + (double)d_a[i2
        + 3] * dv7[1 + 3 * i3]) + (double)d_a[i2 + 6] * dv7[2 + 3 * i3];
      a_tmp = i3 + 3 * i2;
      dxdot[(i3 + b_dxdot_tmp) + 1] = 0.5 * -wdot_tmp[a_tmp];
      dxdot[(i3 + c_dxdot_tmp) + 1] = 0.5 * qdot_tmp[a_tmp];
    }
  }

  for (i2 = 0; i2 < 4; i2++) {
    dxdot[7 * i2 + 4] = 0.0;
    dxdot[7 * i2 + 5] = 0.0;
    dxdot[7 * i2 + 6] = 0.0;
  }

  for (i2 = 0; i2 < 3; i2++) {
    b_dxdot_tmp = 7 * (i2 + 4);
    dxdot[b_dxdot_tmp + 4] = 0.5 * b_a[3 * i2];
    c_dxdot_tmp = 1 + 3 * i2;
    dxdot[b_dxdot_tmp + 5] = 0.5 * b_a[c_dxdot_tmp];
    a_tmp = 2 + 3 * i2;
    dxdot[b_dxdot_tmp + 6] = 0.5 * b_a[a_tmp];
    b_dxdot_tmp = 7 * (i2 + 7);
    dxdot[b_dxdot_tmp] = 0.0;
    dxdot[1 + b_dxdot_tmp] = 0.0;
    dxdot[2 + b_dxdot_tmp] = 0.0;
    dxdot[3 + b_dxdot_tmp] = 0.0;
    dxdot[b_dxdot_tmp + 4] = c_wdot_tmp[3 * i2];
    dxdot[b_dxdot_tmp + 5] = c_wdot_tmp[c_dxdot_tmp];
    dxdot[b_dxdot_tmp + 6] = c_wdot_tmp[a_tmp];
  }
}

/*
 * Arguments    : double *lambda
 *                double direction
 * Return Type  : void
 */
static void updateLambda(double *lambda, double direction)
{
  /*  Increases or decreases the regularization parameter according */
  /*  to a non-linear scaling regime. */
  if (direction == 1.0) {
    /*  increase lambda */
    *lambda = fmax(*lambda * 1.6, 1.0E-6);
  } else {
    /*  decrease lambda */
    *lambda = *lambda * 0.625 * (double)(*lambda > 1.0E-6);

    /*  set = 0 if lambda too small */
  }
}

/*
 * Arguments    : const double x0[7]
 *                const double xg[7]
 *                const double u0[297]
 *                const double u_lims[6]
 *                double dt
 *                const double B_ECI[297]
 *                double x[700]
 *                double u[297]
 *                double K[1782]
 *                bool *result
 * Return Type  : void
 */
void milqr_efficient(const double x0[7], const double xg[7], const double u0[297],
                     const double u_lims[6], double dt, const double B_ECI[297],
                     double x[700], double u[297], double K[1782], bool *r)
{
  double alphas[11];
  double lambda;
  double cost_n;
  int i0;
  int j;
  static double x_n[700];
  static double u_n[297];
  double dV[2];
  static double l[297];
  double dv1[297];
  static double dv2[1782];
  static double b_x0[700];
  double cost;
  int iter;
  int b_iter;
  bool exitg1;
  bool backPassCheck;
  int exitg2;
  bool diverged;
  double varargin_1[297];
  double maxval[99];
  double d0;
  bool lineSearchCheck;
  bool exitg3;
  double expected_change;
  double dcost;
  double c_ratio;

  /*  Solves finite horizon optimal control problem using a */
  /*  multiplicative iterative linear quadratic regulator */
  /*  Note that this solver will not work without control limits */
  /*  Inputs */
  /*  =========================================== */
  /*  x0 - The intial state (n, 1) */
  /*  */
  /*  xg - The goal state (n, 1) */
  /*  */
  /*  u0 - The initial control sequeunce (m, N-1) */
  /*  */
  /*  u_lims - The control limits (m, 2) (lower, upper) */
  /*  */
  /*  B_ECI - A time sequence of Earth magnetic field vectors in ECI (3, N) */
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
  /*  dt = 0.03;            % Timestep (Should be highest we can get away with) */
  /*  maximum iterations */
  /*  cost reduction exit tolerance */
  /*  feedforward control change exit criterion */
  /*  max regularizaion param allowed for exit */
  /*  minimum accepted cost reduction ratio */
  /*  maximum regularization parameter */
  /*  set lambda = 0 below this value */
  /*  amount to scale dlambda by */
  /*  Init optimisation params */
  power(alphas);

  /*  line search param vector */
  lambda = 1.0;

  /*  error state size (3 param. error representation for attitude) */
  /*  variable initialization */
  cost_n = 0.0;

  /*  initialize x and u */
  /*  initialize the "next" variables */
  for (i0 = 0; i0 < 100; i0++) {
    for (j = 0; j < 7; j++) {
      x_n[j + 7 * i0] = x0[j];
    }
  }

  /*  Initialize Forward rollout */
  for (i0 = 0; i0 < 297; i0++) {
    u_n[i0] = u0[i0];
    l[i0] = 0.0;
  }

  memset(&K[0], 0, 1782U * sizeof(double));
  dV[0] = 0.0;
  dV[1] = 0.0;

  /*  Throw out everything after x and u, then in backward pass evalute on the */
  /*  spot */
  for (i0 = 0; i0 < 100; i0++) {
    for (j = 0; j < 7; j++) {
      b_x0[j + 7 * i0] = x0[j];
    }
  }

  memset(&dv1[0], 0, 297U * sizeof(double));
  memset(&dv2[0], 0, 1782U * sizeof(double));
  forwardRollout(b_x0, xg, u0, dv1, dv2, 0.0, u_lims, dt, B_ECI, x, u, &cost);

  /*  Convergence check params */
  /*  Expected cost change */
  /*  Ratio of cost change to expected cost change */
  // printf("%d\n",*r);
  *r = false;

  printf("\n=====Running MILQR Optimisation====\n");
  iter = 1;
  b_iter = 0;
  exitg1 = false;
  while ((!exitg1) && (b_iter < 300)) {
    iter = b_iter + 1;

    printf("\n---New Iteration---");
    printf("\n lambda: %f\n",lambda);
    /*      fprintf(string(lambda)); */
    /*      if exist('dcost') */
    /*          fprintf("\n cost change: "); */
    /*          fprintf(string(dcost)); */
    /*      end */
    /*      fprintf("\n"); */
    /*  Backward Pass */
    /* ======================================= */
    backPassCheck = false;
    do {
      exitg2 = 0;
      if (!backPassCheck) {
        /*  don't pass in gradient info, calculate in place */
        backwardPass(lambda, u_lims, u, x, xg, dt, B_ECI, l, K, dV, &diverged);
        //printf("%d\n",diverged);
        if (diverged) {
          printf("---Warning: Cholesky factorizaton failed---\n");
          /*  Increase regularization parameter (lambda) */
          updateLambda(&lambda, 1.0);
          if (lambda > 1.0E+7) {
            exitg2 = 1;
          }
        } else {
          backPassCheck = true;
        }
      } else {
        exitg2 = 1;
      }
    } while (exitg2 == 0);

    /*  Check relative change of feedforward control */
    /*  Terminate if sufficiently small (success) */
    b_abs(l, varargin_1);
    b_abs(u, dv1);
    for (i0 = 0; i0 < 297; i0++) {
      varargin_1[i0] /= dv1[i0] + 1.0;
    }

    for (j = 0; j < 99; j++) {
      maxval[j] = varargin_1[3 * j];
      d0 = varargin_1[3 * j + 1];
      if (maxval[j] < d0) {
        maxval[j] = d0;
      }

      d0 = varargin_1[3 * j + 2];
      if (maxval[j] < d0) {
        maxval[j] = d0;
      }
    }

    /*  Avg over time of max change */
    if ((lambda < 1.0E-5) && (mean(maxval) < 0.0001)) {
      printf("\n---Success: Control change decreased below tolerance---\n");
      *r = true;
      exitg1 = true;
    } else {
      /*  Forward Line-Search */
      /* =========================================== */
      lineSearchCheck = false;
      if (backPassCheck) {
        //printf("hello");
        j = 0;
        exitg3 = false;
        while ((!exitg3) && (j < 11)) {
          forwardRollout(x, xg, u, l, K, alphas[j], u_lims, dt, B_ECI, x_n, u_n,
                         &cost_n);
          expected_change = alphas[j] * dV[0] + pow(alphas[j], 2.0) * dV[1];
          if (expected_change < 0.0) {
            c_ratio = (cost_n - cost) / expected_change;
          } else {
            /*  Non-positive expected cost reduction */
            /*  actual cost change must be negative to accept the step */
            d0 = cost_n - cost;
            b_sign(&d0);
            c_ratio = -d0;
          }

          //printf("c_ratio: %f\n",c_ratio);

          if (c_ratio > 0.0) {
            printf("Cost ratio is good ... moving to line search.\n");
            lineSearchCheck = true;
            exitg3 = true;
          } else {
            j++;
          }
        }
      }

      /*  Parameter Updates */
      /* ============================================= */
      if (lineSearchCheck) {
        /*  Decrease Lambda */
        updateLambda(&lambda, -1.0);
        dcost = cost - cost_n;

        /*  Update the trajectory and controls */
        memcpy(&x[0], &x_n[0], 700U * sizeof(double));
        memcpy(&u[0], &u_n[0], 297U * sizeof(double));
        cost = cost_n;

        /*  Change in cost small enough to terminate? */
        if (dcost < 1.0E-7) {
          *r = true;

          printf("\n---Success: cost change < tolerance---\n");
          exitg1 = true;
        } else {
          b_iter++;
        }
      } else {
        /*  No cost reduction (based on cost change ratio) */
        /*  Increase lambda */
        printf("No cost reduction. \n");
        updateLambda(&lambda, 1.0);
        if (lambda > 1.0E+7) {
          printf("\n---Diverged: new lambda > lambda_max---\n");
          exitg1 = true;
        } else {
          b_iter++;
        }
      }
    }
  }

  if (iter == 300) {
    /*  Ddin't converge completely */
    *r = false;

    /*      % fprintf("\n---Warning: Max iterations exceeded---\n"); */
  }
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void milqr_efficient_initialize(void)
{
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void milqr_efficient_terminate(void)
{
  /* (no terminate code required) */
}

/*
 * File trailer for milqr_efficient.c
 *
 * [EOF]
 */
