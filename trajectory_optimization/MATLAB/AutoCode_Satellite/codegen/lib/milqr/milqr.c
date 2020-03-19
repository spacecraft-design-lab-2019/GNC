/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: milqr.c
 *
 * MATLAB Coder version            : 4.3
 * C/C++ source code generated on  : 17-Mar-2020 19:06:03
 */

/* Include Files */
#include "milqr.h"
#include "rt_nonfinite.h"
#include <math.h>
#include <string.h>

/* Variable Definitions */
static const signed char iv[9] = { 1, 0, 0, 0, 1, 0, 0, 0, 1 };

static const double dv[9] = { 0.01, 0.0, 0.0, 0.0, 0.01, 0.0, 0.0, 0.0, 0.01 };

static const double dv1[9] = { 0.05, 0.0, 0.0, 0.0, 0.05, 0.0, 0.0, 0.0, 0.05 };

/* Function Declarations */
static void b_mldivide(const double A_data[], const int A_size[2], double
  B_data[], int B_size[2]);
static void backwardPass(const double fx[17964], const double fu[8982], const
  double cx[3000], const double cu[1497], const double cxx[18000], const double
  cuu[4491], double lambda, const double u_lims[6], const double u[1497], double
  l[1497], double K[8982], double dV[2], boolean_T *diverge);
static void boxQPsolve(const double Quu[9], const double Qu[3], const double
  lower[3], const double upper[3], const double u0[3], double u[3], double
  *result, double Luu_data[], int Luu_size[2], boolean_T b_free[3]);
static void forwardRollout(const double x[3500], const double xg[7], const
  double u[1497], const double l[1497], const double K[8982], double alpha,
  const double u_lims[6], double xnew[3500], double unew[1497], double fx[17964],
  double fu[8982], double cx[3000], double cu[1497], double cxx[18000], double
  cuu[4491], double *cost);
static void mldivide(const double A_data[], const int A_size[2], double B_data[],
                     int B_size[1]);
static void qrpf(double A_data[], const int A_size[2], int m, int n, double
                 tau_data[], int jpvt_data[]);
static int rankFromQR(const double A_data[], const int A_size[2]);
static double rt_hypotd_snf(double u0, double u1);
static double rt_powd_snf(double u0, double u1);
static void satellite_cost(const double x[7], const double xg[7], const double
  u[3], double terminal, double *cost, double cx[6], double cu[3], double cxx[36]);
static void xgeqp3(double A_data[], const int A_size[2], double tau_data[], int
                   tau_size[1], int jpvt_data[], int jpvt_size[2]);
static void xgetrf(int m, int n, double A_data[], const int A_size[2], int lda,
                   int ipiv_data[], int ipiv_size[2], int *info);
static double xnrm2(int n, const double x_data[], int ix0);

/* Function Definitions */

/*
 * Arguments    : const double A_data[]
 *                const int A_size[2]
 *                double B_data[]
 *                int B_size[2]
 * Return Type  : void
 */
static void b_mldivide(const double A_data[], const int A_size[2], double
  B_data[], int B_size[2])
{
  int b_A_size[2];
  int mn;
  int jBcol;
  int i;
  double b_A_data[9];
  int i1;
  double tau_data[3];
  int tau_size[1];
  int jpvt_data[3];
  int jpvt_size[2];
  int rankA;
  int B_size_idx_0;
  int m;
  double b_B_data[18];
  int j;
  int xj;
  int k;
  int wj_tmp;
  double wj;
  int b_i;
  if ((A_size[0] == 0) || (A_size[1] == 0) || (B_size[0] == 0)) {
    B_size[0] = (signed char)A_size[1];
    B_size[1] = 6;
    jBcol = (signed char)A_size[1];
    for (i = 0; i < 6; i++) {
      for (i1 = 0; i1 < jBcol; i1++) {
        B_data[i1 + B_size[0] * i] = 0.0;
      }
    }
  } else if (A_size[0] == A_size[1]) {
    mn = A_size[1];
    b_A_size[0] = A_size[0];
    b_A_size[1] = A_size[1];
    jBcol = A_size[0] * A_size[1];
    if (0 <= jBcol - 1) {
      memcpy(&b_A_data[0], &A_data[0], jBcol * sizeof(double));
    }

    xgetrf(A_size[1], A_size[1], b_A_data, b_A_size, A_size[1], jpvt_data,
           jpvt_size, &jBcol);
    i = A_size[1];
    for (m = 0; m <= i - 2; m++) {
      if (jpvt_data[m] != m + 1) {
        for (xj = 0; xj < 6; xj++) {
          jBcol = B_size[0] * xj;
          wj_tmp = m + jBcol;
          wj = B_data[wj_tmp];
          i1 = (jpvt_data[m] + jBcol) - 1;
          B_data[wj_tmp] = B_data[i1];
          B_data[i1] = wj;
        }
      }
    }

    if (B_size[0] != 0) {
      for (j = 0; j < 6; j++) {
        jBcol = mn * j;
        for (k = 0; k < mn; k++) {
          m = mn * k;
          i = k + jBcol;
          if (B_data[i] != 0.0) {
            i1 = k + 2;
            for (b_i = i1; b_i <= mn; b_i++) {
              xj = (b_i + jBcol) - 1;
              B_data[xj] -= B_data[i] * b_A_data[(b_i + m) - 1];
            }
          }
        }
      }
    }

    if (B_size[0] != 0) {
      for (j = 0; j < 6; j++) {
        jBcol = mn * j - 1;
        for (k = mn; k >= 1; k--) {
          m = mn * (k - 1) - 1;
          i = k + jBcol;
          if (B_data[i] != 0.0) {
            B_data[i] /= b_A_data[k + m];
            for (b_i = 0; b_i <= k - 2; b_i++) {
              i1 = (b_i + jBcol) + 1;
              B_data[i1] -= B_data[i] * b_A_data[(b_i + m) + 1];
            }
          }
        }
      }
    }
  } else {
    b_A_size[0] = A_size[0];
    b_A_size[1] = A_size[1];
    jBcol = A_size[0] * A_size[1];
    if (0 <= jBcol - 1) {
      memcpy(&b_A_data[0], &A_data[0], jBcol * sizeof(double));
    }

    xgeqp3(b_A_data, b_A_size, tau_data, tau_size, jpvt_data, jpvt_size);
    rankA = rankFromQR(b_A_data, b_A_size);
    B_size_idx_0 = B_size[0];
    jBcol = B_size[0] * B_size[1];
    if (0 <= jBcol - 1) {
      memcpy(&b_B_data[0], &B_data[0], jBcol * sizeof(double));
    }

    B_size[0] = (signed char)b_A_size[1];
    B_size[1] = 6;
    jBcol = (signed char)b_A_size[1];
    for (i = 0; i < 6; i++) {
      for (i1 = 0; i1 < jBcol; i1++) {
        B_data[i1 + B_size[0] * i] = 0.0;
      }
    }

    m = b_A_size[0];
    jBcol = b_A_size[0];
    mn = b_A_size[1];
    if (jBcol < mn) {
      mn = jBcol;
    }

    for (j = 0; j < mn; j++) {
      if (tau_data[j] != 0.0) {
        i = j + 2;
        for (k = 0; k < 6; k++) {
          wj_tmp = B_size_idx_0 * k;
          xj = j + wj_tmp;
          wj = b_B_data[xj];
          for (b_i = i; b_i <= m; b_i++) {
            wj += b_A_data[(b_i + b_A_size[0] * j) - 1] * b_B_data[(b_i + wj_tmp)
              - 1];
          }

          wj *= tau_data[j];
          if (wj != 0.0) {
            b_B_data[xj] -= wj;
            i1 = j + 2;
            for (b_i = i1; b_i <= m; b_i++) {
              xj = (b_i + wj_tmp) - 1;
              b_B_data[xj] -= b_A_data[(b_i + b_A_size[0] * j) - 1] * wj;
            }
          }
        }
      }
    }

    for (k = 0; k < 6; k++) {
      for (b_i = 0; b_i < rankA; b_i++) {
        B_data[(jpvt_data[b_i] + B_size[0] * k) - 1] = b_B_data[b_i +
          B_size_idx_0 * k];
      }

      for (j = rankA; j >= 1; j--) {
        i = B_size[0] * k;
        i1 = (jpvt_data[j - 1] + i) - 1;
        xj = b_A_size[0] * (j - 1);
        B_data[i1] /= b_A_data[(j + xj) - 1];
        for (b_i = 0; b_i <= j - 2; b_i++) {
          i1 = (jpvt_data[b_i] + i) - 1;
          B_data[i1] -= B_data[(jpvt_data[j - 1] + B_size[0] * k) - 1] *
            b_A_data[b_i + xj];
        }
      }
    }
  }
}

/*
 * Perfoms the LQR backward pass to find the optimal controls
 *  Solves a quadratic program (QP) at each timestep for the optimal
 *  controls given the control limits
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
  double y;
  int i2;
  double Qu[3];
  double dV_idx_0;
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
  boolean_T b_y;
  boolean_T exitg2;
  int trueCount;
  signed char tmp_data[3];
  double b_cx[6];
  double Vx_tmp[6];
  double b_cxx[36];
  double b_fx[36];
  double b_Vx_tmp[18];
  signed char b_tmp_data[3];
  int tmp_size[2];
  double c_cxx[36];
  double d;
  int i3;
  int b_Luu_size[2];
  int i4;
  double c_tmp_data[18];

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
      y = 0.0;
      for (i2 = 0; i2 < 6; i2++) {
        y += fu[(i2 + 6 * i1) + i] * Vx[i2];
        dV_idx_0 = 0.0;
        for (b_k = 0; b_k < 6; b_k++) {
          dV_idx_0 += fu[(b_k + 6 * i1) + i] * Vxx[b_k + 6 * i2];
        }

        b_fu[i1 + 3 * i2] = dV_idx_0;
      }

      Qu[i1] = cu[i1 + 3 * (498 - k)] + y;
      for (i2 = 0; i2 < 3; i2++) {
        y = 0.0;
        for (b_k = 0; b_k < 6; b_k++) {
          y += b_fu[i1 + 3 * b_k] * fu[(b_k + 6 * i2) + 18 * (498 - k)];
        }

        b_k = i1 + 3 * i2;
        b_Quu[b_k] = cuu[b_k + 9 * (498 - k)] + y;
      }
    }

    for (i = 0; i < 3; i++) {
      for (i1 = 0; i1 < 6; i1++) {
        y = 0.0;
        for (i2 = 0; i2 < 6; i2++) {
          y += fu[(i2 + 6 * i) + 18 * (498 - k)] * Vxx[i2 + 6 * i1];
        }

        b_fu[i + 3 * i1] = y;
      }

      for (i1 = 0; i1 < 6; i1++) {
        y = 0.0;
        for (i2 = 0; i2 < 6; i2++) {
          y += b_fu[i + 3 * i2] * fx[(i2 + 6 * i1) + 36 * (498 - k)];
        }

        Qux[i + 3 * i1] = y;
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
      b_y = false;
      b_k = 0;
      exitg2 = false;
      while ((!exitg2) && (b_k < 3)) {
        if (b_free[b_k]) {
          b_y = true;
          exitg2 = true;
        } else {
          b_k++;
        }
      }

      if (b_y) {
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

        tmp_size[0] = trueCount;
        tmp_size[1] = 6;
        for (i = 0; i < 6; i++) {
          for (i1 = 0; i1 < trueCount; i1++) {
            c_tmp_data[i1 + trueCount * i] = Qux[(tmp_data[i1] + 3 * i) - 1];
          }
        }

        b_Luu_size[0] = Luu_size[1];
        b_Luu_size[1] = Luu_size[0];
        b_k = Luu_size[0];
        for (i = 0; i < b_k; i++) {
          trueCount = Luu_size[1];
          for (i1 = 0; i1 < trueCount; i1++) {
            Quu[i1 + b_Luu_size[0] * i] = Luu_data[i + Luu_size[0] * i1];
          }
        }

        b_mldivide(Quu, b_Luu_size, c_tmp_data, tmp_size);
        b_Luu_size[0] = Luu_size[0];
        b_Luu_size[1] = Luu_size[1];
        b_k = Luu_size[0] * Luu_size[1];
        for (i = 0; i < b_k; i++) {
          Quu[i] = -Luu_data[i];
        }

        b_mldivide(Quu, b_Luu_size, c_tmp_data, tmp_size);
        b_k = tmp_size[0];
        for (i = 0; i < 6; i++) {
          for (i1 = 0; i1 < b_k; i1++) {
            Kk[(b_tmp_data[i1] + 3 * i) - 1] = c_tmp_data[i1 + tmp_size[0] * i];
          }
        }
      }

      /*  Update Cost to Go */
      y = 0.0;
      for (i = 0; i < 3; i++) {
        y += ((0.5 * lk[0] * b_Quu[3 * i] + 0.5 * lk[1] * b_Quu[3 * i + 1]) +
              0.5 * lk[2] * b_Quu[3 * i + 2]) * lk[i];
      }

      dV_idx_0 = dV[0] + ((lk[0] * Qu[0] + lk[1] * Qu[1]) + lk[2] * Qu[2]);
      y += dV[1];
      dV[0] = dV_idx_0;
      dV[1] = y;
      for (i = 0; i < 6; i++) {
        y = Kk[3 * i];
        i1 = 3 * i + 1;
        i2 = 3 * i + 2;
        for (b_k = 0; b_k < 3; b_k++) {
          b_Vx_tmp[i + 6 * b_k] = (y * b_Quu[3 * b_k] + Kk[i1] * b_Quu[3 * b_k +
            1]) + Kk[i2] * b_Quu[3 * b_k + 2];
        }

        dV_idx_0 = 0.0;
        for (b_k = 0; b_k < 6; b_k++) {
          dV_idx_0 += fx[(b_k + 6 * i) + 36 * (498 - k)] * Vx[b_k];
        }

        b_cx[i] = ((cx[i + 6 * (498 - k)] + dV_idx_0) + ((b_Vx_tmp[i] * lk[0] +
          b_Vx_tmp[i + 6] * lk[1]) + b_Vx_tmp[i + 12] * lk[2])) + ((y * Qu[0] +
          Kk[i1] * Qu[1]) + Kk[i2] * Qu[2]);
        Vx_tmp[i] = (Qux[3 * i] * lk[0] + Qux[i1] * lk[1]) + Qux[i2] * lk[2];
      }

      for (i = 0; i < 6; i++) {
        Vx[i] = b_cx[i] + Vx_tmp[i];
        for (i1 = 0; i1 < 6; i1++) {
          y = 0.0;
          for (i2 = 0; i2 < 6; i2++) {
            y += fx[(i2 + 6 * i) + 36 * (498 - k)] * Vxx[i2 + 6 * i1];
          }

          b_fx[i + 6 * i1] = y;
        }

        for (i1 = 0; i1 < 6; i1++) {
          y = 0.0;
          for (i2 = 0; i2 < 6; i2++) {
            y += b_fx[i + 6 * i2] * fx[(i2 + 6 * i1) + 36 * (498 - k)];
          }

          trueCount = i + 6 * i1;
          c_cxx[trueCount] = cxx[trueCount + 36 * (498 - k)] + y;
        }

        y = b_Vx_tmp[i + 6];
        dV_idx_0 = b_Vx_tmp[i + 12];
        i1 = 3 * i + 1;
        i2 = 3 * i + 2;
        for (b_k = 0; b_k < 6; b_k++) {
          d = Kk[3 * b_k];
          i3 = 3 * b_k + 1;
          i4 = 3 * b_k + 2;
          trueCount = i + 6 * b_k;
          b_cxx[trueCount] = (c_cxx[trueCount] + ((b_Vx_tmp[i] * d + y * Kk[i3])
            + dV_idx_0 * Kk[i4])) + ((Kk[3 * i] * Qux[3 * b_k] + Kk[i1] * Qux[i3])
            + Kk[i2] * Qux[i4]);
          b_fx[trueCount] = (Qux[3 * i] * d + Qux[i1] * Kk[i3]) + Qux[i2] *
            Kk[i4];
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
 * Finds the optimal control with limits to minimize a quadratic cost
 *  Minimize 0.5*u'*Quu*u + u'*Qu  s.t. lower <= u <= upper
 *
 *   inputs:
 *      Quu       - positive definite matrix              (m * m)
 *      Qu        - bias vector                           (m)
 *      lower     - lower bounds                          (m)
 *      upper     - upper bounds                          (m)
 *      u0        - initial control input for warm-start  (m)
 *
 *   outputs:
 *      u         - solution                   (m)
 *      result    - result type (roughly, higher is better, see below)
 *      Luu       - cholesky factor            (m * m)
 *      free      - set of free dimensions     (m)
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
  double ssq;
  double z1[3];
  double c;
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
  int A_size[2];
  double gnorm;
  int nmj;
  int i1;
  double deltaX_idx_0;
  int info;
  double A_data[9];
  int j;
  double deltaX_idx_1;
  int iac;
  int idxA1j;
  double deltaX_idx_2;
  int idxAjj;
  int ix;
  int iy;
  signed char c_tmp_data[3];
  int ia0;
  double b_A_data[9];
  int tmp_size[1];
  double d_tmp_data[3];
  double sdotg;
  double step;
  double uc_idx_0;
  double uc_idx_1;
  double uc_idx_2;
  double vc;

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
  ssq = fmax(lower[0], fmin(upper[0], u0[0]));
  z1[0] = ssq;
  u[0] = ssq;
  c = ssq * Qu[0];
  ssq = fmax(lower[1], fmin(upper[1], u0[1]));
  z1[1] = ssq;
  u[1] = ssq;
  c += ssq * Qu[1];
  ssq = fmax(lower[2], fmin(upper[2], u0[2]));
  z1[2] = ssq;
  u[2] = ssq;
  c += ssq * Qu[2];
  t = 0.0;
  for (i = 0; i < 3; i++) {
    t += ((0.5 * z1[0] * Quu[3 * i] + 0.5 * z1[1] * Quu[3 * i + 1]) + 0.5 * ssq *
          Quu[3 * i + 2]) * z1[i];
  }

  value = c + t;

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
          ssq = Qu[k] + ((Quu[k] * u[0] + Quu[k + 3] * u[1]) + Quu[k + 6] * u[2]);
          grad[k] = ssq;
          prev_clamped[k] = clamped[k];
          y = false;
          clamped[k] = false;
          if ((u[k] == lower[k]) && (ssq > 0.0)) {
            y = true;
            clamped[k] = true;
          }

          if ((u[k] == upper[k]) && (ssq < 0.0)) {
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

            /*  Wrapper for MATLAB chol for use with auto coder */
            /*  Inputs: */
            /* =========== */
            /*  A - positive semi-definite matrix */
            A_size[0] = trueCount;
            A_size[1] = trueCount;
            for (i = 0; i < trueCount; i++) {
              for (i1 = 0; i1 < trueCount; i1++) {
                A_data[i1 + trueCount * i] = Quu[(tmp_data[i1] + 3 * (tmp_data[i]
                  - 1)) - 1];
              }
            }

            nmj = 0;
            if (trueCount != 0) {
              info = -1;
              j = 0;
              exitg2 = false;
              while ((!exitg2) && (j <= trueCount - 1)) {
                idxA1j = j * trueCount;
                idxAjj = idxA1j + j;
                ssq = 0.0;
                if (j >= 1) {
                  ix = idxA1j;
                  iy = idxA1j;
                  for (k = 0; k < j; k++) {
                    ssq += A_data[ix] * A_data[iy];
                    ix++;
                    iy++;
                  }
                }

                ssq = A_data[idxAjj] - ssq;
                if (ssq > 0.0) {
                  ssq = sqrt(ssq);
                  A_data[idxAjj] = ssq;
                  if (j + 1 < trueCount) {
                    nmj = (trueCount - j) - 2;
                    ia0 = (idxA1j + trueCount) + 1;
                    idxAjj += trueCount;
                    if ((j != 0) && (nmj + 1 != 0)) {
                      iy = idxAjj;
                      i = ia0 + trueCount * nmj;
                      for (iac = ia0; trueCount < 0 ? iac >= i : iac <= i; iac +=
                           trueCount) {
                        ix = idxA1j;
                        c = 0.0;
                        i1 = (iac + j) - 1;
                        for (k = iac; k <= i1; k++) {
                          c += A_data[k - 1] * A_data[ix];
                          ix++;
                        }

                        A_data[iy] += -c;
                        iy += trueCount;
                      }
                    }

                    ssq = 1.0 / ssq;
                    i = (idxAjj + trueCount * nmj) + 1;
                    for (k = idxAjj + 1; trueCount < 0 ? k >= i : k <= i; k +=
                         trueCount) {
                      A_data[k - 1] *= ssq;
                    }
                  }

                  j++;
                } else {
                  A_data[idxAjj] = ssq;
                  info = j;
                  exitg2 = true;
                }
              }

              nmj = info + 1;
              if (info + 1 == 0) {
                idxAjj = trueCount;
              } else {
                idxAjj = info;
              }

              for (j = 0; j < idxAjj; j++) {
                i = j + 2;
                for (k = i; k <= idxAjj; k++) {
                  A_data[(k + trueCount * j) - 1] = 0.0;
                }
              }

              if (1 > idxAjj) {
                iac = 0;
                k = 0;
              } else {
                iac = idxAjj;
                k = idxAjj;
              }

              for (i = 0; i < k; i++) {
                for (i1 = 0; i1 < iac; i1++) {
                  b_A_data[i1 + iac * i] = A_data[i1 + trueCount * i];
                }
              }

              A_size[0] = iac;
              A_size[1] = k;
              iac *= k;
              if (0 <= iac - 1) {
                memcpy(&A_data[0], &b_A_data[0], iac * sizeof(double));
              }
            }

            Luu_size[0] = A_size[0];
            Luu_size[1] = A_size[1];
            iac = A_size[0] * A_size[1];
            if (0 <= iac - 1) {
              memcpy(&Luu_data[0], &A_data[0], iac * sizeof(double));
            }

            if (nmj != 0) {
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
                ssq = 3.3121686421112381E-170;
                for (k = 0; k < trueCount; k++) {
                  c = fabs(grad[b_tmp_data[k] - 1]);
                  if (c > ssq) {
                    t = ssq / c;
                    gnorm = gnorm * t * t + 1.0;
                    ssq = c;
                  } else {
                    t = c / ssq;
                    gnorm += t * t;
                  }
                }

                gnorm = ssq * sqrt(gnorm);
              }
            }

            if (gnorm < 1.0E-8) {
              *result = 4.0;
              exitg1 = true;
            } else {
              /*  get search direction */
              trueCount = 0;
              deltaX_idx_0 = 0.0;
              if (b_free[0]) {
                trueCount = 1;
              }

              deltaX_idx_1 = 0.0;
              if (b_free[1]) {
                trueCount++;
              }

              deltaX_idx_2 = 0.0;
              if (b_free[2]) {
                trueCount++;
              }

              k = 0;
              if (b_free[0]) {
                c_tmp_data[0] = 1;
                k = 1;
              }

              ssq = u[0] * (double)clamped[0];
              if (b_free[1]) {
                c_tmp_data[k] = 2;
                k++;
              }

              c = u[1] * (double)clamped[1];
              if (b_free[2]) {
                c_tmp_data[k] = 3;
              }

              t = u[2] * (double)clamped[2];
              for (i = 0; i < 3; i++) {
                z1[i] = Qu[i] + ((Quu[i] * ssq + Quu[i + 3] * c) + Quu[i + 6] *
                                 t);
              }

              tmp_size[0] = trueCount;
              for (i = 0; i < trueCount; i++) {
                d_tmp_data[i] = z1[c_tmp_data[i] - 1];
              }

              A_size[0] = Luu_size[1];
              A_size[1] = Luu_size[0];
              iac = Luu_size[0];
              for (i = 0; i < iac; i++) {
                k = Luu_size[1];
                for (i1 = 0; i1 < k; i1++) {
                  A_data[i1 + A_size[0] * i] = Luu_data[i + Luu_size[0] * i1];
                }
              }

              mldivide(A_data, A_size, d_tmp_data, tmp_size);
              A_size[0] = Luu_size[0];
              A_size[1] = Luu_size[1];
              iac = Luu_size[0] * Luu_size[1];
              for (i = 0; i < iac; i++) {
                A_data[i] = -Luu_data[i];
              }

              mldivide(A_data, A_size, d_tmp_data, tmp_size);
              k = 0;

              /*  cholesky solver */
              /*  check for descent direction */
              if (b_free[0]) {
                deltaX_idx_0 = d_tmp_data[0] - u[0];
                k = 1;
              }

              if (b_free[1]) {
                deltaX_idx_1 = d_tmp_data[k] - u[1];
                k++;
              }

              if (b_free[2]) {
                deltaX_idx_2 = d_tmp_data[k] - u[2];
              }

              sdotg = (deltaX_idx_0 * grad[0] + deltaX_idx_1 * grad[1]) +
                deltaX_idx_2 * grad[2];
              if (sdotg >= 0.0) {
                /*  (should not happen) */
                exitg1 = true;
              } else {
                /*  Armijo linesearch */
                step = 1.0;

                /*  Returns array x with all values clamped between lower and upper */
                ssq = fmax(lower[0], fmin(upper[0], u[0] + deltaX_idx_0));
                z1[0] = ssq;
                uc_idx_0 = ssq;
                c = ssq * Qu[0];
                ssq = fmax(lower[1], fmin(upper[1], u[1] + deltaX_idx_1));
                z1[1] = ssq;
                uc_idx_1 = ssq;
                c += ssq * Qu[1];
                ssq = fmax(lower[2], fmin(upper[2], u[2] + deltaX_idx_2));
                z1[2] = ssq;
                uc_idx_2 = ssq;
                c += ssq * Qu[2];
                t = 0.0;
                for (i = 0; i < 3; i++) {
                  t += ((0.5 * z1[0] * Quu[3 * i] + 0.5 * z1[1] * Quu[3 * i + 1])
                        + 0.5 * ssq * Quu[3 * i + 2]) * z1[i];
                }

                vc = c + t;
                exitg2 = false;
                while ((!exitg2) && ((vc - value) / (step * sdotg) < 0.1)) {
                  step *= 0.6;

                  /*  Returns array x with all values clamped between lower and upper */
                  ssq = fmax(lower[0], fmin(upper[0], u[0] + step * deltaX_idx_0));
                  z1[0] = ssq;
                  uc_idx_0 = ssq;
                  c = ssq * Qu[0];
                  ssq = fmax(lower[1], fmin(upper[1], u[1] + step * deltaX_idx_1));
                  z1[1] = ssq;
                  uc_idx_1 = ssq;
                  c += ssq * Qu[1];
                  ssq = fmax(lower[2], fmin(upper[2], u[2] + step * deltaX_idx_2));
                  z1[2] = ssq;
                  uc_idx_2 = ssq;
                  c += ssq * Qu[2];
                  t = 0.0;
                  for (i = 0; i < 3; i++) {
                    t += ((0.5 * z1[0] * Quu[3 * i] + 0.5 * z1[1] * Quu[3 * i +
                           1]) + 0.5 * ssq * Quu[3 * i + 2]) * z1[i];
                  }

                  vc = c + t;
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
 * Uses an rk method to roll out a trajectory
 *  Returns the new trajectory, cost and the derivatives along the trajectory
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
  int b_tmp;
  double b[6];
  int b_b_tmp;
  double c;
  double b_z1[3];
  int c_b_tmp;
  double q1[4];
  static const signed char b_iv[4] = { 1, -1, -1, -1 };

  double qdot_tmp[9];
  double b_q1[16];
  double y;
  double absxk;
  double scale;
  int b_k;
  double d;
  int q1_tmp;
  double quat_error[4];
  double t;
  double xnew_tmp;
  int qdot_tmp_tmp;
  int b_qdot_tmp_tmp;
  int c_qdot_tmp_tmp;
  double wdot_tmp[9];
  double b_wdot_tmp[9];
  double c_wdot_tmp[9];
  int i1;
  double dxdot1[70];
  double a[9];
  double b_dv[12];
  signed char i2;
  static const signed char b_a[9] = { -100, 0, 0, 0, -100, 0, 0, 0, -100 };

  signed char i3;
  int dxdot1_tmp;
  double c_a[3];
  static const signed char d_a[9] = { 100, 0, 0, 0, 100, 0, 0, 0, 100 };

  double b_x[7];
  static const signed char iv1[21] = { 0, 0, 0, 0, 100, 0, 0, 0, 0, 0, 0, 0, 100,
    0, 0, 0, 0, 0, 0, 0, 100 };

  double dxdot2[70];
  double x1[7];
  double E1[42];
  static const signed char iv2[18] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1,
    0, 0, 0, 1 };

  signed char b_I[49];
  double fx_tmp[49];
  double c_I[49];
  double b_xnew[42];
  double b_E1[42];
  double b_dv1[21];

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
    b_tmp = 7 * k + 4;
    b[3] = xnew[b_tmp] - x[b_tmp];
    b_b_tmp = 7 * k + 5;
    b[4] = xnew[b_b_tmp] - x[b_b_tmp];
    c_b_tmp = 7 * k + 6;
    b[5] = xnew[c_b_tmp] - x[c_b_tmp];

    /*  Calculate error between q1 and q2 */
    /*  Defined as conj(q1)*q2 */
    for (i = 0; i < 4; i++) {
      q1[i] = (double)b_iv[i] * xnew[i + 7 * k];
    }

    /*  conjugate */
    /*  Returns skew symmetric - cross product matrix of a 3x1 vector */
    qdot_tmp[0] = 0.0;
    qdot_tmp[3] = -q1[3];
    qdot_tmp[6] = q1[2];
    qdot_tmp[1] = q1[3];
    qdot_tmp[4] = 0.0;
    qdot_tmp[7] = -q1[1];
    qdot_tmp[2] = -q1[2];
    qdot_tmp[5] = q1[1];
    qdot_tmp[8] = 0.0;
    b_q1[0] = q1[0];
    for (i = 0; i < 3; i++) {
      absxk = q1[i + 1];
      b_k = (i + 1) << 2;
      b_q1[b_k] = -absxk;
      b_q1[i + 1] = absxk;
      b_q1[b_k + 1] = q1[0] * (double)iv[3 * i] + qdot_tmp[3 * i];
      q1_tmp = 3 * i + 1;
      b_q1[b_k + 2] = q1[0] * (double)iv[q1_tmp] + qdot_tmp[q1_tmp];
      q1_tmp = 3 * i + 2;
      b_q1[b_k + 3] = q1[0] * (double)iv[q1_tmp] + qdot_tmp[q1_tmp];
    }

    y = 0.0;
    scale = 3.3121686421112381E-170;
    for (b_k = 0; b_k < 4; b_k++) {
      d = 0.0;
      for (i = 0; i < 4; i++) {
        d += b_q1[b_k + (i << 2)] * x[i + 7 * k];
      }

      quat_error[b_k] = d;
      absxk = fabs(d);
      if (absxk > scale) {
        t = scale / absxk;
        y = y * t * t + 1.0;
        scale = absxk;
      } else {
        t = absxk / scale;
        y += t * t;
      }
    }

    y = scale * sqrt(y);
    for (i = 0; i < 4; i++) {
      quat_error[i] /= y;
    }

    /*  re-normalzie */
    b[0] = quat_error[1] / quat_error[0];
    b[1] = quat_error[2] / quat_error[0];
    b[2] = quat_error[3] / quat_error[0];

    /*  inverse Cayley Map */
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
        d += K[(b_k + 3 * i) + 18 * k] * b[i];
      }

      i = b_k + 3 * k;
      d = fmin(u_lims[b_k + 3], fmax(u_lims[b_k], (u[i] - alpha * l[i]) - d));
      unew[i] = d;
      b_z1[b_k] = -xnew[(b_k + 7 * k) + 1];
    }

    /*  Returns skew symmetric - cross product matrix of a 3x1 vector */
    xnew_tmp = xnew[7 * k];
    qdot_tmp[0] = 0.0;
    qdot_tmp_tmp = 7 * k + 3;
    qdot_tmp[3] = -xnew[qdot_tmp_tmp];
    b_qdot_tmp_tmp = 7 * k + 2;
    qdot_tmp[6] = xnew[b_qdot_tmp_tmp];
    qdot_tmp[1] = xnew[qdot_tmp_tmp];
    qdot_tmp[4] = 0.0;
    c_qdot_tmp_tmp = 7 * k + 1;
    qdot_tmp[7] = -xnew[c_qdot_tmp_tmp];
    qdot_tmp[2] = -xnew[b_qdot_tmp_tmp];
    qdot_tmp[5] = xnew[c_qdot_tmp_tmp];
    qdot_tmp[8] = 0.0;
    for (i = 0; i < 9; i++) {
      qdot_tmp[i] += xnew_tmp * (double)iv[i];
    }

    /*  Returns skew symmetric - cross product matrix of a 3x1 vector */
    wdot_tmp[0] = 0.0;
    wdot_tmp[3] = -xnew[c_b_tmp];
    wdot_tmp[6] = xnew[b_b_tmp];
    wdot_tmp[1] = xnew[c_b_tmp];
    wdot_tmp[4] = 0.0;
    wdot_tmp[7] = -xnew[b_tmp];
    wdot_tmp[2] = -xnew[b_b_tmp];
    wdot_tmp[5] = xnew[b_tmp];
    wdot_tmp[8] = 0.0;

    /*  Jacobians  */
    for (i = 0; i < 3; i++) {
      d = wdot_tmp[i + 3];
      absxk = wdot_tmp[i + 6];
      t = 0.0;
      for (i1 = 0; i1 < 3; i1++) {
        q1_tmp = i + 3 * i1;
        c_wdot_tmp[q1_tmp] = (wdot_tmp[i] * dv[3 * i1] + d * dv[3 * i1 + 1]) +
          absxk * dv[3 * i1 + 2];
        t += dv[q1_tmp] * xnew[(i1 + 7 * k) + 4];
      }

      z1[i] = t;
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
      q1_tmp = 7 * (i + 1);
      dxdot1[q1_tmp] = -d;
      b_k = 7 * (i + 4);
      dxdot1[b_k] = b_z1[i];
      dxdot1[i + 1] = d;
      i2 = b_a[i + 3];
      i3 = b_a[i + 6];
      for (i1 = 0; i1 < 3; i1++) {
        a[i + 3 * i1] = ((double)b_a[i] * b_wdot_tmp[3 * i1] + (double)i2 *
                         b_wdot_tmp[3 * i1 + 1]) + (double)i3 * b_wdot_tmp[3 *
          i1 + 2];
        dxdot1_tmp = i1 + 3 * i;
        dxdot1[(i1 + q1_tmp) + 1] = -wdot_tmp[dxdot1_tmp];
        dxdot1[(i1 + b_k) + 1] = qdot_tmp[dxdot1_tmp];
      }
    }

    for (i = 0; i < 4; i++) {
      dxdot1[7 * i + 4] = 0.0;
      dxdot1[7 * i + 5] = 0.0;
      dxdot1[7 * i + 6] = 0.0;
    }

    for (i = 0; i < 3; i++) {
      q1_tmp = 7 * (i + 4);
      dxdot1[q1_tmp + 4] = a[3 * i];
      b_k = 3 * i + 1;
      dxdot1[q1_tmp + 5] = a[b_k];
      dxdot1_tmp = 3 * i + 2;
      dxdot1[q1_tmp + 6] = a[dxdot1_tmp];
      for (i1 = 0; i1 < 7; i1++) {
        dxdot1[i1 + 7 * (i + 7)] = iv1[i1 + 7 * i];
      }

      i1 = i << 2;
      b_dv[i1] = 0.5 * b_z1[i];
      b_dv[i1 + 1] = 0.5 * qdot_tmp[3 * i];
      b_dv[i1 + 2] = 0.5 * qdot_tmp[b_k];
      b_dv[i1 + 3] = 0.5 * qdot_tmp[dxdot1_tmp];
    }

    for (i = 0; i < 4; i++) {
      q1[i] = (b_dv[i] * xnew[b_tmp] + b_dv[i + 4] * xnew[b_b_tmp]) + b_dv[i + 8]
        * xnew[c_b_tmp];
    }

    for (i = 0; i < 3; i++) {
      z1[i] = unew[i + 3 * k] - ((c_wdot_tmp[i] * xnew[b_tmp] + c_wdot_tmp[i + 3]
        * xnew[b_b_tmp]) + c_wdot_tmp[i + 6] * xnew[c_b_tmp]);
    }

    for (i = 0; i < 3; i++) {
      c_a[i] = ((double)d_a[i] * z1[0] + (double)d_a[i + 3] * z1[1]) + (double)
        d_a[i + 6] * z1[2];
    }

    for (i = 0; i < 4; i++) {
      b_x[i] = xnew[i + 7 * k] + 0.015 * q1[i];
    }

    /*  Calculates the continuous time state derivative and Jacobians */
    /*  kgm^2 */
    /*  Angular velocity */
    /*  Quaternion components */
    /*  Non-linear dynamics */
    b_x[4] = xnew[b_tmp] + 0.015 * c_a[0];
    b_z1[0] = -b_x[1];
    b_x[5] = xnew[b_b_tmp] + 0.015 * c_a[1];
    b_z1[1] = -b_x[2];
    b_x[6] = xnew[c_b_tmp] + 0.015 * c_a[2];
    b_z1[2] = -b_x[3];

    /*  Returns skew symmetric - cross product matrix of a 3x1 vector */
    qdot_tmp[0] = 0.0;
    qdot_tmp[3] = -b_x[3];
    qdot_tmp[6] = b_x[2];
    qdot_tmp[1] = b_x[3];
    qdot_tmp[4] = 0.0;
    qdot_tmp[7] = -b_x[1];
    qdot_tmp[2] = -b_x[2];
    qdot_tmp[5] = b_x[1];
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
      absxk = wdot_tmp[i + 6];
      t = 0.0;
      for (i1 = 0; i1 < 3; i1++) {
        q1_tmp = i + 3 * i1;
        c_wdot_tmp[q1_tmp] = (wdot_tmp[i] * dv[3 * i1] + d * dv[3 * i1 + 1]) +
          absxk * dv[3 * i1 + 2];
        t += dv[q1_tmp] * b_x[i1 + 4];
      }

      z1[i] = t;
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
      absxk = b_x[i + 4];
      q1_tmp = 7 * (i + 1);
      dxdot2[q1_tmp] = -absxk;
      b_k = 7 * (i + 4);
      dxdot2[b_k] = b_z1[i];
      dxdot2[i + 1] = absxk;
      i2 = b_a[i + 3];
      i3 = b_a[i + 6];
      for (i1 = 0; i1 < 3; i1++) {
        a[i + 3 * i1] = ((double)b_a[i] * b_wdot_tmp[3 * i1] + (double)i2 *
                         b_wdot_tmp[3 * i1 + 1]) + (double)i3 * b_wdot_tmp[3 *
          i1 + 2];
        dxdot1_tmp = i1 + 3 * i;
        dxdot2[(i1 + q1_tmp) + 1] = -wdot_tmp[dxdot1_tmp];
        dxdot2[(i1 + b_k) + 1] = qdot_tmp[dxdot1_tmp];
      }
    }

    for (i = 0; i < 4; i++) {
      dxdot2[7 * i + 4] = 0.0;
      dxdot2[7 * i + 5] = 0.0;
      dxdot2[7 * i + 6] = 0.0;
    }

    for (i = 0; i < 3; i++) {
      q1_tmp = 7 * (i + 4);
      dxdot2[q1_tmp + 4] = a[3 * i];
      b_k = 3 * i + 1;
      dxdot2[q1_tmp + 5] = a[b_k];
      dxdot1_tmp = 3 * i + 2;
      dxdot2[q1_tmp + 6] = a[dxdot1_tmp];
      for (i1 = 0; i1 < 7; i1++) {
        dxdot2[i1 + 7 * (i + 7)] = iv1[i1 + 7 * i];
      }

      i1 = i << 2;
      b_dv[i1] = 0.5 * b_z1[i];
      b_dv[i1 + 1] = 0.5 * qdot_tmp[3 * i];
      b_dv[i1 + 2] = 0.5 * qdot_tmp[b_k];
      b_dv[i1 + 3] = 0.5 * qdot_tmp[dxdot1_tmp];
    }

    d = b_x[4];
    absxk = b_x[5];
    t = b_x[6];
    for (i = 0; i < 4; i++) {
      q1[i] = (b_dv[i] * d + b_dv[i + 4] * absxk) + b_dv[i + 8] * t;
    }

    for (i = 0; i < 3; i++) {
      z1[i] = unew[i + 3 * k] - ((c_wdot_tmp[i] * d + c_wdot_tmp[i + 3] * absxk)
        + c_wdot_tmp[i + 6] * t);
    }

    for (i = 0; i < 3; i++) {
      c_a[i] = ((double)d_a[i] * z1[0] + (double)d_a[i + 3] * z1[1]) + (double)
        d_a[i + 6] * z1[2];
    }

    x1[4] = xnew[b_tmp] + 0.03 * c_a[0];
    x1[5] = xnew[b_b_tmp] + 0.03 * c_a[1];
    x1[6] = xnew[c_b_tmp] + 0.03 * c_a[2];

    /*  Normalize the quaternion */
    y = 0.0;
    scale = 3.3121686421112381E-170;
    for (b_k = 0; b_k < 4; b_k++) {
      d = xnew[b_k + 7 * k] + 0.03 * q1[b_k];
      x1[b_k] = d;
      absxk = fabs(d);
      if (absxk > scale) {
        t = scale / absxk;
        y = y * t * t + 1.0;
        scale = absxk;
      } else {
        t = absxk / scale;
        y += t * t;
      }
    }

    y = scale * sqrt(y);

    /*  Continuous time Jacobians */
    /*  Discrete time Jacobians */
    /*  First form the state error Jacobians E(x_k) & E(x_k+1) */
    /*  Returns skew symmetric - cross product matrix of a 3x1 vector */
    /*  Returns skew symmetric - cross product matrix of a 3x1 vector */
    qdot_tmp[0] = 0.0;
    qdot_tmp[3] = -x1[3];
    qdot_tmp[6] = x1[2];
    qdot_tmp[1] = x1[3];
    qdot_tmp[4] = 0.0;
    qdot_tmp[7] = -x1[1];
    qdot_tmp[2] = -x1[2];
    qdot_tmp[5] = x1[1];
    qdot_tmp[8] = 0.0;
    for (i = 0; i < 3; i++) {
      E1[7 * i] = -x1[i + 1];
      b_k = 7 * (i + 3);
      E1[b_k] = 0.0;
      E1[7 * i + 1] = x1[0] * (double)iv[3 * i] + qdot_tmp[3 * i];
      E1[b_k + 1] = 0.0;
      q1_tmp = 3 * i + 1;
      E1[7 * i + 2] = x1[0] * (double)iv[q1_tmp] + qdot_tmp[q1_tmp];
      E1[b_k + 2] = 0.0;
      q1_tmp = 3 * i + 2;
      E1[7 * i + 3] = x1[0] * (double)iv[q1_tmp] + qdot_tmp[q1_tmp];
      E1[b_k + 3] = 0.0;
    }

    for (i = 0; i < 6; i++) {
      E1[7 * i + 4] = iv2[3 * i];
      E1[7 * i + 5] = iv2[3 * i + 1];
      E1[7 * i + 6] = iv2[3 * i + 2];
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
        d = 0.0;
        for (b_k = 0; b_k < 7; b_k++) {
          d += fx_tmp[i + 7 * b_k] * dxdot1[b_k + 7 * i1];
        }

        b_k = i + 7 * i1;
        c_I[b_k] = ((double)b_I[b_k] + 0.03 * dxdot2[b_k]) + d;
      }
    }

    qdot_tmp[0] = 0.0;
    qdot_tmp[3] = -xnew[qdot_tmp_tmp];
    qdot_tmp[6] = xnew[b_qdot_tmp_tmp];
    qdot_tmp[1] = xnew[qdot_tmp_tmp];
    qdot_tmp[4] = 0.0;
    qdot_tmp[7] = -xnew[c_qdot_tmp_tmp];
    qdot_tmp[2] = -xnew[b_qdot_tmp_tmp];
    qdot_tmp[5] = xnew[c_qdot_tmp_tmp];
    qdot_tmp[8] = 0.0;
    for (i = 0; i < 6; i++) {
      for (i1 = 0; i1 < 7; i1++) {
        d = 0.0;
        for (b_k = 0; b_k < 7; b_k++) {
          d += E1[b_k + 7 * i] * c_I[b_k + 7 * i1];
        }

        b_E1[i + 6 * i1] = d;
      }
    }

    for (i = 0; i < 3; i++) {
      b_xnew[7 * i] = -xnew[(i + 7 * k) + 1];
      b_k = 7 * (i + 3);
      b_xnew[b_k] = 0.0;
      b_xnew[7 * i + 1] = xnew_tmp * (double)iv[3 * i] + qdot_tmp[3 * i];
      b_xnew[b_k + 1] = 0.0;
      q1_tmp = 3 * i + 1;
      b_xnew[7 * i + 2] = xnew_tmp * (double)iv[q1_tmp] + qdot_tmp[q1_tmp];
      b_xnew[b_k + 2] = 0.0;
      q1_tmp = 3 * i + 2;
      b_xnew[7 * i + 3] = xnew_tmp * (double)iv[q1_tmp] + qdot_tmp[q1_tmp];
      b_xnew[b_k + 3] = 0.0;
    }

    for (i = 0; i < 6; i++) {
      b_xnew[7 * i + 4] = iv2[3 * i];
      b_xnew[7 * i + 5] = iv2[3 * i + 1];
      b_xnew[7 * i + 6] = iv2[3 * i + 2];
    }

    for (i = 0; i < 6; i++) {
      for (i1 = 0; i1 < 6; i1++) {
        d = 0.0;
        for (b_k = 0; b_k < 7; b_k++) {
          d += b_E1[i + 6 * b_k] * b_xnew[b_k + 7 * i1];
        }

        fx[(i + 6 * i1) + 36 * k] = d;
      }
    }

    for (i = 0; i < 7; i++) {
      for (i1 = 0; i1 < 3; i1++) {
        d = 0.0;
        for (b_k = 0; b_k < 7; b_k++) {
          d += fx_tmp[i + 7 * b_k] * dxdot1[b_k + 7 * (i1 + 7)];
        }

        b_dv1[i + 7 * i1] = 0.03 * dxdot2[i + 7 * (i1 + 7)] + d;
      }
    }

    for (i = 0; i < 6; i++) {
      for (i1 = 0; i1 < 3; i1++) {
        d = 0.0;
        for (b_k = 0; b_k < 7; b_k++) {
          d += E1[b_k + 7 * i] * b_dv1[b_k + 7 * i1];
        }

        fu[(i + 6 * i1) + 18 * k] = d;
      }
    }

    for (i = 0; i < 4; i++) {
      xnew[i + 7 * (k + 1)] = x1[i] / y;
    }

    /*  Calculate the cost */
    for (i = 0; i < 3; i++) {
      xnew[(i + 7 * (k + 1)) + 4] = x1[i + 4];
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
                 cx[2994], b_z1, *(double (*)[36])&cxx[17964]);
  *cost += c;
}

/*
 * Arguments    : const double A_data[]
 *                const int A_size[2]
 *                double B_data[]
 *                int B_size[1]
 * Return Type  : void
 */
static void mldivide(const double A_data[], const int A_size[2], double B_data[],
                     int B_size[1])
{
  int b_A_size[2];
  int mn;
  int info;
  double b_A_data[9];
  double tau_data[3];
  int tau_size[1];
  int jpvt_data[3];
  int jpvt_size[2];
  int rankA;
  int i;
  double b_B_data[3];
  double wj;
  int m;
  int b_i;
  int j;
  if ((A_size[0] == 0) || (A_size[1] == 0) || (B_size[0] == 0)) {
    B_size[0] = (signed char)A_size[1];
    info = (signed char)A_size[1];
    if (0 <= info - 1) {
      memset(&B_data[0], 0, info * sizeof(double));
    }
  } else if (A_size[0] == A_size[1]) {
    mn = A_size[1];
    b_A_size[0] = A_size[0];
    b_A_size[1] = A_size[1];
    info = A_size[0] * A_size[1];
    if (0 <= info - 1) {
      memcpy(&b_A_data[0], &A_data[0], info * sizeof(double));
    }

    xgetrf(A_size[1], A_size[1], b_A_data, b_A_size, A_size[1], jpvt_data,
           jpvt_size, &info);
    i = A_size[1];
    for (info = 0; info <= i - 2; info++) {
      if (jpvt_data[info] != info + 1) {
        wj = B_data[info];
        B_data[info] = B_data[jpvt_data[info] - 1];
        B_data[jpvt_data[info] - 1] = wj;
      }
    }

    for (info = 0; info < mn; info++) {
      m = mn * info;
      if (B_data[info] != 0.0) {
        i = info + 2;
        for (b_i = i; b_i <= mn; b_i++) {
          B_data[b_i - 1] -= B_data[info] * b_A_data[(b_i + m) - 1];
        }
      }
    }

    for (info = mn; info >= 1; info--) {
      m = mn * (info - 1);
      wj = B_data[info - 1];
      if (wj != 0.0) {
        B_data[info - 1] = wj / b_A_data[(info + m) - 1];
        for (b_i = 0; b_i <= info - 2; b_i++) {
          B_data[b_i] -= B_data[info - 1] * b_A_data[b_i + m];
        }
      }
    }
  } else {
    b_A_size[0] = A_size[0];
    b_A_size[1] = A_size[1];
    info = A_size[0] * A_size[1];
    if (0 <= info - 1) {
      memcpy(&b_A_data[0], &A_data[0], info * sizeof(double));
    }

    xgeqp3(b_A_data, b_A_size, tau_data, tau_size, jpvt_data, jpvt_size);
    rankA = rankFromQR(b_A_data, b_A_size);
    info = B_size[0];
    if (0 <= info - 1) {
      memcpy(&b_B_data[0], &B_data[0], info * sizeof(double));
    }

    B_size[0] = (signed char)b_A_size[1];
    info = (signed char)b_A_size[1];
    if (0 <= info - 1) {
      memset(&B_data[0], 0, info * sizeof(double));
    }

    m = b_A_size[0];
    info = b_A_size[0];
    mn = b_A_size[1];
    if (info < mn) {
      mn = info;
    }

    for (j = 0; j < mn; j++) {
      if (tau_data[j] != 0.0) {
        wj = b_B_data[j];
        i = j + 2;
        for (b_i = i; b_i <= m; b_i++) {
          wj += b_A_data[(b_i + b_A_size[0] * j) - 1] * b_B_data[b_i - 1];
        }

        wj *= tau_data[j];
        if (wj != 0.0) {
          b_B_data[j] -= wj;
          i = j + 2;
          for (b_i = i; b_i <= m; b_i++) {
            b_B_data[b_i - 1] -= b_A_data[(b_i + b_A_size[0] * j) - 1] * wj;
          }
        }
      }
    }

    for (b_i = 0; b_i < rankA; b_i++) {
      B_data[jpvt_data[b_i] - 1] = b_B_data[b_i];
    }

    for (j = rankA; j >= 1; j--) {
      i = jpvt_data[j - 1];
      info = b_A_size[0] * (j - 1);
      B_data[i - 1] /= b_A_data[(j + info) - 1];
      for (b_i = 0; b_i <= j - 2; b_i++) {
        B_data[jpvt_data[b_i] - 1] -= B_data[jpvt_data[j - 1] - 1] *
          b_A_data[b_i + info];
      }
    }
  }
}

/*
 * Arguments    : double A_data[]
 *                const int A_size[2]
 *                int m
 *                int n
 *                double tau_data[]
 *                int jpvt_data[]
 * Return Type  : void
 */
static void qrpf(double A_data[], const int A_size[2], int m, int n, double
                 tau_data[], int jpvt_data[])
{
  int ma;
  int minmn;
  int itemp;
  double work_data[3];
  double vn1_data[3];
  double vn2_data[3];
  int j;
  int i;
  double smax;
  int ip1;
  int iy;
  int ii;
  int nmi;
  int mmi;
  int ix;
  int pvt;
  double temp2;
  double s;
  int jA;
  int lastv;
  int b_i;
  boolean_T exitg2;
  int exitg1;
  ma = A_size[0];
  if (m < n) {
    minmn = m;
  } else {
    minmn = n;
  }

  itemp = A_size[1];
  if (0 <= itemp - 1) {
    memset(&work_data[0], 0, itemp * sizeof(double));
  }

  itemp = A_size[1];
  if (0 <= itemp - 1) {
    memset(&vn1_data[0], 0, itemp * sizeof(double));
  }

  itemp = A_size[1];
  if (0 <= itemp - 1) {
    memset(&vn2_data[0], 0, itemp * sizeof(double));
  }

  for (j = 0; j < n; j++) {
    smax = xnrm2(m, A_data, j * ma + 1);
    vn1_data[j] = smax;
    vn2_data[j] = smax;
  }

  for (i = 0; i < minmn; i++) {
    ip1 = i + 2;
    iy = i * ma;
    ii = iy + i;
    nmi = n - i;
    mmi = m - i;
    if (nmi < 1) {
      itemp = -1;
    } else {
      itemp = 0;
      if (nmi > 1) {
        ix = i;
        smax = fabs(vn1_data[i]);
        for (j = 2; j <= nmi; j++) {
          ix++;
          s = fabs(vn1_data[ix]);
          if (s > smax) {
            itemp = j - 1;
            smax = s;
          }
        }
      }
    }

    pvt = i + itemp;
    if (pvt + 1 != i + 1) {
      ix = pvt * ma;
      for (j = 0; j < m; j++) {
        smax = A_data[ix];
        A_data[ix] = A_data[iy];
        A_data[iy] = smax;
        ix++;
        iy++;
      }

      itemp = jpvt_data[pvt];
      jpvt_data[pvt] = jpvt_data[i];
      jpvt_data[i] = itemp;
      vn1_data[pvt] = vn1_data[i];
      vn2_data[pvt] = vn2_data[i];
    }

    if (i + 1 < m) {
      temp2 = A_data[ii];
      itemp = ii + 2;
      tau_data[i] = 0.0;
      if (mmi > 0) {
        smax = xnrm2(mmi - 1, A_data, ii + 2);
        if (smax != 0.0) {
          s = rt_hypotd_snf(A_data[ii], smax);
          if (A_data[ii] >= 0.0) {
            s = -s;
          }

          if (fabs(s) < 1.0020841800044864E-292) {
            pvt = -1;
            b_i = ii + mmi;
            do {
              pvt++;
              for (j = itemp; j <= b_i; j++) {
                A_data[j - 1] *= 9.9792015476736E+291;
              }

              s *= 9.9792015476736E+291;
              temp2 *= 9.9792015476736E+291;
            } while (!(fabs(s) >= 1.0020841800044864E-292));

            s = rt_hypotd_snf(temp2, xnrm2(mmi - 1, A_data, ii + 2));
            if (temp2 >= 0.0) {
              s = -s;
            }

            tau_data[i] = (s - temp2) / s;
            smax = 1.0 / (temp2 - s);
            for (j = itemp; j <= b_i; j++) {
              A_data[j - 1] *= smax;
            }

            for (j = 0; j <= pvt; j++) {
              s *= 1.0020841800044864E-292;
            }

            temp2 = s;
          } else {
            tau_data[i] = (s - A_data[ii]) / s;
            smax = 1.0 / (A_data[ii] - s);
            b_i = ii + mmi;
            for (j = itemp; j <= b_i; j++) {
              A_data[j - 1] *= smax;
            }

            temp2 = s;
          }
        }
      }

      A_data[ii] = temp2;
    } else {
      tau_data[i] = 0.0;
    }

    if (i + 1 < n) {
      temp2 = A_data[ii];
      A_data[ii] = 1.0;
      jA = (ii + ma) + 1;
      if (tau_data[i] != 0.0) {
        lastv = mmi - 1;
        itemp = (ii + mmi) - 1;
        while ((lastv + 1 > 0) && (A_data[itemp] == 0.0)) {
          lastv--;
          itemp--;
        }

        nmi -= 2;
        exitg2 = false;
        while ((!exitg2) && (nmi + 1 > 0)) {
          itemp = jA + nmi * ma;
          j = itemp;
          do {
            exitg1 = 0;
            if (j <= itemp + lastv) {
              if (A_data[j - 1] != 0.0) {
                exitg1 = 1;
              } else {
                j++;
              }
            } else {
              nmi--;
              exitg1 = 2;
            }
          } while (exitg1 == 0);

          if (exitg1 == 1) {
            exitg2 = true;
          }
        }
      } else {
        lastv = -1;
        nmi = -1;
      }

      if (lastv + 1 > 0) {
        if (nmi + 1 != 0) {
          if (0 <= nmi) {
            memset(&work_data[0], 0, (nmi + 1) * sizeof(double));
          }

          iy = 0;
          b_i = jA + ma * nmi;
          for (itemp = jA; ma < 0 ? itemp >= b_i : itemp <= b_i; itemp += ma) {
            ix = ii;
            smax = 0.0;
            pvt = itemp + lastv;
            for (j = itemp; j <= pvt; j++) {
              smax += A_data[j - 1] * A_data[ix];
              ix++;
            }

            work_data[iy] += smax;
            iy++;
          }
        }

        if (!(-tau_data[i] == 0.0)) {
          itemp = 0;
          for (j = 0; j <= nmi; j++) {
            if (work_data[itemp] != 0.0) {
              smax = work_data[itemp] * -tau_data[i];
              ix = ii;
              b_i = lastv + jA;
              for (pvt = jA; pvt <= b_i; pvt++) {
                A_data[pvt - 1] += A_data[ix] * smax;
                ix++;
              }
            }

            itemp++;
            jA += ma;
          }
        }
      }

      A_data[ii] = temp2;
    }

    for (j = ip1; j <= n; j++) {
      itemp = i + (j - 1) * ma;
      smax = vn1_data[j - 1];
      if (smax != 0.0) {
        s = fabs(A_data[itemp]) / smax;
        s = 1.0 - s * s;
        if (s < 0.0) {
          s = 0.0;
        }

        temp2 = smax / vn2_data[j - 1];
        temp2 = s * (temp2 * temp2);
        if (temp2 <= 1.4901161193847656E-8) {
          if (i + 1 < m) {
            smax = xnrm2(mmi - 1, A_data, itemp + 2);
            vn1_data[j - 1] = smax;
            vn2_data[j - 1] = smax;
          } else {
            vn1_data[j - 1] = 0.0;
            vn2_data[j - 1] = 0.0;
          }
        } else {
          vn1_data[j - 1] = smax * sqrt(s);
        }
      }
    }
  }
}

/*
 * Arguments    : const double A_data[]
 *                const int A_size[2]
 * Return Type  : int
 */
static int rankFromQR(const double A_data[], const int A_size[2])
{
  int r;
  int minmn;
  int maxmn;
  double tol;
  r = 0;
  if (A_size[0] < A_size[1]) {
    minmn = A_size[0];
    maxmn = A_size[1];
  } else {
    minmn = A_size[1];
    maxmn = A_size[0];
  }

  if (minmn > 0) {
    tol = 2.2204460492503131E-15 * (double)maxmn * fabs(A_data[0]);
    while ((r < minmn) && (!(fabs(A_data[r + A_size[0] * r]) <= tol))) {
      r++;
    }
  }

  return r;
}

/*
 * Arguments    : double u0
 *                double u1
 * Return Type  : double
 */
static double rt_hypotd_snf(double u0, double u1)
{
  double y;
  double a;
  a = fabs(u0);
  y = fabs(u1);
  if (a < y) {
    a /= y;
    y *= sqrt(a * a + 1.0);
  } else if (a > y) {
    y /= a;
    y = a * sqrt(y * y + 1.0);
  } else {
    if (!rtIsNaN(y)) {
      y = a * 1.4142135623730951;
    }
  }

  return y;
}

/*
 * Arguments    : double u0
 *                double u1
 * Return Type  : double
 */
static double rt_powd_snf(double u0, double u1)
{
  double y;
  double d;
  double d1;
  if (rtIsNaN(u0) || rtIsNaN(u1)) {
    y = rtNaN;
  } else {
    d = fabs(u0);
    d1 = fabs(u1);
    if (rtIsInf(u1)) {
      if (d == 1.0) {
        y = 1.0;
      } else if (d > 1.0) {
        if (u1 > 0.0) {
          y = rtInf;
        } else {
          y = 0.0;
        }
      } else if (u1 > 0.0) {
        y = 0.0;
      } else {
        y = rtInf;
      }
    } else if (d1 == 0.0) {
      y = 1.0;
    } else if (d1 == 1.0) {
      if (u1 > 0.0) {
        y = u0;
      } else {
        y = 1.0 / u0;
      }
    } else if (u1 == 2.0) {
      y = u0 * u0;
    } else if ((u1 == 0.5) && (u0 >= 0.0)) {
      y = sqrt(u0);
    } else if ((u0 < 0.0) && (u1 > floor(u1))) {
      y = rtNaN;
    } else {
      y = pow(u0, u1);
    }
  }

  return y;
}

/*
 * Calculates the cost contribution of a given state and control
 *  Also calculates the 2nd order expansion of the cost function
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
  double b_xg;
  double quat_cost;
  int b_sign;
  double b[3];
  double b_dv[3];
  int cxx_tmp;
  int sign_tmp;
  double b_dv1[9];
  double c_sign[12];
  int b_cxx_tmp;

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
  b_xg = 0.0;
  for (i = 0; i < 4; i++) {
    b_xg += xg[i] * x[i];
  }

  if (b_xg + 1.0 < 1.0 - b_xg) {
    quat_cost = b_xg + 1.0;
    b_sign = 1;
  } else {
    quat_cost = 1.0 - b_xg;
    b_sign = -1;
  }

  b[0] = x[4] - xg[4];
  b[1] = x[5] - xg[5];
  b[2] = x[6] - xg[6];
  b_xg = 0.0;
  for (i = 0; i < 3; i++) {
    cxx_tmp = 3 * i + 1;
    sign_tmp = 3 * i + 2;
    b_xg += ((0.5 * b[0] * Qw[3 * i] + 0.5 * b[1] * Qw[cxx_tmp]) + 0.5 * b[2] *
             Qw[sign_tmp]) * b[i];
    b_dv[i] = (0.5 * u[0] * dv1[3 * i] + 0.5 * u[1] * dv1[cxx_tmp]) + 0.5 * u[2]
      * dv1[sign_tmp];
  }

  *cost = (quat_cost + b_xg) + ((b_dv[0] * u[0] + b_dv[1] * u[1]) + b_dv[2] * u
    [2]);

  /*  State cost Hessian */
  b_xg = 0.0;
  for (i = 0; i < 4; i++) {
    b_xg += xg[i] * x[i];
  }

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
    cxx[6 * i] = (double)(b_sign * iv[3 * i]) * b_xg;
    cxx_tmp = 6 * (i + 3);
    cxx[cxx_tmp] = 0.0;
    cxx[6 * i + 3] = 0.0;
    cxx[cxx_tmp + 3] = Qw[3 * i];
    sign_tmp = 3 * (i + 1);
    c_sign[sign_tmp] = (double)b_sign * (x[0] * (double)iv[i] + b_dv1[i]);
    b_cxx_tmp = 3 * i + 1;
    cxx[6 * i + 1] = (double)(b_sign * iv[b_cxx_tmp]) * b_xg;
    cxx[cxx_tmp + 1] = 0.0;
    cxx[6 * i + 4] = 0.0;
    cxx[cxx_tmp + 4] = Qw[b_cxx_tmp];
    c_sign[sign_tmp + 1] = (double)b_sign * (x[0] * (double)iv[i + 3] + b_dv1[i
      + 3]);
    b_cxx_tmp = 3 * i + 2;
    cxx[6 * i + 2] = (double)(b_sign * iv[b_cxx_tmp]) * b_xg;
    cxx[cxx_tmp + 2] = 0.0;
    cxx[6 * i + 5] = 0.0;
    cxx[cxx_tmp + 5] = Qw[b_cxx_tmp];
    c_sign[sign_tmp + 2] = (double)b_sign * (x[0] * (double)iv[i + 6] + b_dv1[i
      + 6]);
  }

  /*  Control cost Hessian & Jacobian */
  for (i = 0; i < 3; i++) {
    b_xg = 0.0;
    for (cxx_tmp = 0; cxx_tmp < 4; cxx_tmp++) {
      b_xg += c_sign[i + 3 * cxx_tmp] * xg[cxx_tmp];
    }

    cx[i] = b_xg;
    cx[i + 3] = (Qw[i] * b[0] + Qw[i + 3] * b[1]) + Qw[i + 6] * b[2];
    cu[i] = (dv1[i] * u[0] + dv1[i + 3] * u[1]) + dv1[i + 6] * u[2];
  }
}

/*
 * Arguments    : double A_data[]
 *                const int A_size[2]
 *                double tau_data[]
 *                int tau_size[1]
 *                int jpvt_data[]
 *                int jpvt_size[2]
 * Return Type  : void
 */
static void xgeqp3(double A_data[], const int A_size[2], double tau_data[], int
                   tau_size[1], int jpvt_data[], int jpvt_size[2])
{
  int n;
  int k;
  int u1;
  n = A_size[1] - 1;
  k = A_size[0];
  u1 = A_size[1];
  if (k < u1) {
    u1 = k;
  }

  tau_size[0] = u1;
  if (0 <= u1 - 1) {
    memset(&tau_data[0], 0, u1 * sizeof(double));
  }

  if ((A_size[0] == 0) || (A_size[1] == 0) || (u1 < 1)) {
    jpvt_size[0] = 1;
    jpvt_size[1] = A_size[1];
    k = A_size[1];
    if (0 <= k - 1) {
      memset(&jpvt_data[0], 0, k * sizeof(int));
    }

    for (k = 0; k <= n; k++) {
      jpvt_data[k] = k + 1;
    }
  } else {
    jpvt_size[0] = 1;
    jpvt_size[1] = A_size[1];
    k = A_size[1];
    if (0 <= k - 1) {
      memset(&jpvt_data[0], 0, k * sizeof(int));
    }

    for (k = 0; k <= n; k++) {
      jpvt_data[k] = k + 1;
    }

    qrpf(A_data, A_size, A_size[0], A_size[1], tau_data, jpvt_data);
  }
}

/*
 * Arguments    : int m
 *                int n
 *                double A_data[]
 *                const int A_size[2]
 *                int lda
 *                int ipiv_data[]
 *                int ipiv_size[2]
 *                int *info
 * Return Type  : void
 */
static void xgetrf(int m, int n, double A_data[], const int A_size[2], int lda,
                   int ipiv_data[], int ipiv_size[2], int *info)
{
  int yk;
  int b_n;
  int jA;
  int u0;
  int j;
  int mmj;
  int b_tmp;
  int jp1j;
  int ix;
  double smax;
  int i;
  double s;
  int i1;
  int ijA;
  if (m < n) {
    yk = m;
  } else {
    yk = n;
  }

  if (yk < 1) {
    b_n = 0;
  } else {
    b_n = yk;
  }

  ipiv_size[0] = 1;
  ipiv_size[1] = b_n;
  if (b_n > 0) {
    ipiv_data[0] = 1;
    yk = 1;
    for (jA = 2; jA <= b_n; jA++) {
      yk++;
      ipiv_data[jA - 1] = yk;
    }
  }

  *info = 0;
  if ((m >= 1) && (n >= 1)) {
    u0 = m - 1;
    if (u0 >= n) {
      u0 = n;
    }

    for (j = 0; j < u0; j++) {
      mmj = m - j;
      b_tmp = j * (lda + 1);
      jp1j = b_tmp + 2;
      if (mmj < 1) {
        yk = -1;
      } else {
        yk = 0;
        if (mmj > 1) {
          ix = b_tmp;
          smax = fabs(A_data[b_tmp]);
          for (jA = 2; jA <= mmj; jA++) {
            ix++;
            s = fabs(A_data[ix]);
            if (s > smax) {
              yk = jA - 1;
              smax = s;
            }
          }
        }
      }

      if (A_data[b_tmp + yk] != 0.0) {
        if (yk != 0) {
          yk += j;
          ipiv_data[j] = yk + 1;
          ix = j;
          for (jA = 0; jA < n; jA++) {
            smax = A_data[ix];
            A_data[ix] = A_data[yk];
            A_data[yk] = smax;
            ix += lda;
            yk += lda;
          }
        }

        i = b_tmp + mmj;
        for (yk = jp1j; yk <= i; yk++) {
          A_data[yk - 1] /= A_data[b_tmp];
        }
      } else {
        *info = j + 1;
      }

      b_n = n - j;
      yk = b_tmp + lda;
      jA = yk;
      for (jp1j = 0; jp1j <= b_n - 2; jp1j++) {
        smax = A_data[yk];
        if (A_data[yk] != 0.0) {
          ix = b_tmp + 1;
          i = jA + 2;
          i1 = mmj + jA;
          for (ijA = i; ijA <= i1; ijA++) {
            A_data[ijA - 1] += A_data[ix] * -smax;
            ix++;
          }
        }

        yk += lda;
        jA += lda;
      }
    }

    if ((*info == 0) && (m <= n) && (!(A_data[(m + A_size[0] * (m - 1)) - 1] !=
          0.0))) {
      *info = m;
    }
  }
}

/*
 * Arguments    : int n
 *                const double x_data[]
 *                int ix0
 * Return Type  : double
 */
static double xnrm2(int n, const double x_data[], int ix0)
{
  double y;
  double scale;
  int kend;
  int k;
  double absxk;
  double t;
  y = 0.0;
  if (n >= 1) {
    if (n == 1) {
      y = fabs(x_data[ix0 - 1]);
    } else {
      scale = 3.3121686421112381E-170;
      kend = (ix0 + n) - 1;
      for (k = ix0; k <= kend; k++) {
        absxk = fabs(x_data[k - 1]);
        if (absxk > scale) {
          t = scale / absxk;
          y = y * t * t + 1.0;
          scale = absxk;
        } else {
          t = absxk / scale;
          y += t * t;
        }
      }

      y = scale * sqrt(y);
    }
  }

  return y;
}

/*
 * Solves finite horizon optimal control problem using a
 *  multiplicative iterative linear quadratic regulator
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
  static double l[1497];
  double dV[2];
  static double y[1497];
  static double b_dv[8982];
  static double fx[17964];
  static double fu[8982];
  static double cx[3000];
  double cu[1497];
  static double cxx[18000];
  static double cuu[4491];
  double cost;
  double cost_n;
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
    Alphas[k] = rt_powd_snf(10.0, -0.3 * (double)k);
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

  /*  Initial Forward rollout */
  memset(&l[0], 0, 1497U * sizeof(double));
  memset(&K[0], 0, 8982U * sizeof(double));
  dV[0] = 0.0;
  dV[1] = 0.0;
  memset(&y[0], 0, 1497U * sizeof(double));
  memset(&b_dv[0], 0, 8982U * sizeof(double));
  forwardRollout(x0, xg, u0, y, b_dv, 0.0, u_lims, x, u, fx, fu, cx, cu, cxx,
                 cuu, &cost);
  cost_n = 0.0;

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
      if ((!rtIsNaN(d)) && (rtIsNaN(b_y) || (b_y < d))) {
        b_y = d;
      }

      d = y[3 * k + 2];
      if ((!rtIsNaN(d)) && (rtIsNaN(b_y) || (b_y < d))) {
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
            } else if (z > 0.0) {
              z = 1.0;
            } else {
              if (z == 0.0) {
                z = 0.0;
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

        /*  Terminate ? */
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
 * Arguments    : void
 * Return Type  : void
 */
void milqr_initialize(void)
{
  rt_InitInfAndNaN();
}

/*
 * File trailer for milqr.c
 *
 * [EOF]
 */
