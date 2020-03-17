/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: ilqrCar.c
 *
 * MATLAB Coder version            : 4.3
 * C/C++ source code generated on  : 01-Mar-2020 20:37:52
 */

/* Include Files */
#include "ilqrCar.h"
#include <math.h>
#include <string.h>

/* Variable Definitions */
static const double dv[4] = { 0.01, 0.0, 0.0, 0.0001 };

/* Function Declarations */
static void b_mldivide(const double A_data[], const int A_size[2], double
  B_data[], int B_size[2]);

static void backwardPass(const double fx[8000], const double fu[4000], const
  double cx[2004], const double cu[1000], const double cxx[8016], const double
  cuu[2000], double lambda, const double u_lims[4], const double u[1000], double
  l[1000], double K[4000], double dV[2], bool *diverge);

static void boxQPsolve(const double Quu[4], const double Qu[2], const double
  lower[2], const double upper[2], const double u0[2], double u[2], double
  *result, double Luu[4], bool b_free[2]);

static void car_cost(const double x[4], const double xg[4], const double u[2],
                     double terminal, double *cost, double cx[4], double cu[2],
                     double cxx[16]);

static void forwardRollout(const double x[2004], const double xg[4], const
  double u[1000], const double l[1000], const double K[4000], double alpha,
  const double u_lims[4], double xnew[2004], double unew[1000], double fx[8000],
  double fu[4000], double cx[2004], double cu[1000], double cxx[8016], double
  cuu[2000], double *cost);

static void mldivide(const double A_data[], const int A_size[2], double B_data[],
                     int B_size[1]);

static void qrpf(double A_data[], const int A_size[2], int m, int n, double
                 tau_data[], int jpvt_data[]);

static int rankFromQR(const double A_data[], const int A_size[2]);
static double rt_hypotd(double u0, double u1);

static void xgeqp3(double A_data[], const int A_size[2], double tau_data[], int
                   tau_size[1], int jpvt_data[], int jpvt_size[2]);

static void xgetrf(int m, int n, double A_data[], const int A_size[2], int lda,
                   int ipiv_data[], int ipiv_size[2], int *info);

static double xnrm2(int n, const double x_data[], int ix0);




static void b_mldivide(const double A_data[], const int A_size[2], double
  B_data[], int B_size[2])
{
  int b_A_size[2];
  int mn;
  int jBcol;
  int i;
  double b_A_data[4];
  int i1;
  double tau_data[2];
  int tau_size[1];
  int jpvt_data[2];
  int jpvt_size[2];
  int rankA;
  int B_size_idx_0;
  double b_B_data[8];
  int j;
  int m;
  int wj_tmp;
  int k;
  double wj;
  int b_i;
  if ((A_size[0] == 0) || (A_size[1] == 0) || (B_size[0] == 0)) {
    B_size[0] = (signed char)A_size[1];
    B_size[1] = 4;
    jBcol = (signed char)A_size[1];
    for (i = 0; i < 4; i++) {
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
    for (jBcol = 0; jBcol <= i - 2; jBcol++) {
      if (jpvt_data[0] != 1) {
        for (m = 0; m < 4; m++) {
          wj_tmp = B_size[0] * m;
          wj = B_data[wj_tmp];
          i1 = (jpvt_data[0] + wj_tmp) - 1;
          B_data[wj_tmp] = B_data[i1];
          B_data[i1] = wj;
        }
      }
    }

    if (B_size[0] != 0) {
      for (j = 0; j < 4; j++) {
        jBcol = mn * j;
        for (k = 0; k < mn; k++) {
          m = mn * k;
          i = k + jBcol;
          if (B_data[i] != 0.0) {
            i1 = k + 2;
            for (b_i = i1; b_i <= mn; b_i++) {
              B_data[jBcol + 1] -= B_data[i] * b_A_data[m + 1];
            }
          }
        }
      }
    }

    if (B_size[0] != 0) {
      for (j = 0; j < 4; j++) {
        jBcol = mn * j - 1;
        for (k = mn; k >= 1; k--) {
          m = mn * (k - 1);
          i = k + jBcol;
          if (B_data[i] != 0.0) {
            B_data[i] /= b_A_data[(k + m) - 1];
            for (b_i = 0; b_i <= k - 2; b_i++) {
              B_data[jBcol + 1] -= B_data[i] * b_A_data[m];
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
    B_size[1] = 4;
    jBcol = (signed char)b_A_size[1];
    for (i = 0; i < 4; i++) {
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
        for (k = 0; k < 4; k++) {
          wj_tmp = B_size_idx_0 * k;
          jBcol = j + wj_tmp;
          wj = b_B_data[jBcol];
          for (b_i = i; b_i <= m; b_i++) {
            wj += b_A_data[b_A_size[0] * j + 1] * b_B_data[wj_tmp + 1];
          }

          wj *= tau_data[j];
          if (wj != 0.0) {
            b_B_data[jBcol] -= wj;
            i1 = j + 2;
            for (b_i = i1; b_i <= m; b_i++) {
              jBcol = B_size_idx_0 * k + 1;
              b_B_data[jBcol] -= b_A_data[b_A_size[0] * j + 1] * wj;
            }
          }
        }
      }
    }

    for (k = 0; k < 4; k++) {
      for (b_i = 0; b_i < rankA; b_i++) {
        B_data[(jpvt_data[b_i] + B_size[0] * k) - 1] = b_B_data[b_i +
          B_size_idx_0 * k];
      }

      for (j = rankA; j >= 1; j--) {
        i = B_size[0] * k;
        i1 = (jpvt_data[j - 1] + i) - 1;
        jBcol = b_A_size[0] * (j - 1);
        B_data[i1] /= b_A_data[(j + jBcol) - 1];
        for (b_i = 0; b_i <= j - 2; b_i++) {
          i1 = (jpvt_data[0] + i) - 1;
          B_data[i1] -= B_data[(jpvt_data[j - 1] + B_size[0] * k) - 1] *
            b_A_data[jBcol];
        }
      }
    }
  }
}




static void backwardPass(const double fx[8000], const double fu[4000], const
  double cx[2004], const double cu[1000], const double cxx[8016], const double
  cuu[2000], double lambda, const double u_lims[4], const double u[1000], double
  l[1000], double K[4000], double dV[2], bool *diverge)
{
  int i;
  double Vx[4];
  int k;
  int i1;
  bool exitg1;
  int Vxx_tmp;
  double Vxx[16];
  double result;
  int i2;
  int partialTrueCount;
  double Qu[2];
  double d;
  int u_lims_tmp;
  double Quu[4];
  double b_Quu[4];
  static const signed char a[4] = { 1, 0, 0, 1 };

  double b_u_lims[2];
  double Kk[8];
  double c_u_lims[2];
  double Qux[8];
  double b_dv[2];
  int b_u_lims_tmp;
  double lk[2];
  double Luu[4];
  bool b_free[2];
  bool y;
  bool exitg2;
  signed char tmp_data[2];
  double a_tmp[8];
  signed char b_tmp_data[2];
  double b_a_tmp[8];
  double Vx_tmp[8];
  int tmp_size[2];
  double Luu_data[4];
  int Luu_size[2];
  double b_fx[16];
  double c_tmp_data[8];
  double b_cxx[16];

  /*  Perfoms the LQR backward pass to find the optimal controls */
  /*  Solves a quadratic program (QP) at each timestep for the optimal */
  /*  controls given the control limits */
  /*  Initialize matrices (for C) */
  memset(&l[0], 0, 1000U * sizeof(double));
  memset(&K[0], 0, 4000U * sizeof(double));

  /*  Change in cost */
  dV[0] = 0.0;
  dV[1] = 0.0;

  /*  Set cost-to-go Jacobain and Hessian equal to final costs */
  for (i = 0; i < 4; i++) {
    Vx[i] = cx[i + 2000];
    for (i1 = 0; i1 < 4; i1++) {
      Vxx_tmp = i1 + (i << 2);
      Vxx[Vxx_tmp] = cxx[Vxx_tmp + 8000];
    }
  }

  *diverge = false;
  k = 0;
  exitg1 = false;
  while ((!exitg1) && (k < 500)) {
    /*  Define cost gradients */
    i = (499 - k) << 3;
    for (i1 = 0; i1 < 2; i1++) {
      result = 0.0;
      i2 = i1 << 2;
      for (partialTrueCount = 0; partialTrueCount < 4; partialTrueCount++) {
        result += fu[(partialTrueCount + i2) + i] * Vx[partialTrueCount];
        d = 0.0;
        for (Vxx_tmp = 0; Vxx_tmp < 4; Vxx_tmp++) {
          d += fu[(Vxx_tmp + i2) + i] * Vxx[Vxx_tmp + (partialTrueCount << 2)];
        }

        Kk[i1 + (partialTrueCount << 1)] = d;
      }

      Qu[i1] = cu[i1 + ((499 - k) << 1)] + result;
      for (i2 = 0; i2 < 2; i2++) {
        result = 0.0;
        for (partialTrueCount = 0; partialTrueCount < 4; partialTrueCount++) {
          result += Kk[i1 + (partialTrueCount << 1)] * fu[(partialTrueCount +
            (i2 << 2)) + ((499 - k) << 3)];
        }

        Vxx_tmp = i1 + (i2 << 1);
        b_Quu[Vxx_tmp] = cuu[Vxx_tmp + ((499 - k) << 2)] + result;
      }
    }

    for (i = 0; i < 2; i++) {
      for (i1 = 0; i1 < 4; i1++) {
        result = 0.0;
        for (i2 = 0; i2 < 4; i2++) {
          result += fu[(i2 + (i << 2)) + ((499 - k) << 3)] * Vxx[i2 + (i1 << 2)];
        }

        Kk[i + (i1 << 1)] = result;
      }

      for (i1 = 0; i1 < 4; i1++) {
        result = 0.0;
        for (i2 = 0; i2 < 4; i2++) {
          result += Kk[i + (i2 << 1)] * fx[(i2 + (i1 << 2)) + ((499 - k) << 4)];
        }

        Qux[i + (i1 << 1)] = result;
      }
    }

    /*  Regularization (for Cholesky positive definiteness) */
    /*  Solve the Quadratic program with control limits */
    for (i = 0; i < 4; i++) {
      Quu[i] = b_Quu[i] + (double)a[i] * lambda;
    }

    u_lims_tmp = (499 - k) << 1;
    b_u_lims[0] = u_lims[0] - u[u_lims_tmp];
    c_u_lims[0] = u_lims[2] - u[u_lims_tmp];
    i = ((int)fmin(500.0, (-(double)k + 500.0) + 1.0) - 1) << 1;
    b_dv[0] = -l[i];
    b_u_lims_tmp = u_lims_tmp + 1;
    b_u_lims[1] = u_lims[1] - u[b_u_lims_tmp];
    c_u_lims[1] = u_lims[3] - u[b_u_lims_tmp];
    b_dv[1] = -l[i + 1];
    boxQPsolve(Quu, Qu, b_u_lims, c_u_lims, b_dv, lk, &result, Luu, b_free);
    if (result < 1.0) {
      *diverge = true;
      exitg1 = true;
    } else {
      /*  Solve for feedback gains in non-clamped rows of u */
      /*  (using cholesky factor of Quu) */
      memset(&Kk[0], 0, 8U * sizeof(double));
      y = false;
      Vxx_tmp = 0;
      exitg2 = false;
      while ((!exitg2) && (Vxx_tmp < 2)) {
        if (b_free[Vxx_tmp]) {
          y = true;
          exitg2 = true;
        } else {
          Vxx_tmp++;
        }
      }

      if (y) {
        Vxx_tmp = 0;
        if (b_free[0]) {
          Vxx_tmp = 1;
        }

        if (b_free[1]) {
          Vxx_tmp++;
        }

        partialTrueCount = 0;
        if (b_free[0]) {
          tmp_data[0] = 1;
          partialTrueCount = 1;
        }

        if (b_free[1]) {
          tmp_data[partialTrueCount] = 2;
        }

        partialTrueCount = 0;
        if (b_free[0]) {
          b_tmp_data[0] = 1;
          partialTrueCount = 1;
        }

        if (b_free[1]) {
          b_tmp_data[partialTrueCount] = 2;
        }

        tmp_size[0] = Vxx_tmp;
        tmp_size[1] = 4;
        for (i = 0; i < 4; i++) {
          for (i1 = 0; i1 < Vxx_tmp; i1++) {
            c_tmp_data[i1 + Vxx_tmp * i] = Qux[(tmp_data[i1] + (i << 1)) - 1];
          }
        }

        Luu_size[0] = Vxx_tmp;
        Luu_size[1] = Vxx_tmp;
        for (i = 0; i < Vxx_tmp; i++) {
          for (i1 = 0; i1 < Vxx_tmp; i1++) {
            Luu_data[i1 + Vxx_tmp * i] = Luu[(tmp_data[i] + ((tmp_data[i1] - 1) <<
              1)) - 1];
          }
        }

        b_mldivide(Luu_data, Luu_size, c_tmp_data, tmp_size);
        Luu_size[0] = Vxx_tmp;
        Luu_size[1] = Vxx_tmp;
        for (i = 0; i < Vxx_tmp; i++) {
          for (i1 = 0; i1 < Vxx_tmp; i1++) {
            Luu_data[i1 + Vxx_tmp * i] = -Luu[(tmp_data[i1] + ((tmp_data[i] - 1)
              << 1)) - 1];
          }
        }

        b_mldivide(Luu_data, Luu_size, c_tmp_data, tmp_size);
        Vxx_tmp = tmp_size[0];
        for (i = 0; i < 4; i++) {
          for (i1 = 0; i1 < Vxx_tmp; i1++) {
            Kk[(b_tmp_data[i1] + (i << 1)) - 1] = c_tmp_data[i1 + tmp_size[0] *
              i];
          }
        }
      }

      /*  Update Cost to Go */
      result = 0.0;
      for (i = 0; i < 2; i++) {
        i1 = i << 1;
        result += (0.5 * lk[0] * b_Quu[i1] + 0.5 * lk[1] * b_Quu[i1 + 1]) * lk[i];
      }

      b_u_lims[0] = dV[0] + (lk[0] * Qu[0] + lk[1] * Qu[1]);
      b_u_lims[1] = dV[1] + result;
      for (i = 0; i < 2; i++) {
        dV[i] = b_u_lims[i];
        for (i1 = 0; i1 < 4; i1++) {
          a_tmp[i1 + (i << 2)] = Kk[i + (i1 << 1)];
        }
      }

      for (i = 0; i < 4; i++) {
        result = a_tmp[i + 4];
        for (i1 = 0; i1 < 2; i1++) {
          i2 = i1 << 1;
          Vx_tmp[i + (i1 << 2)] = a_tmp[i] * b_Quu[i2] + result * b_Quu[i2 + 1];
        }
      }

      for (i = 0; i < 2; i++) {
        for (i1 = 0; i1 < 4; i1++) {
          b_a_tmp[i1 + (i << 2)] = Qux[i + (i1 << 1)];
        }
      }

      for (i = 0; i < 4; i++) {
        result = 0.0;
        for (i1 = 0; i1 < 4; i1++) {
          result += fx[(i1 + (i << 2)) + ((499 - k) << 4)] * Vx[i1];
        }

        Quu[i] = cx[i + ((499 - k) << 2)] + result;
        Luu_data[i] = Vx_tmp[i] * lk[0] + Vx_tmp[i + 4] * lk[1];
      }

      for (i = 0; i < 4; i++) {
        result = (Quu[i] + Luu_data[i]) + (a_tmp[i] * Qu[0] + a_tmp[i + 4] * Qu
          [1]);
        d = b_a_tmp[i] * lk[0] + b_a_tmp[i + 4] * lk[1];
        Quu[i] = d;
        result += d;
        Vx[i] = result;
        for (i1 = 0; i1 < 4; i1++) {
          result = 0.0;
          for (i2 = 0; i2 < 4; i2++) {
            result += fx[(i2 + (i << 2)) + ((499 - k) << 4)] * Vxx[i2 + (i1 << 2)];
          }

          b_fx[i + (i1 << 2)] = result;
        }

        for (i1 = 0; i1 < 4; i1++) {
          result = 0.0;
          for (i2 = 0; i2 < 4; i2++) {
            result += b_fx[i + (i2 << 2)] * fx[(i2 + (i1 << 2)) + ((499 - k) <<
              4)];
          }

          Vxx_tmp = i + (i1 << 2);
          b_cxx[Vxx_tmp] = cxx[Vxx_tmp + ((499 - k) << 4)] + result;
        }

        result = Vx_tmp[i + 4];
        for (i1 = 0; i1 < 4; i1++) {
          i2 = i1 << 1;
          b_fx[i + (i1 << 2)] = Vx_tmp[i] * Kk[i2] + result * Kk[i2 + 1];
        }
      }

      for (i = 0; i < 4; i++) {
        result = a_tmp[i + 4];
        d = b_a_tmp[i + 4];
        for (i1 = 0; i1 < 4; i1++) {
          i2 = i1 << 1;
          partialTrueCount = i2 + 1;
          Vxx_tmp = i + (i1 << 2);
          Vxx[Vxx_tmp] = (b_cxx[Vxx_tmp] + b_fx[Vxx_tmp]) + (a_tmp[i] * Qux[i2]
            + result * Qux[partialTrueCount]);
          b_fx[Vxx_tmp] = b_a_tmp[i] * Kk[i2] + d * Kk[partialTrueCount];
        }
      }

      for (i = 0; i < 16; i++) {
        Vxx[i] += b_fx[i];
      }

      for (i = 0; i < 4; i++) {
        for (i1 = 0; i1 < 4; i1++) {
          Vxx_tmp = i1 + (i << 2);
          b_fx[Vxx_tmp] = 0.5 * (Vxx[Vxx_tmp] + Vxx[i + (i1 << 2)]);
        }
      }

      memcpy(&Vxx[0], &b_fx[0], 16U * sizeof(double));

      /*  Make sure Hessian is symmetric */
      /*  Update Control Vectors */
      l[u_lims_tmp] = -lk[0];
      l[b_u_lims_tmp] = -lk[1];
      for (i = 0; i < 4; i++) {
        Vxx_tmp = i << 1;
        partialTrueCount = Vxx_tmp + ((499 - k) << 3);
        K[partialTrueCount] = -Kk[Vxx_tmp];
        K[partialTrueCount + 1] = -Kk[Vxx_tmp + 1];
      }

      k++;
    }
  }
}




static void boxQPsolve(const double Quu[4], const double Qu[2], const double
  lower[2], const double upper[2], const double u0[2], double u[2], double
  *result, double Luu[4], bool b_free[2])
{
  bool clamped[2];
  int idxA1j;
  double oldvalue;
  double absxk;
  double ssq;
  double t;
  double value;
  int idxAjj;
  int iter;
  int b_iter;
  bool exitg1;
  bool factorize;
  double grad[2];
  bool prev_clamped[2];
  bool exitg2;
  bool guard1 = false;
  int trueCount;
  signed char tmp_data[2];
  signed char b_tmp_data[2];
  int Luu_size[2];
  double gnorm;
  int info;
  double deltaX_idx_0;
  int b_info;
  double Luu_data[4];
  int j;
  double deltaX_idx_1;
  signed char c_tmp_data[2];
  signed char d_tmp_data[2];
  int tmp_size[1];
  double b_Qu[2];
  double e_tmp_data[2];
  double f_tmp_data[4];
  double step;

  /*  Finds the optimal control with limits to minimize a quadratic cost */
  /*  Minimize 0.5*u'*Quu*u + u'*Qu  s.t. lower <= u <= upper */
  /*  */
  /*   inputs: */
  /*      Quu       - positive definite matrix              (n * n) */
  /*      Qu        - bias vector                           (n) */
  /*      lower     - lower bounds                          (n) */
  /*      upper     - upper bounds                          (n) */
  /*      u0        - initial control input for warm-start  (n) */
  /*  */
  /*   outputs: */
  /*      u         - solution                   (n) */
  /*      result    - result type (roughly, higher is better, see below) */
  /*      Luu       - cholesky factor            (n * n) */
  /*      free      - set of free dimensions     (n) */
  /*  Initialize arrays */
  clamped[0] = false;
  b_free[0] = true;
  clamped[1] = false;
  b_free[1] = true;
  for (idxA1j = 0; idxA1j < 4; idxA1j++) {
    Luu[idxA1j] = 0.0;
  }

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
  absxk = fmax(lower[0], fmin(upper[0], u0[0]));
  u[0] = absxk;
  ssq = absxk * Qu[0];
  absxk = fmax(lower[1], fmin(upper[1], u0[1]));
  u[1] = absxk;
  ssq += absxk * Qu[1];
  t = 0.0;
  for (idxA1j = 0; idxA1j < 2; idxA1j++) {
    idxAjj = idxA1j << 1;
    t += (0.5 * u[0] * Quu[idxAjj] + 0.5 * absxk * Quu[idxAjj + 1]) * u[idxA1j];
  }

  value = ssq + t;

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
        for (idxA1j = 0; idxA1j < 2; idxA1j++) {
          absxk = Qu[idxA1j] + (Quu[idxA1j] * u[0] + Quu[idxA1j + 2] * u[1]);
          grad[idxA1j] = absxk;
          prev_clamped[idxA1j] = clamped[idxA1j];
          factorize = false;
          clamped[idxA1j] = false;
          if ((u[idxA1j] == lower[idxA1j]) && (absxk > 0.0)) {
            factorize = true;
            clamped[idxA1j] = true;
          }

          if ((u[idxA1j] == upper[idxA1j]) && (absxk < 0.0)) {
            factorize = true;
            clamped[idxA1j] = true;
          }

          b_free[idxA1j] = !factorize;
        }

        /*  check for all clamped */
        factorize = true;
        idxA1j = 0;
        exitg2 = false;
        while ((!exitg2) && (idxA1j < 2)) {
          if (!clamped[idxA1j]) {
            factorize = false;
            exitg2 = true;
          } else {
            idxA1j++;
          }
        }

        if (factorize) {
          *result = 5.0;
          exitg1 = true;
        } else {
          /*  Cholesky factorize if clamped controls have changed */
          if (b_iter == 1) {
            factorize = true;
          } else {
            factorize = false;
            idxA1j = 0;
            exitg2 = false;
            while ((!exitg2) && (idxA1j < 2)) {
              if (prev_clamped[idxA1j] != clamped[idxA1j]) {
                factorize = true;
                exitg2 = true;
              } else {
                idxA1j++;
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

            idxA1j = 0;
            if (b_free[0]) {
              tmp_data[0] = 1;
              idxA1j = 1;
            }

            if (b_free[1]) {
              tmp_data[idxA1j] = 2;
            }

            /*  Wrapper for MATLAB chol for use with auto coder */
            /*  Inputs: */
            /* =========== */
            /*  A - positive semi-definite matrix */
            Luu_size[0] = trueCount;
            Luu_size[1] = trueCount;
            for (idxA1j = 0; idxA1j < trueCount; idxA1j++) {
              for (idxAjj = 0; idxAjj < trueCount; idxAjj++) {
                Luu_data[idxAjj + trueCount * idxA1j] = Quu[(tmp_data[idxAjj] +
                  ((tmp_data[idxA1j] - 1) << 1)) - 1];
              }
            }

            info = 0;
            if (trueCount != 0) {
              b_info = -1;
              j = 0;
              exitg2 = false;
              while ((!exitg2) && (j <= trueCount - 1)) {
                idxA1j = j * trueCount;
                idxAjj = idxA1j + j;
                ssq = 0.0;
                if (j >= 1) {
                  ssq = Luu_data[idxA1j] * Luu_data[idxA1j];
                }

                ssq = Luu_data[idxAjj] - ssq;
                if (ssq > 0.0) {
                  ssq = sqrt(ssq);
                  Luu_data[idxAjj] = ssq;
                  if (j + 1 < trueCount) {
                    idxAjj += 3;
                    ssq = 1.0 / ssq;
                    for (idxA1j = idxAjj; idxA1j <= idxAjj; idxA1j += 2) {
                      Luu_data[idxA1j - 1] *= ssq;
                    }
                  }

                  j++;
                } else {
                  Luu_data[idxAjj] = ssq;
                  b_info = j;
                  exitg2 = true;
                }
              }

              info = b_info + 1;
              if (b_info + 1 == 0) {
                idxA1j = trueCount;
              } else {
                idxA1j = b_info;
              }

              for (j = 0; j < idxA1j; j++) {
                if (j + 2 <= idxA1j) {
                  Luu_data[trueCount * j + 1] = 0.0;
                }
              }

              if (1 > idxA1j) {
                j = 0;
                b_info = 0;
              } else {
                j = idxA1j;
                b_info = idxA1j;
              }

              for (idxA1j = 0; idxA1j < b_info; idxA1j++) {
                for (idxAjj = 0; idxAjj < j; idxAjj++) {
                  f_tmp_data[idxAjj + j * idxA1j] = Luu_data[idxAjj + trueCount *
                    idxA1j];
                }
              }

              Luu_size[0] = j;
              Luu_size[1] = b_info;
              j *= b_info;
              if (0 <= j - 1) {
                memcpy(&Luu_data[0], &f_tmp_data[0], j * sizeof(double));
              }
            }

            idxA1j = 0;
            if (b_free[0]) {
              c_tmp_data[0] = 1;
              idxA1j = 1;
            }

            if (b_free[1]) {
              c_tmp_data[idxA1j] = 2;
            }

            j = Luu_size[1];
            for (idxA1j = 0; idxA1j < j; idxA1j++) {
              b_info = Luu_size[0];
              for (idxAjj = 0; idxAjj < b_info; idxAjj++) {
                Luu[(c_tmp_data[idxAjj] + ((c_tmp_data[idxA1j] - 1) << 1)) - 1] =
                  Luu_data[idxAjj + Luu_size[0] * idxA1j];
              }
            }

            if (info != 0) {
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

            idxA1j = 0;
            if (b_free[0]) {
              b_tmp_data[0] = 1;
              idxA1j = 1;
            }

            if (b_free[1]) {
              b_tmp_data[idxA1j] = 2;
            }

            if (trueCount == 0) {
              gnorm = 0.0;
            } else if (trueCount == 1) {
              gnorm = fabs(grad[b_tmp_data[0] - 1]);
            } else {
              ssq = 3.3121686421112381E-170;
              absxk = fabs(grad[b_tmp_data[0] - 1]);
              if (absxk > 3.3121686421112381E-170) {
                gnorm = 1.0;
                ssq = absxk;
              } else {
                t = absxk / 3.3121686421112381E-170;
                gnorm = t * t;
              }

              absxk = fabs(grad[b_tmp_data[1] - 1]);
              if (absxk > ssq) {
                t = ssq / absxk;
                gnorm = gnorm * t * t + 1.0;
                ssq = absxk;
              } else {
                t = absxk / ssq;
                gnorm += t * t;
              }

              gnorm = ssq * sqrt(gnorm);
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

              idxA1j = 0;
              if (b_free[0]) {
                d_tmp_data[0] = 1;
                idxA1j = 1;
              }

              ssq = u[0] * (double)clamped[0];
              if (b_free[1]) {
                d_tmp_data[idxA1j] = 2;
              }

              absxk = u[1] * (double)clamped[1];
              for (idxA1j = 0; idxA1j < 2; idxA1j++) {
                b_Qu[idxA1j] = Qu[idxA1j] + (Quu[idxA1j] * ssq + Quu[idxA1j + 2]
                  * absxk);
              }

              tmp_size[0] = trueCount;
              for (idxA1j = 0; idxA1j < trueCount; idxA1j++) {
                e_tmp_data[idxA1j] = b_Qu[d_tmp_data[idxA1j] - 1];
              }

              Luu_size[0] = trueCount;
              Luu_size[1] = trueCount;
              for (idxA1j = 0; idxA1j < trueCount; idxA1j++) {
                for (idxAjj = 0; idxAjj < trueCount; idxAjj++) {
                  Luu_data[idxAjj + trueCount * idxA1j] = Luu[(d_tmp_data[idxA1j]
                    + ((d_tmp_data[idxAjj] - 1) << 1)) - 1];
                }
              }

              mldivide(Luu_data, Luu_size, e_tmp_data, tmp_size);
              Luu_size[0] = trueCount;
              Luu_size[1] = trueCount;
              for (idxA1j = 0; idxA1j < trueCount; idxA1j++) {
                for (idxAjj = 0; idxAjj < trueCount; idxAjj++) {
                  Luu_data[idxAjj + trueCount * idxA1j] = -Luu
                    [(d_tmp_data[idxAjj] + ((d_tmp_data[idxA1j] - 1) << 1)) - 1];
                }
              }

              mldivide(Luu_data, Luu_size, e_tmp_data, tmp_size);
              idxA1j = 0;

              /*  cholesky solver */
              /*  check for descent direction */
              if (b_free[0]) {
                deltaX_idx_0 = e_tmp_data[0] - u[0];
                idxA1j = 1;
              }

              grad[0] *= deltaX_idx_0;
              if (b_free[1]) {
                deltaX_idx_1 = e_tmp_data[idxA1j] - u[1];
              }

              grad[1] *= deltaX_idx_1;
              gnorm = grad[0] + grad[1];
              if (gnorm >= 0.0) {
                /*  (should not happen) */
                exitg1 = true;
              } else {
                /*  Armijo linesearch */
                step = 1.0;

                /*  Returns array x with all values clamped between lower and upper */
                absxk = fmax(lower[0], fmin(upper[0], u[0] + deltaX_idx_0));
                grad[0] = absxk;
                ssq = absxk * Qu[0];
                absxk = fmax(lower[1], fmin(upper[1], u[1] + deltaX_idx_1));
                grad[1] = absxk;
                ssq += absxk * Qu[1];
                t = 0.0;
                for (idxA1j = 0; idxA1j < 2; idxA1j++) {
                  idxAjj = idxA1j << 1;
                  t += (0.5 * grad[0] * Quu[idxAjj] + 0.5 * absxk * Quu[idxAjj +
                        1]) * grad[idxA1j];
                }

                ssq += t;
                exitg2 = false;
                while ((!exitg2) && ((ssq - value) / (step * gnorm) < 0.1)) {
                  step *= 0.6;

                  /*  Returns array x with all values clamped between lower and upper */
                  absxk = fmax(lower[0], fmin(upper[0], u[0] + step *
                    deltaX_idx_0));
                  grad[0] = absxk;
                  ssq = absxk * Qu[0];
                  absxk = fmax(lower[1], fmin(upper[1], u[1] + step *
                    deltaX_idx_1));
                  grad[1] = absxk;
                  ssq += absxk * Qu[1];
                  t = 0.0;
                  for (idxA1j = 0; idxA1j < 2; idxA1j++) {
                    idxAjj = idxA1j << 1;
                    t += (0.5 * grad[0] * Quu[idxAjj] + 0.5 * absxk * Quu[idxAjj
                          + 1]) * grad[idxA1j];
                  }

                  ssq += t;
                  if (step < 1.0E-22) {
                    *result = 2.0;
                    exitg2 = true;
                  }
                }

                /*  accept candidate */
                u[0] = grad[0];
                u[1] = grad[1];
                value = ssq;
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




static void car_cost(const double x[4], const double xg[4], const double u[2],
                     double terminal, double *cost, double cx[4], double cu[2],
                     double cxx[16])
{
  static const double b_dv[16] = { 1.0E-5, 0.0, 0.0, 0.0, 0.0, 1.0E-5, 0.0, 0.0,
    0.0, 0.0, 1.0E-5, 0.0, 0.0, 0.0, 0.0, 1.0E-5 };

  double Qf[16];
  static const double dv1[16] = { 0.3, 0.0, 0.0, 0.0, 0.0, 0.3, 0.0, 0.0, 0.0,
    0.0, 0.3, 0.0, 0.0, 0.0, 0.0, 0.3 };

  int i;
  double d;
  double b[4];
  double d1;
  int b_i;
  static const double b_b[16] = { 0.001, 0.0, 0.0, 0.0, 0.0, 0.001, 0.0, 0.0,
    0.0, 0.0, 1.0E-5, 0.0, 0.0, 0.0, 0.0, 1.0E-5 };

  static const double c_b[16] = { 0.3, 0.0, 0.0, 0.0, 0.0, 0.3, 0.0, 0.0, 0.0,
    0.0, 1.5, 0.0, 0.0, 0.0, 0.0, 0.3 };

  /*  Calculates the cost component from the input state/control */
  /*  Also calculates the cost derivates */
  /*  Cost matrices */
  /*  Cumulative */
  memcpy(&cxx[0], &b_dv[0], 16U * sizeof(double));
  memcpy(&Qf[0], &dv1[0], 16U * sizeof(double));
  cxx[0] = 0.001;
  cxx[5] = 0.001;

  /*  Terminal */
  Qf[10] = 1.5;
  if (terminal != 0.0) {
    /*  Terminal cost */
    for (i = 0; i < 4; i++) {
      b[i] = x[i] - xg[i];
    }

    *cost = 0.0;
    for (i = 0; i < 4; i++) {
      d = 0.0;
      for (b_i = 0; b_i < 4; b_i++) {
        d += 0.5 * b[b_i] * c_b[b_i + (i << 2)];
      }

      *cost += d * b[i];
    }

    /*  Final cost gradients */
    memcpy(&cxx[0], &Qf[0], 16U * sizeof(double));
    for (i = 0; i < 4; i++) {
      d = 0.0;
      for (b_i = 0; b_i < 4; b_i++) {
        d += c_b[i + (b_i << 2)] * b[b_i];
      }

      cx[i] = d;
    }
  } else {
    /*  Cumulative cost */
    for (i = 0; i < 4; i++) {
      b[i] = x[i] - xg[i];
    }

    d = 0.0;
    for (i = 0; i < 4; i++) {
      d1 = 0.0;
      for (b_i = 0; b_i < 4; b_i++) {
        d1 += 0.5 * b[b_i] * b_b[b_i + (i << 2)];
      }

      d += d1 * b[i];
    }

    d1 = 0.0;
    for (i = 0; i < 2; i++) {
      b_i = i << 1;
      d1 += (0.5 * u[0] * dv[b_i] + 0.5 * u[1] * dv[b_i + 1]) * u[i];
    }

    *cost = d + d1;

    /*  Cost gradients */
    for (i = 0; i < 4; i++) {
      d = 0.0;
      for (b_i = 0; b_i < 4; b_i++) {
        d += b_b[i + (b_i << 2)] * b[b_i];
      }

      cx[i] = d;
    }
  }

  /*  Control cost gradients */
  for (i = 0; i < 2; i++) {
    cu[i] = dv[i] * u[0] + dv[i + 2] * u[1];
  }
}




static void forwardRollout(const double x[2004], const double xg[4], const
  double u[1000], const double l[1000], const double K[4000], double alpha,
  const double u_lims[4], double xnew[2004], double unew[1000], double fx[8000],
  double fu[4000], double cx[2004], double cu[1000], double cxx[8016], double
  cuu[2000], double *cost)
{
  int i;
  int k;
  double b_dv[2];
  double A_tmp;
  double z1[2];
  int xnew_tmp;
  int b_k;
  double x0[4];
  int th_tmp_tmp;
  double th;
  int xdot_tmp;
  double xdot[4];
  double xdot_tmp_tmp;
  double xdot_tmp_tmp_tmp;
  double b_xdot_tmp;
  double c_xdot_tmp;
  double A[16];
  double B[8];
  signed char b_I[16];

  /*  Uses an rk method to roll out a trajectory */
  /*  Returns the new trajectory, cost and the derivatives along the trajectory */
  /*  If feed-forward and feed-back controls l, K are non-zero */
  /*  then the unew returned is the new control sequeunce. Otherwise unew = u */
  /*  Sizes */
  /*  Initialize outputs */
  memset(&xnew[0], 0, 2004U * sizeof(double));
  memset(&unew[0], 0, 1000U * sizeof(double));
  memset(&cx[0], 0, 2004U * sizeof(double));
  memset(&cxx[0], 0, 8016U * sizeof(double));
  *cost = 0.0;
  for (i = 0; i < 4; i++) {
    xnew[i] = x[i];
  }

  for (k = 0; k < 500; k++) {
    /*  Update the control during line-search */
    /*  (During inital forward rollout, l and K will be arrays of zeros) */
    for (i = 0; i < 4; i++) {
      xnew_tmp = i + (k << 2);
      x0[i] = xnew[xnew_tmp] - x[xnew_tmp];
    }

    /*  Ensure control is within limits */
    for (b_k = 0; b_k < 2; b_k++) {
      A_tmp = 0.0;
      for (i = 0; i < 4; i++) {
        A_tmp += K[(b_k + (i << 1)) + (k << 3)] * x0[i];
      }

      xnew_tmp = b_k + (k << 1);
      unew[xnew_tmp] = fmin(u_lims[b_k + 2], fmax(u_lims[b_k], (u[xnew_tmp] -
        alpha * l[xnew_tmp]) - A_tmp));
    }

    /*  Step the dynamics forward */
    /*  Forward Euler step for car */
    /*  returns x, fx, fu */
    /*  === states and controls: */
    /*  x = [x y th v]' = [x; y; car_angle; front_wheel_velocity] */
    /*  u = [delta a]' = [front_wheel_angle; acceleration] */
    /*  constants */
    /*  wheelbase */
    /*  states */
    th_tmp_tmp = k << 2;
    xnew_tmp = th_tmp_tmp + 2;
    th = xnew[xnew_tmp];

    /*  controls */
    /*  front wheel angle */
    xdot_tmp = k << 1;
    xdot[3] = unew[xdot_tmp + 1];

    /*  front wheel acceleration */
    xdot_tmp_tmp = cos(unew[xdot_tmp]);
    xdot_tmp_tmp_tmp = cos(xnew[xnew_tmp]);
    b_k = th_tmp_tmp + 3;
    xdot[0] = xnew[b_k] * xdot_tmp_tmp * xdot_tmp_tmp_tmp;

    /*  xdot */
    b_xdot_tmp = sin(xnew[xnew_tmp]);
    xdot[1] = xnew[b_k] * cos(unew[xdot_tmp]) * b_xdot_tmp;

    /*  ydot */
    c_xdot_tmp = sin(unew[xdot_tmp]);
    xdot[2] = xnew[b_k] * c_xdot_tmp / 2.0;

    /*  thetadot */
    /*  vdot */
    /*  Euler step */
    for (i = 0; i < 4; i++) {
      A_tmp = xnew[i + th_tmp_tmp];
      x0[i] = A_tmp;
      xnew[i + ((k + 1) << 2)] = A_tmp + xdot[i] * 0.03;
    }

    /*  Jacobians */
    memset(&A[0], 0, 16U * sizeof(double));
    A[8] = -x0[3] * xdot_tmp_tmp * b_xdot_tmp;
    A[12] = xdot_tmp_tmp * xdot_tmp_tmp_tmp;
    A_tmp = x0[3] * xdot_tmp_tmp;
    A[9] = A_tmp * cos(th);
    A[13] = xdot_tmp_tmp * b_xdot_tmp;
    A[14] = c_xdot_tmp / 2.0;
    memset(&B[0], 0, 8U * sizeof(double));
    B[0] = -x0[3] * c_xdot_tmp * xdot_tmp_tmp_tmp;
    B[1] = -x0[3] * sin(unew[xdot_tmp]) * b_xdot_tmp;
    B[2] = A_tmp / 2.0;
    B[7] = 1.0;
    for (i = 0; i < 16; i++) {
      b_I[i] = 0;
    }

    for (b_k = 0; b_k < 4; b_k++) {
      b_I[b_k + (b_k << 2)] = 1;
    }

    for (i = 0; i < 4; i++) {
      for (b_k = 0; b_k < 4; b_k++) {
        xnew_tmp = b_k + (i << 2);
        fx[xnew_tmp + (k << 4)] = (double)b_I[xnew_tmp] + 0.03 * A[xnew_tmp];
      }
    }

    /*  Calculate the cost */
    car_cost(*(double (*)[4])&xnew[k << 2], xg, *(double (*)[2])&unew[k << 1],
             0.0, &A_tmp, *(double (*)[4])&cx[k << 2], *(double (*)[2])&cu[k <<
             1], *(double (*)[16])&cxx[k << 4]);
    for (i = 0; i < 2; i++) {
      for (b_k = 0; b_k < 4; b_k++) {
        xnew_tmp = b_k + (i << 2);
        fu[xnew_tmp + (k << 3)] = 0.03 * B[xnew_tmp];
      }

      xnew_tmp = i << 1;
      b_k = xnew_tmp + th_tmp_tmp;
      cuu[b_k] = dv[xnew_tmp];
      cuu[b_k + 1] = dv[xnew_tmp + 1];
    }

    *cost += A_tmp;
  }

  /*  Final cost */
  b_dv[0] = 0.0;
  b_dv[1] = 0.0;
  car_cost(*(double (*)[4])&xnew[2000], xg, b_dv, 1.0, &A_tmp, *(double (*)[4])&
           cx[2000], z1, *(double (*)[16])&cxx[8000]);
  *cost += A_tmp;
}




static void mldivide(const double A_data[], const int A_size[2], double B_data[],
                     int B_size[1])
{
  int b_A_size[2];
  int mn;
  int loop_ub;
  double b_A_data[4];
  double tau_data[2];
  int tau_size[1];
  int jpvt_data[2];
  int jpvt_size[2];
  int rankA;
  int i;
  double b_B_data[2];
  double wj;
  int m;
  int b_i;
  int j;
  if ((A_size[0] == 0) || (A_size[1] == 0) || (B_size[0] == 0)) {
    B_size[0] = (signed char)A_size[1];
    loop_ub = (signed char)A_size[1];
    if (0 <= loop_ub - 1) {
      memset(&B_data[0], 0, loop_ub * sizeof(double));
    }
  } else if (A_size[0] == A_size[1]) {
    mn = A_size[1];
    b_A_size[0] = A_size[0];
    b_A_size[1] = A_size[1];
    loop_ub = A_size[0] * A_size[1];
    if (0 <= loop_ub - 1) {
      memcpy(&b_A_data[0], &A_data[0], loop_ub * sizeof(double));
    }

    xgetrf(A_size[1], A_size[1], b_A_data, b_A_size, A_size[1], jpvt_data,
           jpvt_size, &loop_ub);
    i = A_size[1];
    for (loop_ub = 0; loop_ub <= i - 2; loop_ub++) {
      if (jpvt_data[0] != 1) {
        wj = B_data[0];
        B_data[0] = B_data[jpvt_data[0] - 1];
        B_data[jpvt_data[0] - 1] = wj;
      }
    }

    for (loop_ub = 0; loop_ub < mn; loop_ub++) {
      m = mn * loop_ub;
      if (B_data[loop_ub] != 0.0) {
        i = loop_ub + 2;
        for (b_i = i; b_i <= mn; b_i++) {
          B_data[1] -= B_data[loop_ub] * b_A_data[m + 1];
        }
      }
    }

    for (loop_ub = mn; loop_ub >= 1; loop_ub--) {
      m = mn * (loop_ub - 1);
      wj = B_data[loop_ub - 1];
      if (wj != 0.0) {
        B_data[loop_ub - 1] = wj / b_A_data[(loop_ub + m) - 1];
        for (b_i = 0; b_i <= loop_ub - 2; b_i++) {
          B_data[0] -= B_data[loop_ub - 1] * b_A_data[m];
        }
      }
    }
  } else {
    b_A_size[0] = A_size[0];
    b_A_size[1] = A_size[1];
    loop_ub = A_size[0] * A_size[1];
    if (0 <= loop_ub - 1) {
      memcpy(&b_A_data[0], &A_data[0], loop_ub * sizeof(double));
    }

    xgeqp3(b_A_data, b_A_size, tau_data, tau_size, jpvt_data, jpvt_size);
    rankA = rankFromQR(b_A_data, b_A_size);
    loop_ub = B_size[0];
    if (0 <= loop_ub - 1) {
      memcpy(&b_B_data[0], &B_data[0], loop_ub * sizeof(double));
    }

    B_size[0] = (signed char)b_A_size[1];
    loop_ub = (signed char)b_A_size[1];
    if (0 <= loop_ub - 1) {
      memset(&B_data[0], 0, loop_ub * sizeof(double));
    }

    m = b_A_size[0];
    loop_ub = b_A_size[0];
    mn = b_A_size[1];
    if (loop_ub < mn) {
      mn = loop_ub;
    }

    for (j = 0; j < mn; j++) {
      if (tau_data[j] != 0.0) {
        wj = b_B_data[j];
        i = j + 2;
        for (b_i = i; b_i <= m; b_i++) {
          wj += b_A_data[b_A_size[0] * j + 1] * b_B_data[1];
        }

        wj *= tau_data[j];
        if (wj != 0.0) {
          b_B_data[j] -= wj;
          i = j + 2;
          for (b_i = i; b_i <= m; b_i++) {
            b_B_data[1] -= b_A_data[b_A_size[0] * j + 1] * wj;
          }
        }
      }
    }

    for (b_i = 0; b_i < rankA; b_i++) {
      B_data[jpvt_data[b_i] - 1] = b_B_data[b_i];
    }

    for (j = rankA; j >= 1; j--) {
      i = jpvt_data[j - 1];
      loop_ub = b_A_size[0] * (j - 1);
      B_data[i - 1] /= b_A_data[(j + loop_ub) - 1];
      for (b_i = 0; b_i <= j - 2; b_i++) {
        B_data[jpvt_data[0] - 1] -= B_data[jpvt_data[j - 1] - 1] *
          b_A_data[loop_ub];
      }
    }
  }
}




static void qrpf(double A_data[], const int A_size[2], int m, int n, double
                 tau_data[], int jpvt_data[])
{
  int ma;
  int minmn;
  int itemp;
  double work_data[2];
  double vn1_data[2];
  double vn2_data[2];
  int j;
  int i;
  double temp;
  int ip1;
  int iy;
  int ii;
  int nmi;
  int mmi;
  int pvt;
  int ix;
  double atmp;
  int jA;
  double beta1;
  int lastv;
  int b_i;
  bool exitg2;
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
    temp = xnrm2(m, A_data, j * ma + 1);
    vn1_data[j] = temp;
    vn2_data[j] = temp;
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
      if ((nmi > 1) && (fabs(vn1_data[i + 1]) > fabs(vn1_data[i]))) {
        itemp = 1;
      }
    }

    pvt = i + itemp;
    if (pvt + 1 != i + 1) {
      ix = pvt * ma;
      for (j = 0; j < m; j++) {
        temp = A_data[ix];
        A_data[ix] = A_data[iy];
        A_data[iy] = temp;
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
      atmp = A_data[ii];
      itemp = ii + 2;
      tau_data[0] = 0.0;
      if (mmi > 0) {
        temp = xnrm2(mmi - 1, A_data, ii + 2);
        if (temp != 0.0) {
          beta1 = rt_hypotd(A_data[ii], temp);
          if (A_data[ii] >= 0.0) {
            beta1 = -beta1;
          }

          if (fabs(beta1) < 1.0020841800044864E-292) {
            pvt = -1;
            b_i = ii + mmi;
            do {
              pvt++;
              for (j = itemp; j <= b_i; j++) {
                A_data[j - 1] *= 9.9792015476736E+291;
              }

              beta1 *= 9.9792015476736E+291;
              atmp *= 9.9792015476736E+291;
            } while (!(fabs(beta1) >= 1.0020841800044864E-292));

            beta1 = rt_hypotd(atmp, xnrm2(mmi - 1, A_data, ii + 2));
            if (atmp >= 0.0) {
              beta1 = -beta1;
            }

            tau_data[0] = (beta1 - atmp) / beta1;
            temp = 1.0 / (atmp - beta1);
            for (j = itemp; j <= b_i; j++) {
              A_data[j - 1] *= temp;
            }

            for (j = 0; j <= pvt; j++) {
              beta1 *= 1.0020841800044864E-292;
            }

            atmp = beta1;
          } else {
            tau_data[0] = (beta1 - A_data[ii]) / beta1;
            temp = 1.0 / (A_data[ii] - beta1);
            b_i = ii + mmi;
            for (j = itemp; j <= b_i; j++) {
              A_data[j - 1] *= temp;
            }

            atmp = beta1;
          }
        }
      }

      A_data[ii] = atmp;
    } else {
      tau_data[i] = 0.0;
    }

    if (i + 1 < n) {
      atmp = A_data[ii];
      A_data[ii] = 1.0;
      jA = (ii + ma) + 1;
      if (tau_data[0] != 0.0) {
        lastv = mmi - 1;
        itemp = (ii + mmi) - 1;
        while ((lastv + 1 > 0) && (A_data[itemp] == 0.0)) {
          lastv--;
          itemp--;
        }

        nmi -= 2;
        exitg2 = false;
        while ((!exitg2) && (nmi + 1 > 0)) {
          itemp = jA;
          do {
            exitg1 = 0;
            if (itemp <= jA + lastv) {
              if (A_data[itemp - 1] != 0.0) {
                exitg1 = 1;
              } else {
                itemp++;
              }
            } else {
              nmi = -1;
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
            work_data[0] = 0.0;
          }

          iy = 0;
          b_i = jA + ma * nmi;
          for (pvt = jA; ma < 0 ? pvt >= b_i : pvt <= b_i; pvt += ma) {
            ix = ii;
            temp = 0.0;
            j = pvt + lastv;
            for (itemp = pvt; itemp <= j; itemp++) {
              temp += A_data[itemp - 1] * A_data[ix];
              ix++;
            }

            work_data[iy] += temp;
            iy++;
          }
        }

        if (-tau_data[0] != 0.0) {
          itemp = 0;
          for (j = 0; j <= nmi; j++) {
            if (work_data[itemp] != 0.0) {
              temp = work_data[itemp] * -tau_data[0];
              ix = ii;
              b_i = lastv + jA;
              for (pvt = jA; pvt <= b_i; pvt++) {
                A_data[pvt - 1] += A_data[ix] * temp;
                ix++;
              }
            }

            itemp++;
            jA += ma;
          }
        }
      }

      A_data[ii] = atmp;
    }

    for (j = ip1; j <= n; j++) {
      itemp = i + ma;
      if (vn1_data[1] != 0.0) {
        temp = fabs(A_data[itemp]) / vn1_data[1];
        temp = 1.0 - temp * temp;
        if (temp < 0.0) {
          temp = 0.0;
        }

        beta1 = vn1_data[1] / vn2_data[1];
        beta1 = temp * (beta1 * beta1);
        if (beta1 <= 1.4901161193847656E-8) {
          if (i + 1 < m) {
            temp = xnrm2(mmi - 1, A_data, itemp + 2);
            vn1_data[1] = temp;
            vn2_data[1] = temp;
          } else {
            vn1_data[1] = 0.0;
            vn2_data[1] = 0.0;
          }
        } else {
          vn1_data[1] *= sqrt(temp);
        }
      }
    }
  }
}




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
    while ((r < minmn) && (fabs(A_data[r + A_size[0] * r]) > tol)) {
      r++;
    }
  }

  return r;
}



static double rt_hypotd(double u0, double u1)
{
  double y;
  double a;
  double b;
  a = fabs(u0);
  b = fabs(u1);
  if (a < b) {
    a /= b;
    y = b * sqrt(a * a + 1.0);
  } else if (a > b) {
    b /= a;
    y = a * sqrt(b * b + 1.0);
  } else {
    y = a * 1.4142135623730951;
  }

  return y;
}




static void xgeqp3(double A_data[], const int A_size[2], double tau_data[], int
                   tau_size[1], int jpvt_data[], int jpvt_size[2])
{
  int n;
  int u0;
  int u1;
  n = A_size[1] - 1;
  u0 = A_size[0];
  u1 = A_size[1];
  if (u0 < u1) {
    u1 = u0;
  }

  tau_size[0] = u1;
  if (0 <= u1 - 1) {
    memset(&tau_data[0], 0, u1 * sizeof(double));
  }

  if ((A_size[0] == 0) || (A_size[1] == 0) || (u1 < 1)) {
    jpvt_size[0] = 1;
    jpvt_size[1] = A_size[1];
    u0 = A_size[1];
    if (0 <= u0 - 1) {
      memset(&jpvt_data[0], 0, u0 * sizeof(int));
    }

    for (u0 = 0; u0 <= n; u0++) {
      jpvt_data[u0] = u0 + 1;
    }
  } else {
    jpvt_size[0] = 1;
    jpvt_size[1] = A_size[1];
    u0 = A_size[1];
    if (0 <= u0 - 1) {
      memset(&jpvt_data[0], 0, u0 * sizeof(int));
    }

    for (u0 = 0; u0 <= n; u0++) {
      jpvt_data[u0] = u0 + 1;
    }

    qrpf(A_data, A_size, A_size[0], A_size[1], tau_data, jpvt_data);
  }
}



static void xgetrf(int m, int n, double A_data[], const int A_size[2], int lda,
                   int ipiv_data[], int ipiv_size[2], int *info)
{
  int b_n;
  int yk;
  int k;
  int j;
  int ix;
  double temp;
  int i;
  int i1;
  int ijA;
  if (m < n) {
    b_n = m;
  } else {
    b_n = n;
  }

  if (b_n < 1) {
    b_n = 0;
  }

  ipiv_size[0] = 1;
  ipiv_size[1] = b_n;
  if (b_n > 0) {
    ipiv_data[0] = 1;
    yk = 1;
    for (k = 2; k <= b_n; k++) {
      yk++;
      ipiv_data[1] = yk;
    }
  }

  *info = 0;
  if ((m >= 1) && (n >= 1)) {
    for (j = 0; j <= m - 2; j++) {
      b_n = 0;
      if ((m > 1) && (fabs(A_data[1]) > fabs(A_data[0]))) {
        b_n = 1;
      }

      if (A_data[b_n] != 0.0) {
        if (b_n != 0) {
          ipiv_data[0] = 2;
          ix = 0;
          b_n = 1;
          for (k = 0; k < n; k++) {
            temp = A_data[ix];
            A_data[ix] = A_data[b_n];
            A_data[b_n] = temp;
            ix += lda;
            b_n += lda;
          }
        }

        for (b_n = 2; b_n <= m; b_n++) {
          A_data[1] /= A_data[0];
        }
      } else {
        *info = 1;
      }

      b_n = lda;
      yk = lda;
      for (k = 0; k <= n - 2; k++) {
        temp = A_data[b_n];
        if (A_data[b_n] != 0.0) {
          ix = 1;
          i = yk + 2;
          i1 = m + yk;
          for (ijA = i; ijA <= i1; ijA++) {
            A_data[ijA - 1] += A_data[ix] * -temp;
            ix++;
          }
        }

        b_n += lda;
        yk += lda;
      }
    }

    if ((*info == 0) && (m <= n) && (A_data[(m + A_size[0] * (m - 1)) - 1] ==
         0.0)) {
      *info = m;
    }
  }
}



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
      kend = ix0 + 1;
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



void ilqrCar(const double x0[2004], const double xg[4], const double u0[1000],
             const double u_lims[4], double x[2004], double u[1000], double K
             [4000], bool *result)
{
  int k;
  double Alphas[11];
  double lambda;
  double dlambda;
  static double x_n[2004];
  static double u_n[1000];
  static double fx_n[8000];
  static double fu_n[4000];
  static double cx_n[2004];
  static double cu_n[1000];
  static double cxx_n[8016];
  static double cuu_n[2000];
  double cost_n;
  static double l[1000];
  double dV[2];
  double y[1000];
  static double b_dv[4000];
  static double fx[8000];
  static double fu[4000];
  static double cx[2004];
  double cu[1000];
  static double cxx[8016];
  static double cuu[2000];
  double cost;
  int iter;
  int b_iter;
  bool exitg1;
  bool backPassDone;
  int exitg2;
  bool fwdPassDone;
  double expectedChange;
  double maxval[500];
  int i;
  bool exitg3;

  /*  Solves finite horizon optimal control problem using the iterative */
  /*  linear quadratic redualtor method */
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
  /*  */
  /*  Ops (options): */
  /*  ----------------------- */
  /*  dt, max_iters, exit_tol, grad_tol, lambda_tol, z_min, lambda_max, lambda_min, */
  /*  lambda_scaling */
  /*  Outputs */
  /*  =========================================== */
  /*  x - Final nominal trajectory (n, N) */
  /*  */
  /*  u - Final open-loop controls (m, N-1) */
  /*  */
  /*  K - Feedback control gains (n, m, N-1) */
  /*  Options (pass in as array) */
  /*  Timestep (Should match MCU Hz) */
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

  /*  Init matrices for update */
  memset(&x_n[0], 0, 2004U * sizeof(double));
  memset(&u_n[0], 0, 1000U * sizeof(double));
  memset(&fx_n[0], 0, 8000U * sizeof(double));
  memset(&fu_n[0], 0, 4000U * sizeof(double));
  memset(&cx_n[0], 0, 2004U * sizeof(double));
  memset(&cu_n[0], 0, 1000U * sizeof(double));
  memset(&cxx_n[0], 0, 8016U * sizeof(double));
  memset(&cuu_n[0], 0, 2000U * sizeof(double));
  cost_n = 0.0;

  /*  Initial Forward rollout */
  /*  Returns xtraj, (unew=utraj0), cost */
  memset(&l[0], 0, 1000U * sizeof(double));
  memset(&K[0], 0, 4000U * sizeof(double));
  dV[0] = 0.0;
  dV[1] = 0.0;
  memset(&y[0], 0, 1000U * sizeof(double));
  memset(&b_dv[0], 0, 4000U * sizeof(double));
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
                     &fwdPassDone);
        if (fwdPassDone) {
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
    for (k = 0; k < 1000; k++) {
      y[k] = fabs(l[k]) / (fabs(u[k]) + 1.0);
    }

    for (k = 0; k < 500; k++) {
      i = k << 1;
      maxval[k] = y[i];
      expectedChange = y[i + 1];
      if (y[i] < expectedChange) {
        maxval[k] = expectedChange;
      }
    }

    expectedChange = maxval[0];
    for (k = 0; k < 499; k++) {
      expectedChange += maxval[k + 1];
    }

    /*  Avg of max grad at each time step */
    if ((expectedChange / 500.0 < 0.0001) && (lambda < 1.0E-5)) {
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
            expectedChange = (cost - cost_n) / expectedChange;
          } else {
            expectedChange = cost - cost_n;
            if (expectedChange < 0.0) {
              expectedChange = -1.0;
            } else {
              if (expectedChange > 0.0) {
                expectedChange = 1.0;
              }
            }
          }

          if (expectedChange > 0.0) {
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
        expectedChange = cost - cost_n;

        /*  Update trajectory and controls */
        memcpy(&x[0], &x_n[0], 2004U * sizeof(double));
        memcpy(&u[0], &u_n[0], 1000U * sizeof(double));
        memcpy(&fx[0], &fx_n[0], 8000U * sizeof(double));
        memcpy(&fu[0], &fu_n[0], 4000U * sizeof(double));
        memcpy(&cx[0], &cx_n[0], 2004U * sizeof(double));
        memcpy(&cu[0], &cu_n[0], 1000U * sizeof(double));
        memcpy(&cxx[0], &cxx_n[0], 8016U * sizeof(double));
        memcpy(&cuu[0], &cuu_n[0], 2000U * sizeof(double));
        cost = cost_n;

        /*  Terminate ? */
        if (expectedChange < 1.0E-7) {
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
 * File trailer for ilqrCar.c
 *
 * [EOF]
 */
