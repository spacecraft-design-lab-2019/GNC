/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: iLQRv1.c
 *
 * MATLAB Coder version            : 4.3
 * C/C++ source code generated on  : 06-Feb-2020 15:05:27
 */

/* Include Files */
#include "iLQRv1.h"
#include "rt_nonfinite.h"
#include <math.h>
#include <string.h>

/* Variable Definitions */
static boolean_T isInitialized_iLQRv1 = false;

/* Function Declarations */
static void rkstep(const double x0[2], double u0, double dt, double x1[2],
                   double A[4], double B[2]);

/* Function Definitions */

/*
 * Explicit midpoint step from x_k to x_{k+1}
 * Arguments    : const double x0[2]
 *                double u0
 *                double dt
 *                double x1[2]
 *                double A[4]
 *                double B[2]
 * Return Type  : void
 */
static void rkstep(const double x0[2], double u0, double dt, double x1[2],
                   double A[4], double B[2])
{
  double dxdot1_idx_1;
  double a;
  double x_idx_0;
  double dxdot2[6];
  double A_tmp[4];
  signed char b_I[4];
  int i;

  /*  kg */
  /*  m */
  /*  kg m^2 /s */
  /*  m */
  /* m*l^2; % kg*m^2 */
  /*  m/s^2 */
  dxdot1_idx_1 = -4.905 * cos(x0[0]) / 0.25;
  a = 0.5 * dt;
  x_idx_0 = x0[0] + a * x0[1];
  a = x0[1] + a * (((u0 - 4.905 * sin(x0[0])) - 0.1 * x0[1]) / 0.25);

  /*  kg */
  /*  m */
  /*  kg m^2 /s */
  /*  m */
  /* m*l^2; % kg*m^2 */
  /*  m/s^2 */
  dxdot2[1] = -4.905 * cos(x_idx_0) / 0.25;
  dxdot2[3] = -0.4;
  x1[0] = x0[0] + dt * a;
  x1[1] = x0[1] + dt * (((u0 - 4.905 * sin(x_idx_0)) - 0.1 * a) / 0.25);
  a = 0.5 * dt * dt;
  dxdot2[0] = 0.0;
  dxdot2[4] = 0.0;
  A_tmp[0] = a * 0.0;
  A_tmp[1] = a * dxdot2[1];
  dxdot2[2] = 1.0;
  dxdot2[5] = 4.0;
  A_tmp[2] = a;
  A_tmp[3] = a * -0.4;
  b_I[1] = 0;
  b_I[2] = 0;
  b_I[0] = 1;
  b_I[3] = 1;
  for (i = 0; i < 2; i++) {
    a = A_tmp[i + 2];
    A[i] = ((double)b_I[i] + dt * dxdot2[i]) + (A_tmp[i] * 0.0 + a *
      dxdot1_idx_1);
    A[i + 2] = ((double)b_I[i + 2] + dt * dxdot2[i + 2]) + (A_tmp[i] + a * -0.4);
    B[i] = dt * dxdot2[i + 4] + (A_tmp[i] * 0.0 + a * 4.0);
  }
}

/*
 * iLQR Trajectory Optimization
 * Arguments    : const double x0[2]
 *                const double xg[2]
 *                const double utraj0_data[]
 *                const int utraj0_size[2]
 *                const double Q[4]
 *                double R
 *                const double Qf[4]
 *                double dt
 *                double tol
 *                double xtraj_data[]
 *                int xtraj_size[2]
 *                double utraj_data[]
 *                int utraj_size[2]
 *                double K_data[]
 *                int K_size[3]
 * Return Type  : void
 */
void iLQRv1(const double x0[2], const double xg[2], const double utraj0_data[],
            const int utraj0_size[2], const double Q[4], double R, const double
            Qf[4], double dt, double tol, double xtraj_data[], int xtraj_size[2],
            double utraj_data[], int utraj_size[2], double K_data[], int K_size
            [3])
{
  int N;
  int loop_ub_tmp;
  int loop_ub;
  double A_data[1600];
  int idx;
  double B_data[800];
  double J;
  int i;
  int k;
  double B_idx_0_tmp;
  double B_idx_0;
  int B_idx_1_tmp;
  double b_B_idx_1_tmp;
  double B_idx_1;
  double b_xtraj_data[2];
  double S_data[1604];
  double s_data[802];
  int l_size_idx_1;
  int unew_size_idx_1;
  double l_data[400];
  int exitg1;
  double unew_data[400];
  boolean_T exitg2;
  double d;
  int b_k;
  double r;
  int B_tmp_tmp;
  double xnew_data[802];
  double alpha;
  int i1;
  int i2;
  int i3;
  int i4;
  int i5;
  double b;
  int i6;
  double d1;
  double d2;
  double A[4];
  int K_data_tmp;
  double b_A[4];
  double c_A[4];
  if (isInitialized_iLQRv1 == false) {
    iLQRv1_initialize();
  }

  N = utraj0_size[1];
  xtraj_size[0] = 2;
  xtraj_size[1] = utraj0_size[1] + 1;
  loop_ub_tmp = (utraj0_size[1] + 1) << 1;
  if (0 <= loop_ub_tmp - 1) {
    memset(&xtraj_data[0], 0, loop_ub_tmp * sizeof(double));
  }

  xtraj_data[0] = x0[0];
  xtraj_data[1] = x0[1];
  loop_ub = utraj0_size[1] << 2;
  if (0 <= loop_ub - 1) {
    memset(&A_data[0], 0, loop_ub * sizeof(double));
  }

  idx = utraj0_size[1] << 1;
  if (0 <= idx - 1) {
    memset(&B_data[0], 0, idx * sizeof(double));
  }

  utraj_size[0] = 1;
  utraj_size[1] = utraj0_size[1];
  loop_ub = utraj0_size[0] * utraj0_size[1];
  if (0 <= loop_ub - 1) {
    memcpy(&utraj_data[0], &utraj0_data[0], loop_ub * sizeof(double));
  }

  /* First simulate with utraj0 to get initial matrices */
  J = 0.0;
  i = utraj0_size[1];
  for (k = 0; k < i; k++) {
    B_idx_0_tmp = xtraj_data[2 * k] - xg[0];
    B_idx_0 = 0.5 * B_idx_0_tmp;
    b_B_idx_1_tmp = xtraj_data[2 * k + 1] - xg[1];
    B_idx_1 = 0.5 * b_B_idx_1_tmp;
    J = (J + ((B_idx_0 * Q[0] + B_idx_1 * Q[1]) * B_idx_0_tmp + (B_idx_0 * Q[2]
           + B_idx_1 * Q[3]) * b_B_idx_1_tmp)) + 0.5 * utraj0_data[k] * R *
      utraj0_data[k];
    b_xtraj_data[0] = xtraj_data[2 * k];
    b_xtraj_data[1] = xtraj_data[1 + 2 * k];
    rkstep(b_xtraj_data, utraj0_data[k], dt, *(double (*)[2])&xtraj_data[2 * (k
            + 1)], *(double (*)[4])&A_data[4 * k], *(double (*)[2])&B_data[2 * k]);
  }

  B_idx_0_tmp = xtraj_data[2 * utraj0_size[1]] - xg[0];
  B_idx_0 = 0.5 * B_idx_0_tmp;
  B_idx_1_tmp = 2 * utraj0_size[1] + 1;
  b_B_idx_1_tmp = xtraj_data[B_idx_1_tmp] - xg[1];
  B_idx_1 = 0.5 * b_B_idx_1_tmp;
  J += (B_idx_0 * Qf[0] + B_idx_1 * Qf[1]) * B_idx_0_tmp + (B_idx_0 * Qf[2] +
    B_idx_1 * Qf[3]) * b_B_idx_1_tmp;

  /*  Jhist(1) = J; */
  loop_ub = (utraj0_size[1] + 1) << 2;
  if (0 <= loop_ub - 1) {
    memset(&S_data[0], 0, loop_ub * sizeof(double));
  }

  if (0 <= loop_ub_tmp - 1) {
    memset(&s_data[0], 0, loop_ub_tmp * sizeof(double));
  }

  K_size[0] = 1;
  K_size[1] = 2;
  K_size[2] = utraj0_size[1];
  if (0 <= idx - 1) {
    memset(&K_data[0], 0, idx * sizeof(double));
  }

  l_size_idx_1 = utraj0_size[1];
  loop_ub = utraj0_size[1];
  for (i = 0; i < loop_ub; i++) {
    l_data[i] = tol + 1.0;
  }

  unew_size_idx_1 = (short)l_size_idx_1;
  do {
    exitg1 = 0;
    for (k = 0; k < l_size_idx_1; k++) {
      unew_data[k] = fabs(l_data[k]);
    }

    if (unew_size_idx_1 <= 2) {
      if (unew_size_idx_1 == 1) {
        B_idx_0_tmp = unew_data[0];
      } else if ((unew_data[0] < unew_data[1]) || (rtIsNaN(unew_data[0]) &&
                  (!rtIsNaN(unew_data[1])))) {
        B_idx_0_tmp = unew_data[1];
      } else {
        B_idx_0_tmp = unew_data[0];
      }
    } else {
      if (!rtIsNaN(unew_data[0])) {
        idx = 1;
      } else {
        idx = 0;
        k = 2;
        exitg2 = false;
        while ((!exitg2) && (k <= unew_size_idx_1)) {
          if (!rtIsNaN(unew_data[k - 1])) {
            idx = k;
            exitg2 = true;
          } else {
            k++;
          }
        }
      }

      if (idx == 0) {
        B_idx_0_tmp = unew_data[0];
      } else {
        B_idx_0_tmp = unew_data[idx - 1];
        i = idx + 1;
        for (k = i; k <= unew_size_idx_1; k++) {
          d = unew_data[k - 1];
          if (B_idx_0_tmp < d) {
            B_idx_0_tmp = d;
          }
        }
      }
    }

    if (B_idx_0_tmp > tol) {
      /* Set up backwards LQR pass */
      S_data[4 * N] = Qf[0];
      S_data[4 * N + 1] = Qf[1];
      B_idx_1 = xtraj_data[2 * N] - xg[0];
      S_data[4 * N + 2] = Qf[2];
      S_data[4 * N + 3] = Qf[3];
      B_idx_0 = xtraj_data[B_idx_1_tmp] - xg[1];
      s_data[2 * N] = Qf[0] * B_idx_1 + Qf[2] * B_idx_0;
      s_data[B_idx_1_tmp] = Qf[1] * B_idx_1 + Qf[3] * B_idx_0;
      i = (int)(((-1.0 - ((double)(N + 1) - 1.0)) + 1.0) / -1.0);
      for (k = 0; k < i; k++) {
        b_k = (N - k) - 1;

        /* Calculate cost gradients for this time step */
        r = R * utraj_data[b_k];

        /* Calculate l and K */
        B_tmp_tmp = 2 * (b_k + 1);
        alpha = s_data[B_tmp_tmp + 1];
        i1 = 4 * (b_k + 1);
        d = B_data[2 * b_k];
        i2 = i1 + 1;
        i3 = 2 * b_k + 1;
        i4 = i1 + 2;
        i5 = i1 + 3;
        B_idx_0_tmp = d * S_data[i1] + B_data[i3] * S_data[i2];
        B_idx_0 = d * S_data[i4] + B_data[i3] * S_data[i5];
        b_B_idx_1_tmp = R + (B_idx_0_tmp * d + B_idx_0 * B_data[i3]);
        b = (r + (d * s_data[B_tmp_tmp] + B_data[i3] * alpha)) / b_B_idx_1_tmp;
        l_data[b_k] = b;

        /* Calculate new S and s         */
        for (i6 = 0; i6 < 2; i6++) {
          idx = 2 * i6 + 4 * b_k;
          loop_ub = i6 + 2 * b_k;
          K_data_tmp = idx + 1;
          K_data[loop_ub] = (B_idx_0_tmp * A_data[idx] + B_idx_0 *
                             A_data[K_data_tmp]) / b_B_idx_1_tmp;
          d1 = A_data[idx] - d * K_data[loop_ub];
          d2 = A_data[K_data_tmp] - B_data[i3] * K_data[loop_ub];
          b_A[i6] = d1 * S_data[i1] + d2 * S_data[i2];
          b_A[i6 + 2] = d1 * S_data[i4] + d2 * S_data[i5];
        }

        d1 = A_data[4 * b_k];
        d2 = K_data[2 * b_k];
        A[0] = d1 - d * d2;
        B_idx_0 = d2 * R;
        i6 = 4 * b_k + 2;
        A[2] = A_data[i6] - d * K_data[i3];
        idx = 4 * b_k + 1;
        A[1] = A_data[idx] - B_data[i3] * d2;
        b_B_idx_1_tmp = K_data[i3] * R;
        loop_ub = 4 * b_k + 3;
        A[3] = A_data[loop_ub] - B_data[i3] * K_data[i3];
        for (K_data_tmp = 0; K_data_tmp < 2; K_data_tmp++) {
          B_idx_0_tmp = b_A[K_data_tmp + 2];
          c_A[K_data_tmp] = b_A[K_data_tmp] * A[0] + B_idx_0_tmp * A[1];
          c_A[K_data_tmp + 2] = b_A[K_data_tmp] * A[2] + B_idx_0_tmp * A[3];
        }

        S_data[4 * b_k] = (Q[0] + B_idx_0 * d2) + c_A[0];
        S_data[idx] = (Q[1] + b_B_idx_1_tmp * d2) + c_A[1];
        B_idx_1 = xtraj_data[2 * b_k] - xg[0];
        S_data[i6] = (Q[2] + B_idx_0 * K_data[i3]) + c_A[2];
        S_data[loop_ub] = (Q[3] + b_B_idx_1_tmp * K_data[i3]) + c_A[3];
        B_idx_0 = xtraj_data[i3] - xg[1];
        b_B_idx_1_tmp = s_data[B_tmp_tmp] - (S_data[i1] * d + S_data[i4] *
          B_data[i3]) * b;
        B_idx_0_tmp = alpha - (S_data[i2] * d + S_data[i5] * B_data[i3]) * b;
        s_data[2 * b_k] = (((Q[0] * B_idx_1 + Q[2] * B_idx_0) - d2 * r) + d2 * R
                           * b) + ((d1 - d * d2) * b_B_idx_1_tmp + (A_data[idx]
          - B_data[i3] * d2) * B_idx_0_tmp);
        s_data[i3] = (((Q[1] * B_idx_1 + Q[3] * B_idx_0) - K_data[i3] * r) +
                      K_data[i3] * R * b) + ((A_data[i6] - d * K_data[i3]) *
          b_B_idx_1_tmp + (A_data[loop_ub] - B_data[i3] * K_data[i3]) *
          B_idx_0_tmp);
      }

      /* Now do forward pass line search with new l and K */
      if (0 <= N - 1) {
        memset(&unew_data[0], 0, N * sizeof(double));
      }

      idx = N + 1;
      if (0 <= loop_ub_tmp - 1) {
        memset(&xnew_data[0], 0, loop_ub_tmp * sizeof(double));
      }

      xnew_data[0] = xtraj_data[0];
      xnew_data[1] = xtraj_data[1];
      alpha = 1.0;
      r = J + 1.0;
      while (r > J) {
        r = 0.0;
        for (k = 0; k < N; k++) {
          i = 2 * k + 1;
          d = (utraj_data[k] - alpha * l_data[k]) - (K_data[2 * k] * (xnew_data
            [2 * k] - xtraj_data[2 * k]) + K_data[i] * (xnew_data[i] -
            xtraj_data[i]));
          unew_data[k] = d;
          b_xtraj_data[0] = xnew_data[2 * k];
          b_xtraj_data[1] = xnew_data[1 + 2 * k];
          rkstep(b_xtraj_data, d, dt, *(double (*)[2])&xnew_data[2 * (k + 1)],
                 *(double (*)[4])&A_data[4 * k], *(double (*)[2])&B_data[2 * k]);
          B_idx_0_tmp = xnew_data[2 * k] - xg[0];
          B_idx_0 = 0.5 * B_idx_0_tmp;
          b_B_idx_1_tmp = xnew_data[i] - xg[1];
          B_idx_1 = 0.5 * b_B_idx_1_tmp;
          r = (r + ((B_idx_0 * Q[0] + B_idx_1 * Q[1]) * B_idx_0_tmp + (B_idx_0 *
                 Q[2] + B_idx_1 * Q[3]) * b_B_idx_1_tmp)) + 0.5 * d * R * d;
        }

        B_idx_0_tmp = xnew_data[2 * N] - xg[0];
        B_idx_0 = 0.5 * B_idx_0_tmp;
        b_B_idx_1_tmp = xnew_data[B_idx_1_tmp] - xg[1];
        B_idx_1 = 0.5 * b_B_idx_1_tmp;
        r += (B_idx_0 * Qf[0] + B_idx_1 * Qf[1]) * B_idx_0_tmp + (B_idx_0 * Qf[2]
          + B_idx_1 * Qf[3]) * b_B_idx_1_tmp;
        alpha *= 0.5;
      }

      xtraj_size[0] = 2;
      xtraj_size[1] = idx;
      loop_ub = 2 * idx;
      if (0 <= loop_ub - 1) {
        memcpy(&xtraj_data[0], &xnew_data[0], loop_ub * sizeof(double));
      }

      utraj_size[0] = 1;
      utraj_size[1] = N;
      if (0 <= N - 1) {
        memcpy(&utraj_data[0], &unew_data[0], N * sizeof(double));
      }

      J = r;

      /*      Jhist(iter+1) = J; */
    } else {
      exitg1 = 1;
    }
  } while (exitg1 == 0);

  /*   */
  /*  disp('Iterations: '); */
  /*  disp(iter); */
  /*   */
  /*  disp('Final Cost: '); */
  /*  disp(Jnew); */
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void iLQRv1_initialize(void)
{
  rt_InitInfAndNaN();
  isInitialized_iLQRv1 = true;
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void iLQRv1_terminate(void)
{
  /* (no terminate code required) */
  isInitialized_iLQRv1 = false;
}

/*
 * File trailer for iLQRv1.c
 *
 * [EOF]
 */
