/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: iLQRsimple.c
 *
 * MATLAB Coder version            : 4.3
 * C/C++ source code generated on  : 10-Feb-2020 20:32:25
 */

/* Include Files */
#include "iLQRsimple.h"
#include "rt_nonfinite.h"
#include <math.h>
#include <string.h>

/* Variable Definitions */
static boolean_T isInitialized_iLQRsimple = false;

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


  dxdot1_idx_1 = -4.905 * cos(x0[0]) / 0.25;
  a = 0.5 * dt;
  x_idx_0 = x0[0] + a * x0[1];
  a = x0[1] + a * (((u0 - 4.905 * sin(x0[0])) - 0.1 * x0[1]) / 0.25);


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
 *                const double utraj0[399]
 *                const double Q[4]
 *                double R
 *                const double Qf[4]
 *                double dt
 *                double tol
 *                double xtraj[800]
 *                double utraj[399]
 *                double K[798]
 * Return Type  : void
 */
void iLQRsimple(const double x0[2], const double xg[2], const double utraj0[399],
                const double Q[4], double R, const double Qf[4], double dt,
                double tol, double xtraj[800], double utraj[399], double K[798])
{
  double J;
  int k;
  double snew[2];
  int idx;
  double b_xtraj[2];
  double A[1596];
  double B[798];
  double l[399];
  int exitg1;
  double unew[399];
  boolean_T exitg2;
  double alpha;
  double d;
  double S[4];
  double Jnew;
  double xnew[800];
  double r;
  int B_tmp;
  int b_B_tmp;
  double d1;
  double l_tmp;
  double b_l;
  int i;
  double d2;
  int i1;
  int i2;
  double a_tmp;
  double a[4];
  double b_a_tmp;
  double A_tmp;
  double c_a_tmp;
  double b_Q[2];
  double b_a[4];
  double b_A[4];
  if (isInitialized_iLQRsimple == false) {
    iLQRsimple_initialize();
  }

  memset(&xtraj[0], 0, 800U * sizeof(double));
  xtraj[0] = x0[0];
  xtraj[1] = x0[1];

  /* First simulate with utraj0 to get initial matrices */
  J = 0.0;
  for (k = 0; k < 399; k++) {
    utraj[k] = utraj0[k];
    idx = k << 1;
    snew[0] = xtraj[idx] - xg[0];
    snew[1] = xtraj[idx + 1] - xg[1];
    J = (J + ((0.5 * snew[0] * Q[0] + 0.5 * snew[1] * Q[1]) * snew[0] + (0.5 *
           snew[0] * Q[2] + 0.5 * snew[1] * Q[3]) * snew[1])) + 0.5 * utraj0[k] *
      R * utraj0[k];
    b_xtraj[0] = xtraj[k << 1];
    b_xtraj[1] = xtraj[1 + (k << 1)];
    rkstep(b_xtraj, utraj0[k], dt, *(double (*)[2])&xtraj[(k + 1) << 1],
           *(double (*)[4])&A[k << 2], *(double (*)[2])&B[k << 1]);
  }

  snew[0] = xtraj[798] - xg[0];
  snew[1] = xtraj[799] - xg[1];
  J += (0.5 * snew[0] * Qf[0] + 0.5 * snew[1] * Qf[1]) * snew[0] + (0.5 * snew[0]
    * Qf[2] + 0.5 * snew[1] * Qf[3]) * snew[1];

  /*  Jhist(1) = J; */
  /*  Set up backwards pass matrices */
  memset(&K[0], 0, 798U * sizeof(double));
  for (idx = 0; idx < 399; idx++) {
    l[idx] = tol + 1.0;
  }

  /*  temp matrices */
  /*  Line search temp matrices */
  do {
    exitg1 = 0;
    for (k = 0; k < 399; k++) {
      unew[k] = fabs(l[k]);
    }

    if (!rtIsNaN(unew[0])) {
      idx = 1;
    } else {
      idx = 0;
      k = 2;
      exitg2 = false;
      while ((!exitg2) && (k < 400)) {
        if (!rtIsNaN(unew[k - 1])) {
          idx = k;
          exitg2 = true;
        } else {
          k++;
        }
      }
    }

    if (idx == 0) {
      alpha = unew[0];
    } else {
      alpha = unew[idx - 1];
      idx++;
      for (k = idx; k < 400; k++) {
        d = unew[k - 1];
        if (alpha < d) {
          alpha = d;
        }
      }
    }

    if (alpha > tol) {
      /* Set up backwards LQR pass */
      S[0] = Qf[0];
      S[1] = Qf[1];
      S[2] = Qf[2];
      S[3] = Qf[3];
      alpha = xtraj[798] - xg[0];
      Jnew = xtraj[799] - xg[1];
      snew[0] = Qf[0] * alpha + Qf[2] * Jnew;
      snew[1] = Qf[1] * alpha + Qf[3] * Jnew;

      /* Now do forward pass line search with new l and K */
      for (k = 0; k < 399; k++) {
        /* Calculate cost gradients for this time step */
        r = R * utraj[398 - k];

        /* Calculate l and K */
        B_tmp = (398 - k) << 1;
        b_B_tmp = B_tmp + 1;
        d = B[B_tmp] * S[0];
        d1 = B[b_B_tmp] * S[3];
        alpha = d + B[b_B_tmp] * S[1];
        Jnew = B[B_tmp] * S[2] + d1;
        l_tmp = R + (alpha * B[B_tmp] + Jnew * B[b_B_tmp]);
        b_l = (r + (B[B_tmp] * snew[0] + B[b_B_tmp] * snew[1])) / l_tmp;
        l[398 - k] = b_l;
        idx = (398 - k) << 2;
        i = idx + 1;
        d2 = (alpha * A[idx] + Jnew * A[i]) / l_tmp;
        K[B_tmp] = d2;
        i1 = idx + 2;
        i2 = idx + 3;
        l_tmp = (alpha * A[i1] + Jnew * A[i2]) / l_tmp;
        K[b_B_tmp] = l_tmp;

        /* Calculate new S and s         */
        alpha = xtraj[B_tmp] - xg[0];
        a_tmp = A[idx] - B[B_tmp] * d2;
        a[0] = a_tmp;
        b_a_tmp = A[i1] - B[B_tmp] * l_tmp;
        a[2] = b_a_tmp;
        A_tmp = A[i] - B[b_B_tmp] * d2;
        d = snew[0] - (d + S[2] * B[b_B_tmp]) * b_l;
        Jnew = xtraj[b_B_tmp] - xg[1];
        a[1] = A_tmp;
        c_a_tmp = A[i2] - B[b_B_tmp] * l_tmp;
        a[3] = c_a_tmp;
        d1 = snew[1] - (S[1] * B[B_tmp] + d1) * b_l;
        b_Q[0] = ((Q[0] * alpha + Q[2] * Jnew) - d2 * r) + d2 * R * b_l;
        b_xtraj[0] = a_tmp * d + A_tmp * d1;
        b_Q[1] = ((Q[1] * alpha + Q[3] * Jnew) - l_tmp * r) + l_tmp * R * b_l;
        b_xtraj[1] = b_a_tmp * d + c_a_tmp * d1;
        for (idx = 0; idx < 2; idx++) {
          snew[idx] = b_Q[idx] + b_xtraj[idx];
          i = idx << 1;
          i1 = i + 1;
          b_A[idx] = a[i] * S[0] + a[i1] * S[1];
          b_A[idx + 2] = a[i] * S[2] + a[i1] * S[3];
        }

        for (idx = 0; idx < 2; idx++) {
          alpha = K[idx + B_tmp] * R;
          S[idx] = Q[idx] + alpha * K[B_tmp];
          d = b_A[idx + 2];
          b_a[idx] = b_A[idx] * a_tmp + d * A_tmp;
          S[idx + 2] = Q[idx + 2] + alpha * K[b_B_tmp];
          b_a[idx + 2] = b_A[idx] * b_a_tmp + d * c_a_tmp;
        }

        S[0] += b_a[0];
        S[1] += b_a[1];
        S[2] += b_a[2];
        S[3] += b_a[3];
        unew[k] = 0.0;
      }

      memset(&xnew[0], 0, 800U * sizeof(double));
      xnew[0] = xtraj[0];
      xnew[1] = xtraj[1];
      alpha = 1.0;
      Jnew = J + 1.0;
      while (Jnew > J) {
        Jnew = 0.0;
        for (k = 0; k < 399; k++) {
          idx = k << 1;
          B_tmp = idx + 1;
          d = (utraj[k] - alpha * l[k]) - (K[idx] * (xnew[idx] - xtraj[idx]) +
            K[B_tmp] * (xnew[B_tmp] - xtraj[B_tmp]));
          unew[k] = d;
          b_xtraj[0] = xnew[k << 1];
          b_xtraj[1] = xnew[1 + (k << 1)];
          rkstep(b_xtraj, d, dt, *(double (*)[2])&xnew[(k + 1) << 1], *(double (*)
                  [4])&A[k << 2], *(double (*)[2])&B[k << 1]);
          snew[0] = xnew[idx] - xg[0];
          snew[1] = xnew[B_tmp] - xg[1];
          Jnew = (Jnew + ((0.5 * snew[0] * Q[0] + 0.5 * snew[1] * Q[1]) * snew[0]
                          + (0.5 * snew[0] * Q[2] + 0.5 * snew[1] * Q[3]) *
                          snew[1])) + 0.5 * d * R * d;
        }

        snew[0] = xnew[798] - xg[0];
        snew[1] = xnew[799] - xg[1];
        Jnew += (0.5 * snew[0] * Qf[0] + 0.5 * snew[1] * Qf[1]) * snew[0] + (0.5
          * snew[0] * Qf[2] + 0.5 * snew[1] * Qf[3]) * snew[1];
        alpha *= 0.5;
      }

      memcpy(&xtraj[0], &xnew[0], 800U * sizeof(double));
      memcpy(&utraj[0], &unew[0], 399U * sizeof(double));
      J = Jnew;

      /*      Jhist(iter+1) = J; */
    } else {
      exitg1 = 1;
    }
  } while (exitg1 == 0);
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void iLQRsimple_initialize(void)
{
  rt_InitInfAndNaN();
  isInitialized_iLQRsimple = true;
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void iLQRsimple_terminate(void)
{
  /* (no terminate code required) */
  isInitialized_iLQRsimple = false;
}

/*
 * File trailer for iLQRsimple.c
 *
 * [EOF]
 */
