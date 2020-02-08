/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: iLQRv1.c
 *
 * MATLAB Coder version            : 4.3
 * C/C++ source code generated on  : 07-Feb-2020 14:37:48
 */

/* Include Files */
#include "iLQRv1.h"
#include "rt_nonfinite.h"
#include <math.h>
#include <string.h>

/* Variable Definitions */
static boolean_T isInitialized_iLQRv1 = false;

/* Function Declarations */
static void rkstep(const real_T x0[2], real_T u0, real_T dt, real_T x1[2],
                   real_T A[4], real_T B[2]);

/* Function Definitions */

/*
 * Explicit midpoint step from x_k to x_{k+1}
 * Arguments    : const real_T x0[2]
 *                real_T u0
 *                real_T dt
 *                real_T x1[2]
 *                real_T A[4]
 *                real_T B[2]
 * Return Type  : void
 */
static void rkstep(const real_T x0[2], real_T u0, real_T dt, real_T x1[2],
                   real_T A[4], real_T B[2])
{
  real_T dxdot1_idx_1;
  real_T a;
  real_T x_idx_0;
  real_T dxdot2[6];
  real_T A_tmp[4];
  int8_T b_I[4];
  int32_T i;

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
    A[i] = ((real_T)b_I[i] + dt * dxdot2[i]) + (A_tmp[i] * 0.0 + a *
      dxdot1_idx_1);
    A[i + 2] = ((real_T)b_I[i + 2] + dt * dxdot2[i + 2]) + (A_tmp[i] + a * -0.4);
    B[i] = dt * dxdot2[i + 4] + (A_tmp[i] * 0.0 + a * 4.0);
  }
}

/*
 * iLQR Trajectory Optimization
 * Arguments    : const real_T x0[2]
 *                const real_T xg[2]
 *                const real_T utraj0[399]
 *                const real_T Q[4]
 *                real_T R
 *                const real_T Qf[4]
 *                real_T dt
 *                real_T tol
 *                real_T xtraj[800]
 *                real_T utraj[399]
 *                real_T K[798]
 * Return Type  : void
 */
void iLQRv1(const real_T x0[2], const real_T xg[2], const real_T utraj0[399],
            const real_T Q[4], real_T R, const real_T Qf[4], real_T dt, real_T
            tol, real_T xtraj[800], real_T utraj[399], real_T K[798])
{
  real_T J;
  int32_T k;
  real_T b_idx_0;
  real_T b_idx_1;
  int32_T idx;
  real_T S[1600];
  real_T s[800];
  real_T b_xtraj[2];
  real_T A[1596];
  real_T B[798];
  real_T l[399];
  int32_T exitg1;
  real_T unew[399];
  boolean_T exitg2;
  real_T alpha;
  real_T d;
  real_T xnew[800];
  real_T r;
  int32_T B_tmp;
  real_T Jnew;
  int32_T b_B_tmp;
  int32_T c_B_tmp;
  int32_T d_B_tmp;
  int32_T i;
  int32_T i1;
  int32_T i2;
  real_T b_l;
  int32_T i3;
  int32_T i4;
  int32_T i5;
  int32_T i6;
  real_T a[4];
  int32_T i7;
  real_T b_Q[4];
  real_T b_a[4];
  int32_T i8;
  int32_T i9;
  if (!isInitialized_iLQRv1) {
    iLQRv1_initialize();
  }

  memset(&xtraj[0], 0, 800U * sizeof(real_T));
  xtraj[0] = x0[0];
  xtraj[1] = x0[1];

  /* First simulate with utraj0 to get initial matrices */
  J = 0.0;
  for (k = 0; k < 399; k++) {
    utraj[k] = utraj0[k];
    idx = k << 1;
    b_idx_0 = xtraj[idx] - xg[0];
    b_idx_1 = xtraj[idx + 1] - xg[1];
    J = (J + ((0.5 * b_idx_0 * Q[0] + 0.5 * b_idx_1 * Q[1]) * b_idx_0 + (0.5 *
           b_idx_0 * Q[2] + 0.5 * b_idx_1 * Q[3]) * b_idx_1)) + 0.5 * utraj0[k] *
      R * utraj0[k];
    b_xtraj[0] = xtraj[k << 1];
    b_xtraj[1] = xtraj[1 + (k << 1)];
    rkstep(b_xtraj, utraj0[k], dt, *(real_T (*)[2])&xtraj[(k + 1) << 1],
           *(real_T (*)[4])&A[k << 2], *(real_T (*)[2])&B[k << 1]);
  }

  b_idx_0 = xtraj[798] - xg[0];
  b_idx_1 = xtraj[799] - xg[1];
  J += (0.5 * b_idx_0 * Qf[0] + 0.5 * b_idx_1 * Qf[1]) * b_idx_0 + (0.5 *
    b_idx_0 * Qf[2] + 0.5 * b_idx_1 * Qf[3]) * b_idx_1;

  /*  Jhist(1) = J; */
  memset(&S[0], 0, 1600U * sizeof(real_T));
  memset(&s[0], 0, 800U * sizeof(real_T));
  memset(&K[0], 0, 798U * sizeof(real_T));
  for (idx = 0; idx < 399; idx++) {
    l[idx] = tol + 1.0;
  }

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
      S[1596] = Qf[0];
      S[1597] = Qf[1];
      b_idx_0 = xtraj[798] - xg[0];
      S[1598] = Qf[2];
      S[1599] = Qf[3];
      b_idx_1 = xtraj[799] - xg[1];
      s[798] = Qf[0] * b_idx_0 + Qf[2] * b_idx_1;
      s[799] = Qf[1] * b_idx_0 + Qf[3] * b_idx_1;

      /* Now do forward pass line search with new l and K */
      for (k = 0; k < 399; k++) {
        /* Calculate cost gradients for this time step */
        r = R * utraj[398 - k];

        /* Calculate l and K */
        B_tmp = (398 - k) << 1;
        b_B_tmp = (399 - k) << 1;
        c_B_tmp = B_tmp + 1;
        d_B_tmp = b_B_tmp + 1;
        idx = (399 - k) << 2;
        i = idx + 1;
        i1 = idx + 2;
        i2 = idx + 3;
        alpha = B[B_tmp] * S[idx] + B[c_B_tmp] * S[i];
        Jnew = B[B_tmp] * S[i1] + B[c_B_tmp] * S[i2];
        b_idx_0 = R + (alpha * B[B_tmp] + Jnew * B[c_B_tmp]);
        b_l = (r + (B[B_tmp] * s[b_B_tmp] + B[c_B_tmp] * s[d_B_tmp])) / b_idx_0;
        l[398 - k] = b_l;
        i3 = (398 - k) << 2;
        i4 = i3 + 1;
        d = (alpha * A[i3] + Jnew * A[i4]) / b_idx_0;
        K[B_tmp] = d;
        i5 = i3 + 2;
        i6 = i3 + 3;
        alpha = (alpha * A[i5] + Jnew * A[i6]) / b_idx_0;
        K[c_B_tmp] = alpha;

        /* Calculate new S and s         */
        a[0] = A[i3] - B[B_tmp] * d;
        a[2] = A[i5] - B[B_tmp] * alpha;
        a[1] = A[i4] - B[c_B_tmp] * d;
        a[3] = A[i6] - B[c_B_tmp] * alpha;
        for (i7 = 0; i7 < 2; i7++) {
          i8 = i7 << 1;
          i9 = i8 + 1;
          d = a[i8] * S[idx] + a[i9] * S[i];
          Jnew = K[i7 + B_tmp] * R;
          b_Q[i7] = Q[i7] + Jnew * K[B_tmp];
          alpha = a[i8] * S[i1] + a[i9] * S[i2];
          b_Q[i7 + 2] = Q[i7 + 2] + Jnew * K[c_B_tmp];
          b_a[i7] = d * a[0] + alpha * a[1];
          b_a[i7 + 2] = d * a[2] + alpha * a[3];
        }

        S[i3] = b_Q[0] + b_a[0];
        S[i4] = b_Q[1] + b_a[1];
        b_idx_0 = xtraj[B_tmp] - xg[0];
        S[i5] = b_Q[2] + b_a[2];
        S[i6] = b_Q[3] + b_a[3];
        b_idx_1 = xtraj[c_B_tmp] - xg[1];
        alpha = s[b_B_tmp] - (S[idx] * B[B_tmp] + S[i1] * B[c_B_tmp]) * b_l;
        Jnew = s[d_B_tmp] - (S[i] * B[B_tmp] + S[i2] * B[c_B_tmp]) * b_l;
        s[B_tmp] = (((Q[0] * b_idx_0 + Q[2] * b_idx_1) - K[B_tmp] * r) + K[B_tmp]
                    * R * b_l) + ((A[i3] - B[B_tmp] * K[B_tmp]) * alpha + (A[i4]
          - B[c_B_tmp] * K[B_tmp]) * Jnew);
        s[c_B_tmp] = (((Q[1] * b_idx_0 + Q[3] * b_idx_1) - K[c_B_tmp] * r) +
                      K[c_B_tmp] * R * b_l) + ((A[i5] - B[B_tmp] * K[c_B_tmp]) *
          alpha + (A[i6] - B[c_B_tmp] * K[c_B_tmp]) * Jnew);
        unew[k] = 0.0;
      }

      memset(&xnew[0], 0, 800U * sizeof(real_T));
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
          rkstep(b_xtraj, d, dt, *(real_T (*)[2])&xnew[(k + 1) << 1], *(real_T (*)
                  [4])&A[k << 2], *(real_T (*)[2])&B[k << 1]);
          b_idx_0 = xnew[idx] - xg[0];
          b_idx_1 = xnew[B_tmp] - xg[1];
          Jnew = (Jnew + ((0.5 * b_idx_0 * Q[0] + 0.5 * b_idx_1 * Q[1]) *
                          b_idx_0 + (0.5 * b_idx_0 * Q[2] + 0.5 * b_idx_1 * Q[3])
                          * b_idx_1)) + 0.5 * d * R * d;
        }

        b_idx_0 = xnew[798] - xg[0];
        b_idx_1 = xnew[799] - xg[1];
        Jnew += (0.5 * b_idx_0 * Qf[0] + 0.5 * b_idx_1 * Qf[1]) * b_idx_0 + (0.5
          * b_idx_0 * Qf[2] + 0.5 * b_idx_1 * Qf[3]) * b_idx_1;
        alpha *= 0.5;
      }

      memcpy(&xtraj[0], &xnew[0], 800U * sizeof(real_T));
      memcpy(&utraj[0], &unew[0], 399U * sizeof(real_T));
      J = Jnew;

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
