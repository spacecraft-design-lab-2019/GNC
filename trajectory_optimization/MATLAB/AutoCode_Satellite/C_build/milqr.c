/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: milqr.c
 *
 * MATLAB Coder version            : 4.3
 * C/C++ source code generated on  : 12-Apr-2020 23:15:53
 */

/* Include Files */
#include "milqr.h"
#include <math.h>
#include <string.h>

/* Variable Definitions */
static const signed char iv[3][3] = { { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };

static const signed char iv1[6][3] = { { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, {
    1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };

static const signed char iv2[7][7] = { { 1, 0, 0, 0, 0, 0, 0 }, { 0, 1, 0, 0, 0,
    0, 0 }, { 0, 0, 1, 0, 0, 0, 0 }, { 0, 0, 0, 1, 0, 0, 0 }, { 0, 0, 0, 0, 1, 0,
    0 }, { 0, 0, 0, 0, 0, 1, 0 }, { 0, 0, 0, 0, 0, 0, 1 } };

static const float fv[6][3] = { { 0.0F, 0.0F, 0.0F }, { 0.0F, 0.0F, 0.0F }, {
    0.0F, 0.0F, 0.0F }, { 0.01F, 0.0F, 0.0F }, { 0.0F, 0.01F, 0.0F }, { 0.0F,
    0.0F, 0.01F } };

static const float fv1[3][3] = { { 0.01F, 0.0F, 0.0F }, { 0.0F, 0.01F, 0.0F }, {
    0.0F, 0.0F, 0.01F } };

static const float fv2[3][3] = { { 0.05F, 0.0F, 0.0F }, { 0.0F, 0.05F, 0.0F }, {
    0.0F, 0.0F, 0.05F } };

/* Function Declarations */
static void b_abs(float x[4999][3], float y[4999][3]);
static void b_forwardRollout(float x[5000][7], const float xg[7], float u[4999]
  [3], float l[4999][3], float K[4999][6][3], float alpha, float u_lims[2][3],
  float dt, float B_ECI[5000][3], float Js[6][3], float xnew[5000][7], float
  unew[4999][3], float fx[4999][6][6], float fu[4999][3][6], float cx[5000][6],
  float cu[4999][3], float cxx[5000][6][6], float cuu[4999][3][3], float *cost);
static void backwardPass(float fx[4999][6][6], float fu[4999][3][6], float cx
  [5000][6], float cu[4999][3], float cxx[5000][6][6], float cuu[4999][3][3],
  float lambda, float u_lims[2][3], float u[4999][3], float l[4999][3], float K
  [4999][6][3], float dV[2], bool *diverged);
static void boxQPsolve(float Quu[3][3], const float Qu[3], const float
  lower_lim[3], const float upper_lim[3], const float u0[3], float u[3], float
  *result, float Luu_data[], int Luu_size[2], bool b_free[3]);
static void chol_free(const float A_data[], const int A_size[2], float L_data[],
                      int L_size[2], float *fail);
static void forwardRollout(float x[5000][7], const float xg[7], float u[4999][3],
  float K[4999][6][3], float u_lims[2][3], float dt, float B_ECI[5000][3], float
  Js[6][3], float xnew[5000][7], float unew[4999][3], float fx[4999][6][6],
  float fu[4999][3][6], float cx[5000][6], float cu[4999][3], float cxx[5000][6]
  [6], float cuu[4999][3][3], float *cost);
static void power(float y[11]);
static void satellite_dynamics(const float x[7], const float u[3], const float
  B_ECI[3], float Js[6][3], float xdot[7], float dxdot[10][7]);

/* Function Definitions */

/*
 * Arguments    : float x[4999][3]
 *                float y[4999][3]
 * Return Type  : void
 */
static void b_abs(float x[4999][3], float y[4999][3])
{
  int k;
  for (k = 0; k < 4999; k++) {
    y[k][0] = fabsf(x[k][0]);
    y[k][1] = fabsf(x[k][1]);
    y[k][2] = fabsf(x[k][2]);
  }
}

/*
 * Arguments    : float x[5000][7]
 *                const float xg[7]
 *                float u[4999][3]
 *                float l[4999][3]
 *                float K[4999][6][3]
 *                float alpha
 *                float u_lims[2][3]
 *                float dt
 *                float B_ECI[5000][3]
 *                float Js[6][3]
 *                float xnew[5000][7]
 *                float unew[4999][3]
 *                float fx[4999][6][6]
 *                float fu[4999][3][6]
 *                float cx[5000][6]
 *                float cu[4999][3]
 *                float cxx[5000][6][6]
 *                float cuu[4999][3][3]
 *                float *cost
 * Return Type  : void
 */
static void b_forwardRollout(float x[5000][7], const float xg[7], float u[4999]
  [3], float l[4999][3], float K[4999][6][3], float alpha, float u_lims[2][3],
  float dt, float B_ECI[5000][3], float Js[6][3], float xnew[5000][7], float
  unew[4999][3], float fx[4999][6][6], float fu[4999][3][6], float cx[5000][6],
  float cu[4999][3], float cxx[5000][6][6], float cuu[4999][3][3], float *cost)
{
  int i;
  int i1;
  float a;
  float b_a;
  int k;
  float dx[6];
  float q[4][4];
  float xg_tmp;
  float out;
  int b_sign;
  float cost_tmp[3];
  int c_a;
  int i2;
  int i3;
  float q_error;
  int d_a;
  int i4;
  float b_fv[3][3];
  float y;
  float f;
  float f1;
  float b_q_error[4];
  float f2;
  float f3;
  float f4;
  int i5;
  float e_a[4][3];
  float f5;
  int i6;
  int b_k;
  float f6;
  int i7;
  float b[7];
  float dxdot1[10][7];
  int i8;
  float b_xnew[7];
  float dxdot2[10][7];
  int i9;
  float x1[7];
  float b_y;
  int i10;
  int i11;
  float fx_tmp[7][7];
  int i12;
  int i13;
  float b_fx_tmp[7][6];
  int i14;
  float f7;
  int i15;
  float b_fv1[7][7];
  int i16;
  int i17;
  int i18;
  int i19;
  int i20;
  float c_xnew[6][7];
  float f8;
  int i21;
  int i22;
  float c_fx_tmp[7][6];
  int i23;
  int i24;
  int i25;
  int i26;
  float f9;
  int i27;
  int i28;
  float f10;
  int i29;
  float f11;
  int i30;
  float b_dt[3][7];
  float b_xg_tmp;
  float b_out;
  int c_sign;
  int i31;
  int i32;
  int i33;
  float d_sign[4][3];
  float f12;
  float f13;
  float f14;
  float f15;
  int i34;
  float b_fv2[3];

  /*  Uses an rk method to roll out a trajectory */
  /*  Returns the new trajectory, cost and the derivatives along the trajectory */
  /*  If feed-forward and feed-back controls l, K are non-zero */
  /*  then the unew returned is the new control sequeunce. Otherwise unew = u */
  /*  Sizes */
  /*  Initialize outputs (Don't overwrite x or u) */
  memset(&xnew[0][0], 0, 35000U * sizeof(float));
  for (i = 0; i < 4999; i++) {
    unew[i][0] = 0.0F;
    unew[i][1] = 0.0F;
    unew[i][2] = 0.0F;
  }

  *cost = 0.0F;
  for (i1 = 0; i1 < 7; i1++) {
    xnew[0][i1] = x[0][i1];
  }

  a = 0.5F * dt;
  b_a = 0.5F * dt * dt;
  for (k = 0; k < 4999; k++) {
    /*  Find the state error vector dx */
    dx[3] = xnew[k][4] - x[k][4];
    dx[4] = xnew[k][5] - x[k][5];
    dx[5] = xnew[k][6] - x[k][6];

    /*  Calculate error between qk and q_nom */
    /*  Defined as conj(q_nom)*qnew */
    /*  Returns error as Rodrigues parameters (3x1) */
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
    q[0][0] = x[k][0];
    q[1][0] = x[k][1];
    q[2][0] = x[k][2];
    q[3][0] = x[k][3];
    q[0][1] = -x[k][1];
    q[1][1] = x[k][0];
    q[2][1] = x[k][3];
    q[3][1] = -x[k][2];
    q[0][2] = -x[k][2];
    q[1][2] = -x[k][3];
    q[2][2] = x[k][0];
    q[3][2] = x[k][1];
    q[0][3] = -x[k][3];
    q[1][3] = x[k][2];
    q[2][3] = -x[k][1];
    q[3][3] = x[k][0];
    q_error = 0.0F;
    for (i4 = 0; i4 < 4; i4++) {
      f = ((q[0][i4] * xnew[k][0] + q[1][i4] * xnew[k][1]) + q[2][i4] * xnew[k]
           [2]) + q[3][i4] * xnew[k][3];
      b_q_error[i4] = f;
      q_error += f * f;
    }

    y = sqrtf(q_error);
    f1 = b_q_error[0] / y;
    b_q_error[0] = f1;
    f2 = b_q_error[1] / y;
    b_q_error[1] = f2;
    f3 = b_q_error[2] / y;
    b_q_error[2] = f3;
    f4 = b_q_error[3] / y;
    b_q_error[3] = f4;

    /*  re-normalize */
    /*  inverse Cayley Map */
    dx[0] = f2 / f1;
    dx[1] = f3 / f1;
    dx[2] = f4 / f1;

    /*  Find the new control and ensure it is within the limits */
    for (b_k = 0; b_k < 3; b_k++) {
      f6 = 0.0F;
      for (i7 = 0; i7 < 6; i7++) {
        f6 += K[k][i7][b_k] * dx[i7];
      }

      unew[k][b_k] = fminf(u_lims[1][b_k], fmaxf(u_lims[0][b_k], (u[k][b_k] -
        alpha * l[k][b_k]) - f6));
    }

    /*  Step the dynamics forward */
    /*  Steps the dynamics forward using a 2nd order rk-method */
    /*  Returns: new state, discrete time Jacobians */
    /*  Step dynamics */
    /*  Explicit midpoint step from x_k to x_{k+1} */
    satellite_dynamics(*(float (*)[7])&xnew[k][0], *(float (*)[3])&unew[k][0],
                       *(float (*)[3])&B_ECI[k][0], Js, b, dxdot1);
    for (i8 = 0; i8 < 7; i8++) {
      b_xnew[i8] = xnew[k][i8] + a * b[i8];
    }

    satellite_dynamics(b_xnew, *(float (*)[3])&unew[k][0], *(float (*)[3])&
                       B_ECI[k][0], Js, b, dxdot2);
    for (i9 = 0; i9 < 7; i9++) {
      x1[i9] = xnew[k][i9] + dt * b[i9];
    }

    /*  Re-normalize the quaternion */
    b_y = sqrtf(((x1[0] * x1[0] + x1[1] * x1[1]) + x1[2] * x1[2]) + x1[3] * x1[3]);

    /*  Continuous time Jacobians */
    /*  Discrete time Jacobians */
    /*  First form the state error Jacobians E(x_k) & E(x_k+1) */
    /*  Returns skew symmetric - cross product matrix of a 3x1 vector */
    /*  Returns skew symmetric - cross product matrix of a 3x1 vector */
    for (i10 = 0; i10 < 7; i10++) {
      for (i11 = 0; i11 < 7; i11++) {
        fx_tmp[i10][i11] = b_a * dxdot2[i10][i11];
      }
    }

    b_fv[0][0] = 0.0F;
    b_fv[1][0] = -x1[3];
    b_fv[2][0] = x1[2];
    b_fv[0][1] = x1[3];
    b_fv[1][1] = 0.0F;
    b_fv[2][1] = -x1[1];
    b_fv[0][2] = -x1[2];
    b_fv[1][2] = x1[1];
    b_fv[2][2] = 0.0F;
    for (i12 = 0; i12 < 3; i12++) {
      b_fx_tmp[0][i12] = -x1[i12 + 1];
      b_fx_tmp[0][i12 + 3] = 0.0F;
      b_fx_tmp[i12 + 1][0] = x1[0] * (float)iv[0][i12] + b_fv[0][i12];
      b_fx_tmp[i12 + 1][3] = 0.0F;
      b_fx_tmp[i12 + 1][1] = x1[0] * (float)iv[1][i12] + b_fv[1][i12];
      b_fx_tmp[i12 + 1][4] = 0.0F;
      b_fx_tmp[i12 + 1][2] = x1[0] * (float)iv[2][i12] + b_fv[2][i12];
      b_fx_tmp[i12 + 1][5] = 0.0F;
      for (i16 = 0; i16 < 6; i16++) {
        b_fx_tmp[i12 + 4][i16] = iv1[i16][i12];
      }
    }

    for (i13 = 0; i13 < 7; i13++) {
      for (i14 = 0; i14 < 7; i14++) {
        f7 = 0.0F;
        for (i15 = 0; i15 < 7; i15++) {
          f7 += fx_tmp[i15][i13] * dxdot1[i14][i15];
        }

        b_fv1[i14][i13] = ((float)iv2[i14][i13] + dt * dxdot2[i14][i13]) + f7;
      }
    }

    b_fv[0][0] = 0.0F;
    b_fv[1][0] = -xnew[k][3];
    b_fv[2][0] = xnew[k][2];
    b_fv[0][1] = xnew[k][3];
    b_fv[1][1] = 0.0F;
    b_fv[2][1] = -xnew[k][1];
    b_fv[0][2] = -xnew[k][2];
    b_fv[1][2] = xnew[k][1];
    b_fv[2][2] = 0.0F;
    for (i17 = 0; i17 < 6; i17++) {
      for (i19 = 0; i19 < 7; i19++) {
        f8 = 0.0F;
        for (i21 = 0; i21 < 7; i21++) {
          f8 += b_fx_tmp[i21][i17] * b_fv1[i19][i21];
        }

        c_fx_tmp[i19][i17] = f8;
      }
    }

    for (i18 = 0; i18 < 3; i18++) {
      c_xnew[i18][0] = -xnew[k][i18 + 1];
      c_xnew[i18 + 3][0] = 0.0F;
      c_xnew[i18][1] = xnew[k][0] * (float)iv[i18][0] + b_fv[i18][0];
      c_xnew[i18 + 3][1] = 0.0F;
      c_xnew[i18][2] = xnew[k][0] * (float)iv[i18][1] + b_fv[i18][1];
      c_xnew[i18 + 3][2] = 0.0F;
      c_xnew[i18][3] = xnew[k][0] * (float)iv[i18][2] + b_fv[i18][2];
      c_xnew[i18 + 3][3] = 0.0F;
    }

    for (i20 = 0; i20 < 6; i20++) {
      c_xnew[i20][4] = iv1[i20][0];
      c_xnew[i20][5] = iv1[i20][1];
      c_xnew[i20][6] = iv1[i20][2];
    }

    for (i22 = 0; i22 < 6; i22++) {
      for (i24 = 0; i24 < 6; i24++) {
        f9 = 0.0F;
        for (i27 = 0; i27 < 7; i27++) {
          f9 += c_fx_tmp[i27][i22] * c_xnew[i24][i27];
        }

        fx[k][i24][i22] = f9;
      }
    }

    for (i23 = 0; i23 < 7; i23++) {
      for (i26 = 0; i26 < 3; i26++) {
        f10 = 0.0F;
        for (i29 = 0; i29 < 7; i29++) {
          f10 += fx_tmp[i29][i23] * dxdot1[i26 + 7][i29];
        }

        b_dt[i26][i23] = dt * dxdot2[i26 + 7][i23] + f10;
      }
    }

    for (i25 = 0; i25 < 6; i25++) {
      for (i28 = 0; i28 < 3; i28++) {
        f11 = 0.0F;
        for (i30 = 0; i30 < 7; i30++) {
          f11 += b_fx_tmp[i30][i25] * b_dt[i28][i30];
        }

        fu[k][i28][i25] = f11;
      }
    }

    xnew[k + 1][0] = x1[0] / b_y;
    xnew[k + 1][1] = x1[1] / b_y;
    xnew[k + 1][2] = x1[2] / b_y;
    xnew[k + 1][3] = x1[3] / b_y;
    xnew[k + 1][4] = x1[4];
    xnew[k + 1][5] = x1[5];
    xnew[k + 1][6] = x1[6];

    /*  Calculate the cost */
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
    /*  cumulative angular velocity cost hessian */
    /*  cumulative geodesic cost weighting */
    /*  Finds the geodesic quaternion-error cost */
    /*  quat_cost = min(1+qd'q, 1-qd'q)) where qd is the desired attitude quaternion  */
    /*  Also records the sign (+ or -) which minimizes quat_cost */
    /*  this is used when calculating the Jacobain and Hessian */
    b_xg_tmp = ((xg[0] * xnew[k][0] + xg[1] * xnew[k][1]) + xg[2] * xnew[k][2])
      + xg[3] * xnew[k][3];
    if (b_xg_tmp + 1.0F < 1.0F - b_xg_tmp) {
      b_out = b_xg_tmp + 1.0F;
      c_sign = 1;
    } else {
      b_out = 1.0F - b_xg_tmp;
      c_sign = -1;
    }

    cost_tmp[0] = xnew[k][4] - xg[4];
    cost_tmp[1] = xnew[k][5] - xg[5];
    cost_tmp[2] = xnew[k][6] - xg[6];

    /*  State cost Hessian */
    for (i31 = 0; i31 < 3; i31++) {
      cxx[k][i31][0] = (float)(-c_sign * iv[i31][0]) * b_xg_tmp;
      cxx[k][i31 + 3][0] = 0.0F;
      cxx[k][i31][1] = (float)(-c_sign * iv[i31][1]) * b_xg_tmp;
      cxx[k][i31 + 3][1] = 0.0F;
      cxx[k][i31][2] = (float)(-c_sign * iv[i31][2]) * b_xg_tmp;
      cxx[k][i31 + 3][2] = 0.0F;
    }

    for (i32 = 0; i32 < 6; i32++) {
      cxx[k][i32][3] = fv[i32][0];
      cxx[k][i32][4] = fv[i32][1];
      cxx[k][i32][5] = fv[i32][2];
    }

    /*  State cost Jacobian */
    /*  Returns skew symmetric - cross product matrix of a 3x1 vector */
    b_fv[0][0] = 0.0F;
    b_fv[1][0] = -xnew[k][3];
    b_fv[2][0] = xnew[k][2];
    b_fv[0][1] = xnew[k][3];
    b_fv[1][1] = 0.0F;
    b_fv[2][1] = -xnew[k][1];
    b_fv[0][2] = -xnew[k][2];
    b_fv[1][2] = xnew[k][1];
    b_fv[2][2] = 0.0F;
    for (i33 = 0; i33 < 3; i33++) {
      d_sign[0][i33] = (float)c_sign * -xnew[k][i33 + 1];
      d_sign[i33 + 1][0] = (float)c_sign * (xnew[k][0] * (float)iv[0][i33] +
        b_fv[0][i33]);
      d_sign[i33 + 1][1] = (float)c_sign * (xnew[k][0] * (float)iv[1][i33] +
        b_fv[1][i33]);
      d_sign[i33 + 1][2] = (float)c_sign * (xnew[k][0] * (float)iv[2][i33] +
        b_fv[2][i33]);
    }

    /*  Control cost Hessian & Jacobian */
    f12 = 0.0F;
    f13 = unew[k][0];
    f14 = unew[k][1];
    f15 = unew[k][2];
    for (i34 = 0; i34 < 3; i34++) {
      cx[k][i34] = ((d_sign[0][i34] * xg[0] + d_sign[1][i34] * xg[1]) + d_sign[2]
                    [i34] * xg[2]) + d_sign[3][i34] * xg[3];
      cx[k][i34 + 3] = (fv1[0][i34] * cost_tmp[0] + fv1[1][i34] * cost_tmp[1]) +
        fv1[2][i34] * cost_tmp[2];
      cuu[k][i34][0] = fv2[i34][0];
      cuu[k][i34][1] = fv2[i34][1];
      cuu[k][i34][2] = fv2[i34][2];
      cu[k][i34] = (fv2[0][i34] * f13 + fv2[1][i34] * f14) + fv2[2][i34] * f15;
      f12 += ((0.5F * cost_tmp[0] * fv1[i34][0] + 0.5F * cost_tmp[1] * fv1[i34]
               [1]) + 0.5F * cost_tmp[2] * fv1[i34][2]) * cost_tmp[i34];
      b_fv2[i34] = (0.5F * f13 * fv2[i34][0] + 0.5F * f14 * fv2[i34][1]) + 0.5F *
        f15 * fv2[i34][2];
    }

    *cost += (b_out + f12) + ((b_fv2[0] * unew[k][0] + b_fv2[1] * unew[k][1]) +
      b_fv2[2] * unew[k][2]);
  }

  /*  Final cost */
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
  xg_tmp = ((xg[0] * xnew[4999][0] + xg[1] * xnew[4999][1]) + xg[2] * xnew[4999]
            [2]) + xg[3] * xnew[4999][3];
  if (xg_tmp + 1.0F < 1.0F - xg_tmp) {
    out = xg_tmp + 1.0F;
    b_sign = 1;
  } else {
    out = 1.0F - xg_tmp;
    b_sign = -1;
  }

  cost_tmp[0] = xnew[4999][4] - xg[4];
  cost_tmp[1] = xnew[4999][5] - xg[5];
  cost_tmp[2] = xnew[4999][6] - xg[6];

  /*  State cost Hessian */
  c_a = -b_sign * 10;
  for (i2 = 0; i2 < 3; i2++) {
    cxx[4999][i2][0] = (float)(c_a * iv[i2][0]) * xg_tmp;
    cxx[4999][i2 + 3][0] = 0.0F;
    cxx[4999][i2][1] = (float)(c_a * iv[i2][1]) * xg_tmp;
    cxx[4999][i2 + 3][1] = 0.0F;
    cxx[4999][i2][2] = (float)(c_a * iv[i2][2]) * xg_tmp;
    cxx[4999][i2 + 3][2] = 0.0F;
  }

  for (i3 = 0; i3 < 6; i3++) {
    cxx[4999][i3][3] = iv1[i3][0];
    cxx[4999][i3][4] = iv1[i3][1];
    cxx[4999][i3][5] = iv1[i3][2];
  }

  /*  State cost Jacobian */
  /*  Returns skew symmetric - cross product matrix of a 3x1 vector */
  d_a = b_sign * 10;
  b_fv[0][0] = 0.0F;
  b_fv[1][0] = -xnew[4999][3];
  b_fv[2][0] = xnew[4999][2];
  b_fv[0][1] = xnew[4999][3];
  b_fv[1][1] = 0.0F;
  b_fv[2][1] = -xnew[4999][1];
  b_fv[0][2] = -xnew[4999][2];
  b_fv[1][2] = xnew[4999][1];
  b_fv[2][2] = 0.0F;
  for (i5 = 0; i5 < 3; i5++) {
    e_a[0][i5] = (float)d_a * -xnew[4999][i5 + 1];
    e_a[i5 + 1][0] = (float)d_a * (xnew[4999][0] * (float)iv[0][i5] + b_fv[0][i5]);
    e_a[i5 + 1][1] = (float)d_a * (xnew[4999][0] * (float)iv[1][i5] + b_fv[1][i5]);
    e_a[i5 + 1][2] = (float)d_a * (xnew[4999][0] * (float)iv[2][i5] + b_fv[2][i5]);
  }

  /*  Control cost Hessian & Jacobian */
  f5 = 0.0F;
  for (i6 = 0; i6 < 3; i6++) {
    cx[4999][i6] = ((e_a[0][i6] * xg[0] + e_a[1][i6] * xg[1]) + e_a[2][i6] * xg
                    [2]) + e_a[3][i6] * xg[3];
    cx[4999][i6 + 3] = ((float)iv[0][i6] * cost_tmp[0] + (float)iv[1][i6] *
                        cost_tmp[1]) + (float)iv[2][i6] * cost_tmp[2];
    f5 += ((0.5F * cost_tmp[0] * (float)iv[i6][0] + 0.5F * cost_tmp[1] * (float)
            iv[i6][1]) + 0.5F * cost_tmp[2] * (float)iv[i6][2]) * cost_tmp[i6];
  }

  *cost += 10.0F * out + f5;
}

/*
 * Arguments    : float fx[4999][6][6]
 *                float fu[4999][3][6]
 *                float cx[5000][6]
 *                float cu[4999][3]
 *                float cxx[5000][6][6]
 *                float cuu[4999][3][3]
 *                float lambda
 *                float u_lims[2][3]
 *                float u[4999][3]
 *                float l[4999][3]
 *                float K[4999][6][3]
 *                float dV[2]
 *                bool *diverged
 * Return Type  : void
 */
static void backwardPass(float fx[4999][6][6], float fu[4999][3][6], float cx
  [5000][6], float cu[4999][3], float cxx[5000][6][6], float cuu[4999][3][3],
  float lambda, float u_lims[2][3], float u[4999][3], float l[4999][3], float K
  [4999][6][3], float dV[2], bool *diverged)
{
  int i;
  int i1;
  float Vxx[6][6];
  int i2;
  float Vx[6];
  int k;
  bool exitg1;
  int i3;
  int i4;
  int i5;
  float f;
  float Qu_tmp[6][3];
  float Qx_tmp[6][6];
  int i6;
  int i7;
  int i8;
  float Qu[3];
  int i9;
  float f1;
  float Quu[3][3];
  float b_u_lims[3];
  float c_u_lims[3];
  float b_fv[3];
  float lk[3];
  float result;
  float Luu_data[9];
  int Luu_size[2];
  bool b_free[3];
  float b_Quu[3][3];
  int i10;
  int i11;
  float f2;
  int i12;
  float Quu_tmp[6][3];
  float u_lims_tmp;
  float f3;
  int i13;
  int i14;
  bool y;
  float Kk[6][3];
  float Qux[6][3];
  int b_k;
  bool exitg2;
  int partialTrueCount;
  int i15;
  signed char tmp_data[3];
  int i16;
  int i17;
  int i18;
  int i19;
  float Vx_tmp[3][6];
  int i20;
  int i21;
  float b_Vx_tmp[3][6];
  int i22;
  float f4;
  float c_Vx_tmp[3][6];
  int i23;
  int i24;
  float b_cx[6];
  float d_Vx_tmp[6];
  int value;
  int i25;
  int i26;
  int i27;
  int Y_size_idx_0;
  int i28;
  float f5;
  int loop_ub;
  float b_fv1[6][6];
  int i29;
  float b_cxx[6][6];
  float e_Vx_tmp[6][6];
  int i30;
  int i31;
  float f6;
  int i32;
  float b_Qx_tmp[6][6];
  int coeffs_size_idx_0;
  int i33;
  float f7;
  int b_loop_ub;
  int i34;
  float Y_data[18];
  float coeffs_data[3];
  float dV_idx_0;
  int i35;
  float dV_idx_1;
  int b_i;
  int i36;
  int X_size_idx_0;
  int j;
  int c_loop_ub;
  int i37;
  int i38;
  float coeffs;
  float b_coeffs;
  int LT_size_idx_0;
  int i39;
  int i40;
  int i41;
  int d_loop_ub;
  int i42;
  float X_data[18];
  float out;
  int e_loop_ub;
  int f_loop_ub;
  int i43;
  float LT_data[9];
  int i44;
  int c_i;
  int b_partialTrueCount;
  int d_i;
  signed char b_tmp_data[3];
  int b_j;
  int i45;
  int i46;
  int i47;
  float c_coeffs;
  float d_coeffs;
  int i48;
  int i49;
  int i50;
  int g_loop_ub;
  int i51;
  int i52;
  float b_out;

  /*  Perfoms the LQR backward pass to find the optimal controls */
  /*  Solves a quadratic program (QP) at each timestep for the optimal */
  /*  controls given the control limits */
  /*  Initialize matrices (for C code, not needed in MATLAB) */
  for (i = 0; i < 4999; i++) {
    l[i][0] = 0.0F;
    l[i][1] = 0.0F;
    l[i][2] = 0.0F;
    for (i1 = 0; i1 < 6; i1++) {
      K[i][i1][0] = 0.0F;
      K[i][i1][1] = 0.0F;
      K[i][i1][2] = 0.0F;
    }
  }

  /*  Change in cost */
  dV[0] = 0.0F;
  dV[1] = 0.0F;

  /*  Set cost-to-go Jacobian and Hessian equal to final costs */
  memcpy(&Vxx[0][0], &cxx[4999][0][0], 36U * sizeof(float));
  for (i2 = 0; i2 < 6; i2++) {
    Vx[i2] = cx[4999][i2];
  }

  *diverged = false;
  k = 0;
  exitg1 = false;
  while ((!exitg1) && (k < 4999)) {
    /*  Define cost gradients and hessians */
    for (i3 = 0; i3 < 6; i3++) {
      for (i5 = 0; i5 < 6; i5++) {
        Qx_tmp[i3][i5] = fx[4998 - k][i5][i3];
      }

      Qu_tmp[i3][0] = fu[4998 - k][0][i3];
      Qu_tmp[i3][1] = fu[4998 - k][1][i3];
      Qu_tmp[i3][2] = fu[4998 - k][2][i3];
    }

    for (i4 = 0; i4 < 3; i4++) {
      f = 0.0F;
      for (i6 = 0; i6 < 6; i6++) {
        f += Qu_tmp[i6][i4] * Vx[i6];
        f1 = 0.0F;
        for (i10 = 0; i10 < 6; i10++) {
          f1 += Qu_tmp[i10][i4] * Vxx[i6][i10];
        }

        Quu_tmp[i6][i4] = f1;
      }

      Qu[i4] = cu[4998 - k][i4] + f;
      for (i9 = 0; i9 < 3; i9++) {
        f2 = 0.0F;
        for (i12 = 0; i12 < 6; i12++) {
          f2 += Quu_tmp[i12][i4] * fu[4998 - k][i9][i12];
        }

        b_Quu[i9][i4] = cuu[4998 - k][i9][i4] + f2;
      }

      for (i11 = 0; i11 < 6; i11++) {
        f3 = 0.0F;
        for (i14 = 0; i14 < 6; i14++) {
          f3 += Quu_tmp[i14][i4] * fx[4998 - k][i11][i14];
        }

        Qux[i11][i4] = f3;
      }
    }

    /*  Regularization (for Cholesky positive definiteness) */
    /*  Solve the Quadratic program with control limits */
    i7 = (int)fminf(4999.0F, (-(float)k + 4999.0F) + 1.0F);
    for (i8 = 0; i8 < 3; i8++) {
      Quu[i8][0] = b_Quu[i8][0] + (float)iv[i8][0] * lambda;
      Quu[i8][1] = b_Quu[i8][1] + (float)iv[i8][1] * lambda;
      Quu[i8][2] = b_Quu[i8][2] + (float)iv[i8][2] * lambda;
      u_lims_tmp = u[4998 - k][i8];
      b_u_lims[i8] = u_lims[0][i8] - u_lims_tmp;
      c_u_lims[i8] = u_lims[1][i8] - u_lims_tmp;
      b_fv[i8] = -l[i7 - 1][i8];
    }

    boxQPsolve(Quu, Qu, b_u_lims, c_u_lims, b_fv, lk, &result, Luu_data,
               Luu_size, b_free);
    if (result < 2.0F) {
      *diverged = true;
      exitg1 = true;
    } else {
      /*  Solve for feedback gains in non-clamped rows of u */
      /*  (using cholesky factor of Quu) */
      for (i13 = 0; i13 < 6; i13++) {
        Kk[i13][0] = 0.0F;
        Kk[i13][1] = 0.0F;
        Kk[i13][2] = 0.0F;
      }

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
        partialTrueCount = 0;
        if (b_free[0]) {
          tmp_data[0] = 1;
          partialTrueCount = 1;
        }

        if (b_free[1]) {
          tmp_data[partialTrueCount] = 2;
          partialTrueCount++;
        }

        if (b_free[2]) {
          tmp_data[partialTrueCount] = 3;
        }

        /*  Solves the linear system AX = B for the unknown X using the Cholesky */
        /*  decomposition L of the matrix A. Where LL' = A */
        /*  X can be a vector or a matrix of size n x m */
        /*  Solution is found in O(nm) time using back-substitution */
        /*  This implementation only works for lower triangular factorisations */
        value = Luu_size[0] - 1;

        /*  Check sizes match and L is lower-triangular */
        Y_size_idx_0 = Luu_size[0];
        loop_ub = Luu_size[0];
        for (i31 = 0; i31 < 6; i31++) {
          for (i33 = 0; i33 < loop_ub; i33++) {
            Y_data[i33 + Y_size_idx_0 * i31] = 0.0F;
          }
        }

        coeffs_size_idx_0 = Luu_size[0];
        b_loop_ub = Luu_size[0];
        if (0 <= b_loop_ub - 1) {
          memset(&coeffs_data[0], 0, b_loop_ub * sizeof(float));
        }

        i35 = Luu_size[0];
        for (b_i = 0; b_i < i35; b_i++) {
          if (b_i + 1 != 1) {
            for (i36 = 0; i36 < b_i; i36++) {
              coeffs_data[i36] = Luu_data[b_i + Luu_size[0] * i36];
            }
          }

          for (j = 0; j < 6; j++) {
            if (((signed char)coeffs_size_idx_0 == 1) || (Y_size_idx_0 == 1)) {
              b_coeffs = 0.0F;
              for (i41 = 0; i41 < coeffs_size_idx_0; i41++) {
                b_coeffs += coeffs_data[i41] * Y_data[i41 + Y_size_idx_0 * j];
              }

              out = b_coeffs;
            } else {
              coeffs = 0.0F;
              for (i40 = 0; i40 < coeffs_size_idx_0; i40++) {
                coeffs += coeffs_data[i40] * Y_data[i40 + Y_size_idx_0 * j];
              }

              out = coeffs;
            }

            Y_data[b_i + Y_size_idx_0 * j] = (Qux[j][tmp_data[b_i] - 1] - out) /
              Luu_data[b_i + Luu_size[0] * b_i];
          }
        }

        /*  Solve L'X = Y for X */
        /* ======================= */
        X_size_idx_0 = Luu_size[0];
        c_loop_ub = Luu_size[0];
        for (i37 = 0; i37 < 6; i37++) {
          for (i39 = 0; i39 < c_loop_ub; i39++) {
            X_data[i39 + X_size_idx_0 * i37] = 0.0F;
          }
        }

        LT_size_idx_0 = Luu_size[1];
        d_loop_ub = Luu_size[0];
        for (i42 = 0; i42 < d_loop_ub; i42++) {
          e_loop_ub = Luu_size[1];
          for (i43 = 0; i43 < e_loop_ub; i43++) {
            LT_data[i43 + LT_size_idx_0 * i42] = Luu_data[i42 + Luu_size[0] *
              i43];
          }
        }

        coeffs_size_idx_0 = Luu_size[0];
        f_loop_ub = Luu_size[0];
        if (0 <= f_loop_ub - 1) {
          memset(&coeffs_data[0], 0, f_loop_ub * sizeof(float));
        }

        i44 = (int)(((-1.0F - (float)Luu_size[0]) + 1.0F) / -1.0F);
        for (c_i = 0; c_i < i44; c_i++) {
          d_i = value - c_i;
          if (d_i + 1 != value + 1) {
            if (d_i + 2 > Luu_size[0]) {
              i45 = 0;
              i46 = -1;
              i47 = 0;
            } else {
              i45 = d_i + 1;
              i46 = value;
              i47 = d_i + 1;
            }

            g_loop_ub = i46 - i45;
            for (i51 = 0; i51 <= g_loop_ub; i51++) {
              coeffs_data[i47 + i51] = LT_data[d_i + LT_size_idx_0 * (i45 + i51)];
            }
          }

          for (b_j = 0; b_j < 6; b_j++) {
            if (((signed char)coeffs_size_idx_0 == 1) || (X_size_idx_0 == 1)) {
              d_coeffs = 0.0F;
              for (i50 = 0; i50 < coeffs_size_idx_0; i50++) {
                d_coeffs += coeffs_data[i50] * X_data[i50 + X_size_idx_0 * b_j];
              }

              b_out = d_coeffs;
            } else {
              c_coeffs = 0.0F;
              for (i49 = 0; i49 < coeffs_size_idx_0; i49++) {
                c_coeffs += coeffs_data[i49] * X_data[i49 + X_size_idx_0 * b_j];
              }

              b_out = c_coeffs;
            }

            X_data[d_i + X_size_idx_0 * b_j] = (Y_data[d_i + Y_size_idx_0 * b_j]
              - b_out) / LT_data[d_i + LT_size_idx_0 * d_i];
          }
        }

        b_partialTrueCount = 0;
        if (b_free[0]) {
          b_tmp_data[0] = 1;
          b_partialTrueCount = 1;
        }

        if (b_free[1]) {
          b_tmp_data[b_partialTrueCount] = 2;
          b_partialTrueCount++;
        }

        if (b_free[2]) {
          b_tmp_data[b_partialTrueCount] = 3;
        }

        for (i48 = 0; i48 < 6; i48++) {
          for (i52 = 0; i52 < X_size_idx_0; i52++) {
            Kk[i48][b_tmp_data[i52] - 1] = -X_data[i52 + X_size_idx_0 * i48];
          }
        }
      }

      /*  Update Cost to Go Jacobian and Hessian */
      for (i15 = 0; i15 < 3; i15++) {
        for (i17 = 0; i17 < 6; i17++) {
          Vx_tmp[i15][i17] = Kk[i17][i15];
        }
      }

      for (i16 = 0; i16 < 6; i16++) {
        for (i19 = 0; i19 < 3; i19++) {
          b_Vx_tmp[i19][i16] = (Vx_tmp[0][i16] * b_Quu[i19][0] + Vx_tmp[1][i16] *
                                b_Quu[i19][1]) + Vx_tmp[2][i16] * b_Quu[i19][2];
        }
      }

      for (i18 = 0; i18 < 3; i18++) {
        for (i21 = 0; i21 < 6; i21++) {
          c_Vx_tmp[i18][i21] = Qux[i21][i18];
        }
      }

      for (i20 = 0; i20 < 6; i20++) {
        f4 = 0.0F;
        for (i23 = 0; i23 < 6; i23++) {
          f4 += Qx_tmp[i23][i20] * Vx[i23];
        }

        b_cx[i20] = ((cx[4998 - k][i20] + f4) + ((b_Vx_tmp[0][i20] * lk[0] +
          b_Vx_tmp[1][i20] * lk[1]) + b_Vx_tmp[2][i20] * lk[2])) + ((Vx_tmp[0]
          [i20] * Qu[0] + Vx_tmp[1][i20] * Qu[1]) + Vx_tmp[2][i20] * Qu[2]);
        d_Vx_tmp[i20] = (c_Vx_tmp[0][i20] * lk[0] + c_Vx_tmp[1][i20] * lk[1]) +
          c_Vx_tmp[2][i20] * lk[2];
      }

      for (i22 = 0; i22 < 6; i22++) {
        Vx[i22] = b_cx[i22] + d_Vx_tmp[i22];
        for (i25 = 0; i25 < 6; i25++) {
          f5 = 0.0F;
          for (i30 = 0; i30 < 6; i30++) {
            f5 += Qx_tmp[i30][i22] * Vxx[i25][i30];
          }

          b_Qx_tmp[i25][i22] = f5;
        }

        for (i28 = 0; i28 < 6; i28++) {
          f6 = 0.0F;
          for (i32 = 0; i32 < 6; i32++) {
            f6 += b_Qx_tmp[i32][i22] * fx[4998 - k][i28][i32];
          }

          b_cxx[i28][i22] = ((cxx[4998 - k][i28][i22] + f6) + ((b_Vx_tmp[0][i22]
            * Kk[i28][0] + b_Vx_tmp[1][i22] * Kk[i28][1]) + b_Vx_tmp[2][i22] *
            Kk[i28][2])) + ((Vx_tmp[0][i22] * Qux[i28][0] + Vx_tmp[1][i22] *
                             Qux[i28][1]) + Vx_tmp[2][i22] * Qux[i28][2]);
          e_Vx_tmp[i28][i22] = (c_Vx_tmp[0][i22] * Kk[i28][0] + c_Vx_tmp[1][i22]
                                * Kk[i28][1]) + c_Vx_tmp[2][i22] * Kk[i28][2];
        }
      }

      for (i24 = 0; i24 < 6; i24++) {
        for (i27 = 0; i27 < 6; i27++) {
          Vxx[i24][i27] = b_cxx[i24][i27] + e_Vx_tmp[i24][i27];
        }
      }

      for (i26 = 0; i26 < 6; i26++) {
        for (i29 = 0; i29 < 6; i29++) {
          b_fv1[i26][i29] = 0.5F * (Vxx[i26][i29] + Vxx[i29][i26]);
        }
      }

      memcpy(&Vxx[0][0], &b_fv1[0][0], 36U * sizeof(float));

      /*  Ensure Hessian is symmetric */
      /*  Record control cost change to check convergence */
      f7 = 0.0F;
      for (i34 = 0; i34 < 3; i34++) {
        f7 += ((0.5F * lk[0] * b_Quu[i34][0] + 0.5F * lk[1] * b_Quu[i34][1]) +
               0.5F * lk[2] * b_Quu[i34][2]) * lk[i34];
      }

      dV_idx_0 = dV[0] + ((lk[0] * Qu[0] + lk[1] * Qu[1]) + lk[2] * Qu[2]);
      dV_idx_1 = dV[1] + f7;
      dV[0] = dV_idx_0;
      dV[1] = dV_idx_1;

      /*  Update Control Vectors */
      l[4998 - k][0] = -lk[0];
      l[4998 - k][1] = -lk[1];
      l[4998 - k][2] = -lk[2];
      for (i38 = 0; i38 < 6; i38++) {
        K[4998 - k][i38][0] = -Kk[i38][0];
        K[4998 - k][i38][1] = -Kk[i38][1];
        K[4998 - k][i38][2] = -Kk[i38][2];
      }

      k++;
    }
  }
}

/*
 * Arguments    : float Quu[3][3]
 *                const float Qu[3]
 *                const float lower_lim[3]
 *                const float upper_lim[3]
 *                const float u0[3]
 *                float u[3]
 *                float *result
 *                float Luu_data[]
 *                int Luu_size[2]
 *                bool b_free[3]
 * Return Type  : void
 */
static void boxQPsolve(float Quu[3][3], const float Qu[3], const float
  lower_lim[3], const float upper_lim[3], const float u0[3], float u[3], float
  *result, float Luu_data[], int Luu_size[2], bool b_free[3])
{
  float old_cost;
  float z1;
  int k;
  float f;
  bool clamped[3];
  int i;
  float cost;
  float b_z1[3];
  int iter;
  float f1;
  int b_iter;
  bool exitg1;
  int b_i;
  float f2;
  bool out;
  float grad[3];
  int b_k;
  bool prev_clamped[3];
  bool exitg2;
  bool b;
  bool factorize;
  int c_k;
  bool guard1 = false;
  int trueCount;
  int b_trueCount;
  int partialTrueCount;
  signed char tmp_data[3];
  int b_partialTrueCount;
  signed char b_tmp_data[3];
  int Quu_size[2];
  int i1;
  float Quu_data[9];
  float indef;
  int i2;
  float fcnOutput;
  float grad_norm;
  float scale;
  int d_k;
  float u_idx_0;
  float absxk;
  float u_idx_1;
  float u_idx_2;
  float t;
  int c_i;
  int c_partialTrueCount;
  float grad_clamped[3];
  float delta_u[3];
  signed char c_tmp_data[3];
  int value;
  int loop_ub;
  float Y_data[3];
  int coeffs_size_idx_0;
  int b_loop_ub;
  float coeffs_data[3];
  int i3;
  int d_i;
  int i4;
  int X_size_idx_0;
  int c_loop_ub;
  float coeffs;
  float b_coeffs;
  int i5;
  int i6;
  float X_data[3];
  int LT_size_idx_0;
  float b_out;
  int d_loop_ub;
  int i7;
  int e_loop_ub;
  int f_loop_ub;
  int i8;
  float LT_data[9];
  int i9;
  int e_i;
  int i10;
  int f_i;
  int d_partialTrueCount;
  float d_tmp_data[3];
  int i11;
  float c_coeffs;
  float d_coeffs;
  int i12;
  int i13;
  int i14;
  int i15;
  float c_out;
  int g_loop_ub;
  int i16;
  float y;
  float step;
  float f3;
  float u_c_idx_0;
  float c_z1;
  float u_c_idx_1;
  float u_c_idx_2;
  float f4;
  int i17;
  float cost_c;
  float f5;
  float d_z1;
  float f6;
  int i18;

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
  Luu_size[0] = 3;
  Luu_size[1] = 3;

  /*  Placeholder to return if Luu not assigned */
  /*  Initialize scalars */
  old_cost = 0.0F;
  *result = 0.0F;

  /*  Solver Options */
  /*  max iterations */
  /*  min norm of non-clamped gradient */
  /*  min relative improvement */
  /*  factor for decreasing stepsize */
  /*  min stepsize for linesearch */
  /*  Armijo tolerance (fraction of linear improvement required) */
  /*  Initial controls */
  /*  Returns array x with all values clamped between lower and upper */
  /*  Initial cost value */
  z1 = 0.0F;
  for (k = 0; k < 3; k++) {
    clamped[k] = false;
    b_free[k] = true;
    Luu_data[3 * k] = 0.0F;
    Luu_data[3 * k + 1] = 0.0F;
    Luu_data[3 * k + 2] = 0.0F;
    f1 = fmaxf(lower_lim[k], fminf(upper_lim[k], u0[k]));
    b_z1[k] = f1;
    u[k] = f1;
    z1 += f1 * Qu[k];
  }

  f = 0.0F;
  for (i = 0; i < 3; i++) {
    f += ((0.5F * b_z1[0] * Quu[i][0] + 0.5F * b_z1[1] * Quu[i][1]) + 0.5F *
          b_z1[2] * Quu[i][2]) * b_z1[i];
  }

  cost = z1 + f;

  /*  Start optimisation */
  iter = 1;
  b_iter = 1;
  exitg1 = false;
  while ((!exitg1) && (b_iter - 1 < 100)) {
    iter = b_iter;
    if (*result != 0.0F) {
      exitg1 = true;
    } else {
      /*  Check relative cost change for convergence */
      if ((b_iter > 1) && (old_cost - cost < 1.0E-8F * fabsf(old_cost))) {
        *result = 4.0F;
        exitg1 = true;
      } else {
        old_cost = cost;

        /*  Gradient of cost function */
        /*  Find clamped controls */
        for (b_i = 0; b_i < 3; b_i++) {
          f2 = Qu[b_i] + ((Quu[0][b_i] * u[0] + Quu[1][b_i] * u[1]) + Quu[2][b_i]
                          * u[2]);
          grad[b_i] = f2;
          prev_clamped[b_i] = clamped[b_i];
          b = false;
          clamped[b_i] = false;
          if ((u[b_i] == lower_lim[b_i]) && (f2 > 0.0F)) {
            b = true;
            clamped[b_i] = true;
          }

          if ((u[b_i] == upper_lim[b_i]) && (f2 < 0.0F)) {
            b = true;
            clamped[b_i] = true;
          }

          b_free[b_i] = !b;
        }

        /*  Check if all controls clamped */
        out = true;
        b_k = 0;
        exitg2 = false;
        while ((!exitg2) && (b_k < 3)) {
          if (!clamped[b_k]) {
            out = false;
            exitg2 = true;
          } else {
            b_k++;
          }
        }

        if (out) {
          *result = 6.0F;
          exitg1 = true;
        } else {
          /*  Cholesky factorize if clamped controls have changed */
          if (b_iter == 1) {
            factorize = true;
          } else {
            factorize = false;
            c_k = 0;
            exitg2 = false;
            while ((!exitg2) && (c_k < 3)) {
              if (prev_clamped[c_k] != clamped[c_k]) {
                factorize = true;
                exitg2 = true;
              } else {
                c_k++;
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

            partialTrueCount = 0;
            if (b_free[0]) {
              tmp_data[0] = 1;
              partialTrueCount = 1;
            }

            if (b_free[1]) {
              tmp_data[partialTrueCount] = 2;
              partialTrueCount++;
            }

            if (b_free[2]) {
              tmp_data[partialTrueCount] = 3;
            }

            Quu_size[0] = trueCount;
            Quu_size[1] = trueCount;
            for (i1 = 0; i1 < trueCount; i1++) {
              for (i2 = 0; i2 < trueCount; i2++) {
                Quu_data[i2 + trueCount * i1] = Quu[tmp_data[i1] - 1]
                  [tmp_data[i2] - 1];
              }
            }

            chol_free(Quu_data, Quu_size, Luu_data, Luu_size, &indef);
            if (indef != 0.0F) {
              *result = 1.0F;
              exitg1 = true;
            } else {
              guard1 = true;
            }
          } else {
            guard1 = true;
          }

          if (guard1) {
            /*  check gradient-norm of free controls */
            b_trueCount = 0;
            if (b_free[0]) {
              b_trueCount = 1;
            }

            if (b_free[1]) {
              b_trueCount++;
            }

            if (b_free[2]) {
              b_trueCount++;
            }

            b_partialTrueCount = 0;
            if (b_free[0]) {
              b_tmp_data[0] = 1;
              b_partialTrueCount = 1;
            }

            if (b_free[1]) {
              b_tmp_data[b_partialTrueCount] = 2;
              b_partialTrueCount++;
            }

            if (b_free[2]) {
              b_tmp_data[b_partialTrueCount] = 3;
            }

            if (b_trueCount == 0) {
              grad_norm = 0.0F;
            } else {
              fcnOutput = 0.0F;
              if (b_trueCount == 1) {
                grad_norm = fabsf(grad[b_tmp_data[0] - 1]);
              } else {
                scale = 1.29246971E-26F;
                for (d_k = 0; d_k < b_trueCount; d_k++) {
                  absxk = fabsf(grad[b_tmp_data[d_k] - 1]);
                  if (absxk > scale) {
                    t = scale / absxk;
                    fcnOutput = fcnOutput * t * t + 1.0F;
                    scale = absxk;
                  } else {
                    t = absxk / scale;
                    fcnOutput += t * t;
                  }
                }

                grad_norm = scale * sqrtf(fcnOutput);
              }
            }

            if (grad_norm < 1.0E-8F) {
              *result = 5.0F;
              exitg1 = true;
            } else {
              /*  get search direction */
              u_idx_0 = u[0] * (float)clamped[0];
              u_idx_1 = u[1] * (float)clamped[1];
              u_idx_2 = u[2] * (float)clamped[2];
              for (c_i = 0; c_i < 3; c_i++) {
                grad_clamped[c_i] = Qu[c_i] + ((Quu[0][c_i] * u_idx_0 + Quu[1]
                  [c_i] * u_idx_1) + Quu[2][c_i] * u_idx_2);
                delta_u[c_i] = 0.0F;
              }

              c_partialTrueCount = 0;
              if (b_free[0]) {
                c_tmp_data[0] = 1;
                c_partialTrueCount = 1;
              }

              if (b_free[1]) {
                c_tmp_data[c_partialTrueCount] = 2;
                c_partialTrueCount++;
              }

              if (b_free[2]) {
                c_tmp_data[c_partialTrueCount] = 3;
              }

              /*  Solves the linear system AX = B for the unknown X using the Cholesky */
              /*  decomposition L of the matrix A. Where LL' = A */
              /*  X can be a vector or a matrix of size n x m */
              /*  Solution is found in O(nm) time using back-substitution */
              /*  This implementation only works for lower triangular factorisations */
              value = Luu_size[0] - 1;

              /*  Check sizes match and L is lower-triangular */
              loop_ub = Luu_size[0];
              if (0 <= loop_ub - 1) {
                memset(&Y_data[0], 0, loop_ub * sizeof(float));
              }

              coeffs_size_idx_0 = Luu_size[0];
              b_loop_ub = Luu_size[0];
              if (0 <= b_loop_ub - 1) {
                memset(&coeffs_data[0], 0, b_loop_ub * sizeof(float));
              }

              i3 = Luu_size[0];
              for (d_i = 0; d_i < i3; d_i++) {
                if (d_i + 1 != 1) {
                  for (i4 = 0; i4 < d_i; i4++) {
                    coeffs_data[i4] = Luu_data[d_i + Luu_size[0] * i4];
                  }
                }

                if (((signed char)coeffs_size_idx_0 == 1) || (Luu_size[0] == 1))
                {
                  b_coeffs = 0.0F;
                  for (i6 = 0; i6 < coeffs_size_idx_0; i6++) {
                    b_coeffs += coeffs_data[i6] * Y_data[i6];
                  }

                  b_out = b_coeffs;
                } else {
                  coeffs = 0.0F;
                  for (i5 = 0; i5 < coeffs_size_idx_0; i5++) {
                    coeffs += coeffs_data[i5] * Y_data[i5];
                  }

                  b_out = coeffs;
                }

                Y_data[d_i] = (grad_clamped[c_tmp_data[d_i] - 1] - b_out) /
                  Luu_data[d_i + Luu_size[0] * d_i];
              }

              /*  Solve L'X = Y for X */
              /* ======================= */
              X_size_idx_0 = Luu_size[0];
              c_loop_ub = Luu_size[0];
              if (0 <= c_loop_ub - 1) {
                memset(&X_data[0], 0, c_loop_ub * sizeof(float));
              }

              LT_size_idx_0 = Luu_size[1];
              d_loop_ub = Luu_size[0];
              for (i7 = 0; i7 < d_loop_ub; i7++) {
                e_loop_ub = Luu_size[1];
                for (i8 = 0; i8 < e_loop_ub; i8++) {
                  LT_data[i8 + LT_size_idx_0 * i7] = Luu_data[i7 + Luu_size[0] *
                    i8];
                }
              }

              coeffs_size_idx_0 = Luu_size[0];
              f_loop_ub = Luu_size[0];
              if (0 <= f_loop_ub - 1) {
                memset(&coeffs_data[0], 0, f_loop_ub * sizeof(float));
              }

              i9 = (int)(((-1.0F - (float)Luu_size[0]) + 1.0F) / -1.0F);
              for (e_i = 0; e_i < i9; e_i++) {
                f_i = value - e_i;
                if (f_i + 1 != value + 1) {
                  if (f_i + 2 > Luu_size[0]) {
                    i11 = 0;
                    i12 = -1;
                    i15 = 0;
                  } else {
                    i11 = f_i + 1;
                    i12 = value;
                    i15 = f_i + 1;
                  }

                  g_loop_ub = i12 - i11;
                  for (i16 = 0; i16 <= g_loop_ub; i16++) {
                    coeffs_data[i15 + i16] = LT_data[f_i + LT_size_idx_0 * (i11
                      + i16)];
                  }
                }

                if (((signed char)coeffs_size_idx_0 == 1) || (X_size_idx_0 == 1))
                {
                  d_coeffs = 0.0F;
                  for (i14 = 0; i14 < coeffs_size_idx_0; i14++) {
                    d_coeffs += coeffs_data[i14] * X_data[i14];
                  }

                  c_out = d_coeffs;
                } else {
                  c_coeffs = 0.0F;
                  for (i13 = 0; i13 < coeffs_size_idx_0; i13++) {
                    c_coeffs += coeffs_data[i13] * X_data[i13];
                  }

                  c_out = c_coeffs;
                }

                X_data[f_i] = (Y_data[f_i] - c_out) / LT_data[f_i +
                  LT_size_idx_0 * f_i];
              }

              for (i10 = 0; i10 < X_size_idx_0; i10++) {
                d_tmp_data[i10] = -X_data[i10] - u[c_tmp_data[i10] - 1];
              }

              d_partialTrueCount = 0;

              /*  cholesky solver */
              /*  check projected change in cost is a reduction */
              if (b_free[0]) {
                delta_u[0] = d_tmp_data[0];
                d_partialTrueCount = 1;
              }

              if (b_free[1]) {
                delta_u[1] = d_tmp_data[d_partialTrueCount];
                d_partialTrueCount++;
              }

              if (b_free[2]) {
                delta_u[2] = d_tmp_data[d_partialTrueCount];
              }

              y = (delta_u[0] * grad[0] + delta_u[1] * grad[1]) + delta_u[2] *
                grad[2];
              if (y >= 0.0F) {
                /*  (should not happen) */
                exitg1 = true;
              } else {
                /*  Armijo linesearch */
                step = 1.0F;

                /*  Returns array x with all values clamped between lower and upper */
                f3 = fmaxf(lower_lim[0], fminf(upper_lim[0], u[0] + delta_u[0]));
                b_z1[0] = f3;
                u_c_idx_0 = f3;
                c_z1 = f3 * Qu[0];
                f3 = fmaxf(lower_lim[1], fminf(upper_lim[1], u[1] + delta_u[1]));
                b_z1[1] = f3;
                u_c_idx_1 = f3;
                c_z1 += f3 * Qu[1];
                f3 = fmaxf(lower_lim[2], fminf(upper_lim[2], u[2] + delta_u[2]));
                b_z1[2] = f3;
                u_c_idx_2 = f3;
                c_z1 += f3 * Qu[2];
                f4 = 0.0F;
                for (i17 = 0; i17 < 3; i17++) {
                  f4 += ((0.5F * b_z1[0] * Quu[i17][0] + 0.5F * b_z1[1] *
                          Quu[i17][1]) + 0.5F * f3 * Quu[i17][2]) * b_z1[i17];
                }

                cost_c = c_z1 + f4;
                exitg2 = false;
                while ((!exitg2) && ((cost_c - cost) / (step * y) < 0.1F)) {
                  step *= 0.6F;

                  /*  Returns array x with all values clamped between lower and upper */
                  f5 = fmaxf(lower_lim[0], fminf(upper_lim[0], u[0] + step *
                              delta_u[0]));
                  b_z1[0] = f5;
                  u_c_idx_0 = f5;
                  d_z1 = f5 * Qu[0];
                  f5 = fmaxf(lower_lim[1], fminf(upper_lim[1], u[1] + step *
                              delta_u[1]));
                  b_z1[1] = f5;
                  u_c_idx_1 = f5;
                  d_z1 += f5 * Qu[1];
                  f5 = fmaxf(lower_lim[2], fminf(upper_lim[2], u[2] + step *
                              delta_u[2]));
                  b_z1[2] = f5;
                  u_c_idx_2 = f5;
                  d_z1 += f5 * Qu[2];
                  f6 = 0.0F;
                  for (i18 = 0; i18 < 3; i18++) {
                    f6 += ((0.5F * b_z1[0] * Quu[i18][0] + 0.5F * b_z1[1] *
                            Quu[i18][1]) + 0.5F * f5 * Quu[i18][2]) * b_z1[i18];
                  }

                  cost_c = d_z1 + f6;
                  if (step < 1.0E-20F) {
                    *result = 3.0F;
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
    *result = 2.0F;
  }
}

/*
 * Arguments    : const float A_data[]
 *                const int A_size[2]
 *                float L_data[]
 *                int L_size[2]
 *                float *fail
 * Return Type  : void
 */
static void chol_free(const float A_data[], const int A_size[2], float L_data[],
                      int L_size[2], float *fail)
{
  int A_size_idx_0;
  int A_size_idx_1;
  int loop_ub;
  int i;
  int n;
  int b_loop_ub;
  int info;
  int i1;
  int b_info;
  float b_A_data[9];
  int j;
  bool  exitg1;
  int i2;
  int i3;
  int idxAjj;
  float ssq;
  int jmax;
  int ix;
  int b_j;
  float ajj;
  int iy;
  int k;
  int b_i;
  int c_loop_ub;
  int d_loop_ub;
  int nmj;
  int i4;
  int ia0;
  int idxAjp1j;
  int i5;
  int b_ix;
  int i6;
  float c_A_data[9];
  float a;
  int i7;
  int i8;
  int iac;
  int i9;
  int b_k;
  float c;
  int b_iy;
  int i10;
  int ia;

  /*  Wrapper for MATLAB chol for use with auto coder */
  /*  Inputs: */
  /* =========== */
  /*  A - positive semi-definite matrix */
  A_size_idx_0 = A_size[0];
  A_size_idx_1 = A_size[1];
  loop_ub = A_size[1];
  for (i = 0; i < loop_ub; i++) {
    b_loop_ub = A_size[0];
    for (i1 = 0; i1 < b_loop_ub; i1++) {
      b_A_data[i1 + A_size_idx_0 * i] = A_data[i1 + A_size[0] * i];
    }
  }

  n = A_size[1];
  info = 0;
  if (A_size[1] != 0) {
    b_info = 0;
    j = 0;
    exitg1 = false;
    while ((!exitg1) && (j <= n - 1)) {
      idxAjj = j + j * n;
      ssq = 0.0F;
      if (j >= 1) {
        ix = j;
        iy = j;
        for (k = 0; k < j; k++) {
          ssq += b_A_data[ix] * b_A_data[iy];
          ix += n;
          iy += n;
        }
      }

      ajj = b_A_data[idxAjj] - ssq;
      if (ajj > 0.0F) {
        ajj = sqrtf(ajj);
        b_A_data[idxAjj] = ajj;
        if (j + 1 < n) {
          nmj = (n - j) - 1;
          ia0 = j + 2;
          idxAjp1j = idxAjj + 2;
          if ((nmj != 0) && (j != 0)) {
            b_ix = j;
            i7 = (j + n * (j - 1)) + 2;
            for (iac = ia0; n < 0 ? iac >= i7 : iac <= i7; iac += n) {
              c = -b_A_data[b_ix];
              b_iy = idxAjj + 1;
              i10 = (iac + nmj) - 1;
              for (ia = iac; ia <= i10; ia++) {
                b_A_data[b_iy] += b_A_data[ia - 1] * c;
                b_iy++;
              }

              b_ix += n;
            }
          }

          a = 1.0F / ajj;
          i8 = idxAjj + nmj;
          for (b_k = idxAjp1j; b_k <= i8 + 1; b_k++) {
            b_A_data[b_k - 1] *= a;
          }
        }

        j++;
      } else {
        b_A_data[idxAjj] = ajj;
        b_info = j + 1;
        exitg1 = true;
      }
    }

    info = b_info;
    if (b_info == 0) {
      jmax = A_size[1];
    } else {
      jmax = b_info - 1;
    }

    for (b_j = 2; b_j <= jmax; b_j++) {
      for (b_i = 0; b_i <= b_j - 2; b_i++) {
        b_A_data[b_i + A_size_idx_0 * (b_j - 1)] = 0.0F;
      }
    }

    if (1 > jmax) {
      c_loop_ub = 0;
      d_loop_ub = 0;
    } else {
      c_loop_ub = jmax;
      d_loop_ub = jmax;
    }

    for (i4 = 0; i4 < d_loop_ub; i4++) {
      for (i5 = 0; i5 < c_loop_ub; i5++) {
        c_A_data[i5 + c_loop_ub * i4] = b_A_data[i5 + A_size_idx_0 * i4];
      }
    }

    A_size_idx_0 = c_loop_ub;
    A_size_idx_1 = d_loop_ub;
    for (i6 = 0; i6 < d_loop_ub; i6++) {
      for (i9 = 0; i9 < c_loop_ub; i9++) {
        b_A_data[i9 + c_loop_ub * i6] = c_A_data[i9 + c_loop_ub * i6];
      }
    }
  }

  L_size[0] = A_size_idx_0;
  L_size[1] = A_size_idx_1;
  for (i2 = 0; i2 < A_size_idx_1; i2++) {
    for (i3 = 0; i3 < A_size_idx_0; i3++) {
      L_data[i3 + A_size_idx_0 * i2] = b_A_data[i3 + A_size_idx_0 * i2];
    }
  }

  *fail = (float)info;
}

/*
 * Arguments    : float x[5000][7]
 *                const float xg[7]
 *                float u[4999][3]
 *                float K[4999][6][3]
 *                float u_lims[2][3]
 *                float dt
 *                float B_ECI[5000][3]
 *                float Js[6][3]
 *                float xnew[5000][7]
 *                float unew[4999][3]
 *                float fx[4999][6][6]
 *                float fu[4999][3][6]
 *                float cx[5000][6]
 *                float cu[4999][3]
 *                float cxx[5000][6][6]
 *                float cuu[4999][3][3]
 *                float *cost
 * Return Type  : void
 */
static void forwardRollout(float x[5000][7], const float xg[7], float u[4999][3],
  float K[4999][6][3], float u_lims[2][3], float dt, float B_ECI[5000][3], float
  Js[6][3], float xnew[5000][7], float unew[4999][3], float fx[4999][6][6],
  float fu[4999][3][6], float cx[5000][6], float cu[4999][3], float cxx[5000][6]
  [6], float cuu[4999][3][3], float *cost)
{
  int i;
  int i1;
  float a;
  float b_a;
  int k;
  float dx[6];
  float q[4][4];
  float xg_tmp;
  float out;
  int b_sign;
  float cost_tmp[3];
  int c_a;
  int i2;
  int i3;
  float q_error;
  int d_a;
  int i4;
  float b_fv[3][3];
  float y;
  float f;
  float f1;
  float b_q_error[4];
  float f2;
  float f3;
  float f4;
  int i5;
  float e_a[4][3];
  float f5;
  int i6;
  int b_k;
  float f6;
  int i7;
  float b[7];
  float dxdot1[10][7];
  int i8;
  float b_xnew[7];
  float dxdot2[10][7];
  int i9;
  float x1[7];
  float b_y;
  int i10;
  int i11;
  float fx_tmp[7][7];
  int i12;
  int i13;
  float b_fx_tmp[7][6];
  int i14;
  float f7;
  int i15;
  float b_fv1[7][7];
  int i16;
  int i17;
  int i18;
  int i19;
  int i20;
  float c_xnew[6][7];
  float f8;
  int i21;
  int i22;
  float c_fx_tmp[7][6];
  int i23;
  int i24;
  int i25;
  int i26;
  float f9;
  int i27;
  int i28;
  float f10;
  int i29;
  float f11;
  int i30;
  float b_dt[3][7];
  float b_xg_tmp;
  float b_out;
  int c_sign;
  int i31;
  int i32;
  int i33;
  float d_sign[4][3];
  float f12;
  float f13;
  float f14;
  float f15;
  int i34;
  float b_fv2[3];

  /*  Uses an rk method to roll out a trajectory */
  /*  Returns the new trajectory, cost and the derivatives along the trajectory */
  /*  If feed-forward and feed-back controls l, K are non-zero */
  /*  then the unew returned is the new control sequeunce. Otherwise unew = u */
  /*  Sizes */
  /*  Initialize outputs (Don't overwrite x or u) */
  memset(&xnew[0][0], 0, 35000U * sizeof(float));
  for (i = 0; i < 4999; i++) {
    unew[i][0] = 0.0F;
    unew[i][1] = 0.0F;
    unew[i][2] = 0.0F;
  }

  *cost = 0.0F;
  for (i1 = 0; i1 < 7; i1++) {
    xnew[0][i1] = x[0][i1];
  }

  a = 0.5F * dt;
  b_a = 0.5F * dt * dt;
  for (k = 0; k < 4999; k++) {
    /*  Find the state error vector dx */
    dx[3] = xnew[k][4] - x[k][4];
    dx[4] = xnew[k][5] - x[k][5];
    dx[5] = xnew[k][6] - x[k][6];

    /*  Calculate error between qk and q_nom */
    /*  Defined as conj(q_nom)*qnew */
    /*  Returns error as Rodrigues parameters (3x1) */
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
    q[0][0] = x[k][0];
    q[1][0] = x[k][1];
    q[2][0] = x[k][2];
    q[3][0] = x[k][3];
    q[0][1] = -x[k][1];
    q[1][1] = x[k][0];
    q[2][1] = x[k][3];
    q[3][1] = -x[k][2];
    q[0][2] = -x[k][2];
    q[1][2] = -x[k][3];
    q[2][2] = x[k][0];
    q[3][2] = x[k][1];
    q[0][3] = -x[k][3];
    q[1][3] = x[k][2];
    q[2][3] = -x[k][1];
    q[3][3] = x[k][0];
    q_error = 0.0F;
    for (i4 = 0; i4 < 4; i4++) {
      f = ((q[0][i4] * xnew[k][0] + q[1][i4] * xnew[k][1]) + q[2][i4] * xnew[k]
           [2]) + q[3][i4] * xnew[k][3];
      b_q_error[i4] = f;
      q_error += f * f;
    }

    y = sqrtf(q_error);
    f1 = b_q_error[0] / y;
    b_q_error[0] = f1;
    f2 = b_q_error[1] / y;
    b_q_error[1] = f2;
    f3 = b_q_error[2] / y;
    b_q_error[2] = f3;
    f4 = b_q_error[3] / y;
    b_q_error[3] = f4;

    /*  re-normalize */
    /*  inverse Cayley Map */
    dx[0] = f2 / f1;
    dx[1] = f3 / f1;
    dx[2] = f4 / f1;

    /*  Find the new control and ensure it is within the limits */
    for (b_k = 0; b_k < 3; b_k++) {
      f6 = 0.0F;
      for (i7 = 0; i7 < 6; i7++) {
        f6 += K[k][i7][b_k] * dx[i7];
      }

      unew[k][b_k] = fminf(u_lims[1][b_k], fmaxf(u_lims[0][b_k], u[k][b_k] - f6));
    }

    /*  Step the dynamics forward */
    /*  Steps the dynamics forward using a 2nd order rk-method */
    /*  Returns: new state, discrete time Jacobians */
    /*  Step dynamics */
    /*  Explicit midpoint step from x_k to x_{k+1} */
    satellite_dynamics(*(float (*)[7])&xnew[k][0], *(float (*)[3])&unew[k][0],
                       *(float (*)[3])&B_ECI[k][0], Js, b, dxdot1);
    for (i8 = 0; i8 < 7; i8++) {
      b_xnew[i8] = xnew[k][i8] + a * b[i8];
    }

    satellite_dynamics(b_xnew, *(float (*)[3])&unew[k][0], *(float (*)[3])&
                       B_ECI[k][0], Js, b, dxdot2);
    for (i9 = 0; i9 < 7; i9++) {
      x1[i9] = xnew[k][i9] + dt * b[i9];
    }

    /*  Re-normalize the quaternion */
    b_y = sqrtf(((x1[0] * x1[0] + x1[1] * x1[1]) + x1[2] * x1[2]) + x1[3] * x1[3]);

    /*  Continuous time Jacobians */
    /*  Discrete time Jacobians */
    /*  First form the state error Jacobians E(x_k) & E(x_k+1) */
    /*  Returns skew symmetric - cross product matrix of a 3x1 vector */
    /*  Returns skew symmetric - cross product matrix of a 3x1 vector */
    for (i10 = 0; i10 < 7; i10++) {
      for (i11 = 0; i11 < 7; i11++) {
        fx_tmp[i10][i11] = b_a * dxdot2[i10][i11];
      }
    }

    b_fv[0][0] = 0.0F;
    b_fv[1][0] = -x1[3];
    b_fv[2][0] = x1[2];
    b_fv[0][1] = x1[3];
    b_fv[1][1] = 0.0F;
    b_fv[2][1] = -x1[1];
    b_fv[0][2] = -x1[2];
    b_fv[1][2] = x1[1];
    b_fv[2][2] = 0.0F;
    for (i12 = 0; i12 < 3; i12++) {
      b_fx_tmp[0][i12] = -x1[i12 + 1];
      b_fx_tmp[0][i12 + 3] = 0.0F;
      b_fx_tmp[i12 + 1][0] = x1[0] * (float)iv[0][i12] + b_fv[0][i12];
      b_fx_tmp[i12 + 1][3] = 0.0F;
      b_fx_tmp[i12 + 1][1] = x1[0] * (float)iv[1][i12] + b_fv[1][i12];
      b_fx_tmp[i12 + 1][4] = 0.0F;
      b_fx_tmp[i12 + 1][2] = x1[0] * (float)iv[2][i12] + b_fv[2][i12];
      b_fx_tmp[i12 + 1][5] = 0.0F;
      for (i16 = 0; i16 < 6; i16++) {
        b_fx_tmp[i12 + 4][i16] = iv1[i16][i12];
      }
    }

    for (i13 = 0; i13 < 7; i13++) {
      for (i14 = 0; i14 < 7; i14++) {
        f7 = 0.0F;
        for (i15 = 0; i15 < 7; i15++) {
          f7 += fx_tmp[i15][i13] * dxdot1[i14][i15];
        }

        b_fv1[i14][i13] = ((float)iv2[i14][i13] + dt * dxdot2[i14][i13]) + f7;
      }
    }

    b_fv[0][0] = 0.0F;
    b_fv[1][0] = -xnew[k][3];
    b_fv[2][0] = xnew[k][2];
    b_fv[0][1] = xnew[k][3];
    b_fv[1][1] = 0.0F;
    b_fv[2][1] = -xnew[k][1];
    b_fv[0][2] = -xnew[k][2];
    b_fv[1][2] = xnew[k][1];
    b_fv[2][2] = 0.0F;
    for (i17 = 0; i17 < 6; i17++) {
      for (i19 = 0; i19 < 7; i19++) {
        f8 = 0.0F;
        for (i21 = 0; i21 < 7; i21++) {
          f8 += b_fx_tmp[i21][i17] * b_fv1[i19][i21];
        }

        c_fx_tmp[i19][i17] = f8;
      }
    }

    for (i18 = 0; i18 < 3; i18++) {
      c_xnew[i18][0] = -xnew[k][i18 + 1];
      c_xnew[i18 + 3][0] = 0.0F;
      c_xnew[i18][1] = xnew[k][0] * (float)iv[i18][0] + b_fv[i18][0];
      c_xnew[i18 + 3][1] = 0.0F;
      c_xnew[i18][2] = xnew[k][0] * (float)iv[i18][1] + b_fv[i18][1];
      c_xnew[i18 + 3][2] = 0.0F;
      c_xnew[i18][3] = xnew[k][0] * (float)iv[i18][2] + b_fv[i18][2];
      c_xnew[i18 + 3][3] = 0.0F;
    }

    for (i20 = 0; i20 < 6; i20++) {
      c_xnew[i20][4] = iv1[i20][0];
      c_xnew[i20][5] = iv1[i20][1];
      c_xnew[i20][6] = iv1[i20][2];
    }

    for (i22 = 0; i22 < 6; i22++) {
      for (i24 = 0; i24 < 6; i24++) {
        f9 = 0.0F;
        for (i27 = 0; i27 < 7; i27++) {
          f9 += c_fx_tmp[i27][i22] * c_xnew[i24][i27];
        }

        fx[k][i24][i22] = f9;
      }
    }

    for (i23 = 0; i23 < 7; i23++) {
      for (i26 = 0; i26 < 3; i26++) {
        f10 = 0.0F;
        for (i29 = 0; i29 < 7; i29++) {
          f10 += fx_tmp[i29][i23] * dxdot1[i26 + 7][i29];
        }

        b_dt[i26][i23] = dt * dxdot2[i26 + 7][i23] + f10;
      }
    }

    for (i25 = 0; i25 < 6; i25++) {
      for (i28 = 0; i28 < 3; i28++) {
        f11 = 0.0F;
        for (i30 = 0; i30 < 7; i30++) {
          f11 += b_fx_tmp[i30][i25] * b_dt[i28][i30];
        }

        fu[k][i28][i25] = f11;
      }
    }

    xnew[k + 1][0] = x1[0] / b_y;
    xnew[k + 1][1] = x1[1] / b_y;
    xnew[k + 1][2] = x1[2] / b_y;
    xnew[k + 1][3] = x1[3] / b_y;
    xnew[k + 1][4] = x1[4];
    xnew[k + 1][5] = x1[5];
    xnew[k + 1][6] = x1[6];

    /*  Calculate the cost */
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
    /*  cumulative angular velocity cost hessian */
    /*  cumulative geodesic cost weighting */
    /*  Finds the geodesic quaternion-error cost */
    /*  quat_cost = min(1+qd'q, 1-qd'q)) where qd is the desired attitude quaternion  */
    /*  Also records the sign (+ or -) which minimizes quat_cost */
    /*  this is used when calculating the Jacobain and Hessian */
    b_xg_tmp = ((xg[0] * xnew[k][0] + xg[1] * xnew[k][1]) + xg[2] * xnew[k][2])
      + xg[3] * xnew[k][3];
    if (b_xg_tmp + 1.0F < 1.0F - b_xg_tmp) {
      b_out = b_xg_tmp + 1.0F;
      c_sign = 1;
    } else {
      b_out = 1.0F - b_xg_tmp;
      c_sign = -1;
    }

    cost_tmp[0] = xnew[k][4] - xg[4];
    cost_tmp[1] = xnew[k][5] - xg[5];
    cost_tmp[2] = xnew[k][6] - xg[6];

    /*  State cost Hessian */
    for (i31 = 0; i31 < 3; i31++) {
      cxx[k][i31][0] = (float)(-c_sign * iv[i31][0]) * b_xg_tmp;
      cxx[k][i31 + 3][0] = 0.0F;
      cxx[k][i31][1] = (float)(-c_sign * iv[i31][1]) * b_xg_tmp;
      cxx[k][i31 + 3][1] = 0.0F;
      cxx[k][i31][2] = (float)(-c_sign * iv[i31][2]) * b_xg_tmp;
      cxx[k][i31 + 3][2] = 0.0F;
    }

    for (i32 = 0; i32 < 6; i32++) {
      cxx[k][i32][3] = fv[i32][0];
      cxx[k][i32][4] = fv[i32][1];
      cxx[k][i32][5] = fv[i32][2];
    }

    /*  State cost Jacobian */
    /*  Returns skew symmetric - cross product matrix of a 3x1 vector */
    b_fv[0][0] = 0.0F;
    b_fv[1][0] = -xnew[k][3];
    b_fv[2][0] = xnew[k][2];
    b_fv[0][1] = xnew[k][3];
    b_fv[1][1] = 0.0F;
    b_fv[2][1] = -xnew[k][1];
    b_fv[0][2] = -xnew[k][2];
    b_fv[1][2] = xnew[k][1];
    b_fv[2][2] = 0.0F;
    for (i33 = 0; i33 < 3; i33++) {
      d_sign[0][i33] = (float)c_sign * -xnew[k][i33 + 1];
      d_sign[i33 + 1][0] = (float)c_sign * (xnew[k][0] * (float)iv[0][i33] +
        b_fv[0][i33]);
      d_sign[i33 + 1][1] = (float)c_sign * (xnew[k][0] * (float)iv[1][i33] +
        b_fv[1][i33]);
      d_sign[i33 + 1][2] = (float)c_sign * (xnew[k][0] * (float)iv[2][i33] +
        b_fv[2][i33]);
    }

    /*  Control cost Hessian & Jacobian */
    f12 = 0.0F;
    f13 = unew[k][0];
    f14 = unew[k][1];
    f15 = unew[k][2];
    for (i34 = 0; i34 < 3; i34++) {
      cx[k][i34] = ((d_sign[0][i34] * xg[0] + d_sign[1][i34] * xg[1]) + d_sign[2]
                    [i34] * xg[2]) + d_sign[3][i34] * xg[3];
      cx[k][i34 + 3] = (fv1[0][i34] * cost_tmp[0] + fv1[1][i34] * cost_tmp[1]) +
        fv1[2][i34] * cost_tmp[2];
      cuu[k][i34][0] = fv2[i34][0];
      cuu[k][i34][1] = fv2[i34][1];
      cuu[k][i34][2] = fv2[i34][2];
      cu[k][i34] = (fv2[0][i34] * f13 + fv2[1][i34] * f14) + fv2[2][i34] * f15;
      f12 += ((0.5F * cost_tmp[0] * fv1[i34][0] + 0.5F * cost_tmp[1] * fv1[i34]
               [1]) + 0.5F * cost_tmp[2] * fv1[i34][2]) * cost_tmp[i34];
      b_fv2[i34] = (0.5F * f13 * fv2[i34][0] + 0.5F * f14 * fv2[i34][1]) + 0.5F *
        f15 * fv2[i34][2];
    }

    *cost += (b_out + f12) + ((b_fv2[0] * unew[k][0] + b_fv2[1] * unew[k][1]) +
      b_fv2[2] * unew[k][2]);
  }

  /*  Final cost */
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
  xg_tmp = ((xg[0] * xnew[4999][0] + xg[1] * xnew[4999][1]) + xg[2] * xnew[4999]
            [2]) + xg[3] * xnew[4999][3];
  if (xg_tmp + 1.0F < 1.0F - xg_tmp) {
    out = xg_tmp + 1.0F;
    b_sign = 1;
  } else {
    out = 1.0F - xg_tmp;
    b_sign = -1;
  }

  cost_tmp[0] = xnew[4999][4] - xg[4];
  cost_tmp[1] = xnew[4999][5] - xg[5];
  cost_tmp[2] = xnew[4999][6] - xg[6];

  /*  State cost Hessian */
  c_a = -b_sign * 10;
  for (i2 = 0; i2 < 3; i2++) {
    cxx[4999][i2][0] = (float)(c_a * iv[i2][0]) * xg_tmp;
    cxx[4999][i2 + 3][0] = 0.0F;
    cxx[4999][i2][1] = (float)(c_a * iv[i2][1]) * xg_tmp;
    cxx[4999][i2 + 3][1] = 0.0F;
    cxx[4999][i2][2] = (float)(c_a * iv[i2][2]) * xg_tmp;
    cxx[4999][i2 + 3][2] = 0.0F;
  }

  for (i3 = 0; i3 < 6; i3++) {
    cxx[4999][i3][3] = iv1[i3][0];
    cxx[4999][i3][4] = iv1[i3][1];
    cxx[4999][i3][5] = iv1[i3][2];
  }

  /*  State cost Jacobian */
  /*  Returns skew symmetric - cross product matrix of a 3x1 vector */
  d_a = b_sign * 10;
  b_fv[0][0] = 0.0F;
  b_fv[1][0] = -xnew[4999][3];
  b_fv[2][0] = xnew[4999][2];
  b_fv[0][1] = xnew[4999][3];
  b_fv[1][1] = 0.0F;
  b_fv[2][1] = -xnew[4999][1];
  b_fv[0][2] = -xnew[4999][2];
  b_fv[1][2] = xnew[4999][1];
  b_fv[2][2] = 0.0F;
  for (i5 = 0; i5 < 3; i5++) {
    e_a[0][i5] = (float)d_a * -xnew[4999][i5 + 1];
    e_a[i5 + 1][0] = (float)d_a * (xnew[4999][0] * (float)iv[0][i5] + b_fv[0][i5]);
    e_a[i5 + 1][1] = (float)d_a * (xnew[4999][0] * (float)iv[1][i5] + b_fv[1][i5]);
    e_a[i5 + 1][2] = (float)d_a * (xnew[4999][0] * (float)iv[2][i5] + b_fv[2][i5]);
  }

  /*  Control cost Hessian & Jacobian */
  f5 = 0.0F;
  for (i6 = 0; i6 < 3; i6++) {
    cx[4999][i6] = ((e_a[0][i6] * xg[0] + e_a[1][i6] * xg[1]) + e_a[2][i6] * xg
                    [2]) + e_a[3][i6] * xg[3];
    cx[4999][i6 + 3] = ((float)iv[0][i6] * cost_tmp[0] + (float)iv[1][i6] *
                        cost_tmp[1]) + (float)iv[2][i6] * cost_tmp[2];
    f5 += ((0.5F * cost_tmp[0] * (float)iv[i6][0] + 0.5F * cost_tmp[1] * (float)
            iv[i6][1]) + 0.5F * cost_tmp[2] * (float)iv[i6][2]) * cost_tmp[i6];
  }

  *cost += 10.0F * out + f5;
}

/*
 * Arguments    : float y[11]
 * Return Type  : void
 */
static void power(float y[11])
{
  int k;
  for (k = 0; k < 11; k++) {
    y[k] = powf(10.0F, -0.3F * (float)k);
  }
}

/*
 * Arguments    : const float x[7]
 *                const float u[3]
 *                const float B_ECI[3]
 *                float Js[6][3]
 *                float xdot[7]
 *                float dxdot[10][7]
 * Return Type  : void
 */
static void satellite_dynamics(const float x[7], const float u[3], const float
  B_ECI[3], float Js[6][3], float xdot[7], float dxdot[10][7])
{
  float fcnOutput;
  float q2_idx_1;
  float q2_idx_2;
  float q2_idx_3;
  float B_B_idx_0;
  float B_B_idx_1;
  float B_B_idx_2;
  float b_fv[3][3];
  float wdot_tmp[3][3];
  float B_B[3];
  int i;
  int i1;
  float b_fv1[3][4];
  float f;
  float b_fv2[4];
  int i2;
  float f1;
  float qdot_tmp[3][3];
  float f2;
  int i3;
  float b_wdot_tmp[3][3];
  float b_x[3];
  float c_wdot_tmp[3][3];
  int i4;
  float fv3[3][3];
  int i5;
  float b_Js[3][3];
  int i6;
  int i7;
  float dxdot_tmp;
  int i8;
  int i9;
  float c_Js[3][3];

  /*  Calculates the continuous time state derivative and Jacobians */
  /*  Inputs */
  /* ======================= */
  /*  x - the current state */
  /*  u - the current control (magentic moment vector) */
  /*  B - Earth magnetic field vector in ECI co-ordinates (3x1) */
  /*  Js - [Inertia, Inertia-inverse] */
  /*  inertia tensor */
  /*  inverse inertia tensor */
  /*  Angular velocity */
  /*  Quaternion */
  /*  Magnetic field section  */
  /*  given B_ECI (ECI magnetic field at the time step) */
  /*  Rotates a vector using a quaternion */
  /*  Multiplies quaternions */
  fcnOutput = (B_ECI[0] * x[1] + B_ECI[1] * x[2]) + B_ECI[2] * x[3];
  q2_idx_1 = x[0] * B_ECI[0] + (B_ECI[1] * x[3] - B_ECI[2] * x[2]);
  q2_idx_2 = x[0] * B_ECI[1] + (B_ECI[2] * x[1] - B_ECI[0] * x[3]);
  q2_idx_3 = x[0] * B_ECI[2] + (B_ECI[0] * x[2] - B_ECI[1] * x[1]);

  /*  Multiplies quaternions */
  B_B_idx_0 = (x[0] * q2_idx_1 + (0.0F - fcnOutput) * -x[1]) + (-x[2] * q2_idx_3
    - -x[3] * q2_idx_2);
  B_B_idx_1 = (x[0] * q2_idx_2 + (0.0F - fcnOutput) * -x[2]) + (-x[3] * q2_idx_1
    - -x[1] * q2_idx_3);
  B_B_idx_2 = (x[0] * q2_idx_3 + (0.0F - fcnOutput) * -x[3]) + (-x[1] * q2_idx_2
    - -x[2] * q2_idx_1);

  /*  Non-linear dynamics */
  /*  Returns skew symmetric - cross product matrix of a 3x1 vector */
  b_fv[0][0] = 0.0F;
  b_fv[1][0] = -x[3];
  b_fv[2][0] = x[2];
  b_fv[0][1] = x[3];
  b_fv[1][1] = 0.0F;
  b_fv[2][1] = -x[1];
  b_fv[0][2] = -x[2];
  b_fv[1][2] = x[1];
  b_fv[2][2] = 0.0F;

  /*  wdot = Jinv*(u - skew_mat(w)*J*w); */
  /*  Returns skew symmetric - cross product matrix of a 3x1 vector */
  wdot_tmp[0][0] = 0.0F;
  wdot_tmp[1][0] = -x[6];
  wdot_tmp[2][0] = x[5];
  wdot_tmp[0][1] = x[6];
  wdot_tmp[1][1] = 0.0F;
  wdot_tmp[2][1] = -x[4];
  wdot_tmp[0][2] = -x[5];
  wdot_tmp[1][2] = x[4];
  wdot_tmp[2][2] = 0.0F;

  /*  with magnetorquers */
  B_B[0] = -(B_B_idx_1 * u[2] - B_B_idx_2 * u[1]);
  B_B[1] = -(B_B_idx_2 * u[0] - B_B_idx_0 * u[2]);
  B_B[2] = -(B_B_idx_0 * u[1] - B_B_idx_1 * u[0]);
  for (i = 0; i < 3; i++) {
    b_fv1[i][0] = 0.5F * -x[i + 1];
    f = 0.0F;
    for (i2 = 0; i2 < 3; i2++) {
      f1 = x[0] * (float)iv[i][i2] + b_fv[i][i2];
      qdot_tmp[i][i2] = f1;
      f2 = (wdot_tmp[0][i] * Js[i2][0] + wdot_tmp[1][i] * Js[i2][1]) + wdot_tmp
        [2][i] * Js[i2][2];
      b_wdot_tmp[i2][i] = f2;
      b_fv1[i][i2 + 1] = 0.5F * f1;
      f += f2 * x[i2 + 4];
    }

    B_B[i] -= f;
  }

  for (i1 = 0; i1 < 4; i1++) {
    b_fv2[i1] = (b_fv1[0][i1] * x[4] + b_fv1[1][i1] * x[5]) + b_fv1[2][i1] * x[6];
  }

  xdot[0] = b_fv2[0];
  xdot[1] = b_fv2[1];
  xdot[2] = b_fv2[2];
  xdot[3] = b_fv2[3];

  /*  Jacobians  */
  for (i3 = 0; i3 < 3; i3++) {
    xdot[i3 + 4] = (Js[3][i3] * B_B[0] + Js[4][i3] * B_B[1]) + Js[5][i3] * B_B[2];
    b_x[i3] = (Js[0][i3] * x[4] + Js[1][i3] * x[5]) + Js[2][i3] * x[6];
  }

  /*  Returns skew symmetric - cross product matrix of a 3x1 vector */
  /*  B = [zeros(4,3); */
  /*       Jinv]; */
  /*  with magnetorquers: */
  /*  Returns skew symmetric - cross product matrix of a 3x1 vector */
  c_wdot_tmp[0][0] = b_wdot_tmp[0][0];
  c_wdot_tmp[1][0] = b_wdot_tmp[1][0] - (-b_x[2]);
  c_wdot_tmp[2][0] = b_wdot_tmp[2][0] - b_x[1];
  c_wdot_tmp[0][1] = b_wdot_tmp[0][1] - b_x[2];
  c_wdot_tmp[1][1] = b_wdot_tmp[1][1];
  c_wdot_tmp[2][1] = b_wdot_tmp[2][1] - (-b_x[0]);
  c_wdot_tmp[0][2] = b_wdot_tmp[0][2] - (-b_x[1]);
  c_wdot_tmp[1][2] = b_wdot_tmp[1][2] - b_x[0];
  c_wdot_tmp[2][2] = b_wdot_tmp[2][2];
  for (i4 = 0; i4 < 3; i4++) {
    for (i5 = 0; i5 < 3; i5++) {
      b_fv[i5][i4] = (-2.0F * Js[3][i4] * c_wdot_tmp[i5][0] + -2.0F * Js[4][i4] *
                      c_wdot_tmp[i5][1]) + -2.0F * Js[5][i4] * c_wdot_tmp[i5][2];
      b_Js[i4][i5] = -Js[i4 + 3][i5];
    }
  }

  fv3[0][0] = 0.0F;
  fv3[1][0] = -B_B_idx_2;
  fv3[2][0] = B_B_idx_1;
  fv3[0][1] = B_B_idx_2;
  fv3[1][1] = 0.0F;
  fv3[2][1] = -B_B_idx_0;
  fv3[0][2] = -B_B_idx_1;
  fv3[1][2] = B_B_idx_0;
  fv3[2][2] = 0.0F;
  dxdot[0][0] = 0.0F;
  for (i6 = 0; i6 < 3; i6++) {
    dxdot_tmp = x[i6 + 4];
    dxdot[i6 + 1][0] = 0.5F * -dxdot_tmp;
    dxdot[i6 + 4][0] = 0.5F * -x[i6 + 1];
    dxdot[0][i6 + 1] = 0.5F * dxdot_tmp;
    for (i9 = 0; i9 < 3; i9++) {
      c_Js[i9][i6] = (b_Js[0][i6] * fv3[i9][0] + b_Js[1][i6] * fv3[i9][1]) +
        b_Js[2][i6] * fv3[i9][2];
      dxdot[i6 + 1][i9 + 1] = 0.5F * -wdot_tmp[i6][i9];
      dxdot[i6 + 4][i9 + 1] = 0.5F * qdot_tmp[i6][i9];
    }
  }

  for (i7 = 0; i7 < 4; i7++) {
    dxdot[i7][4] = 0.0F;
    dxdot[i7][5] = 0.0F;
    dxdot[i7][6] = 0.0F;
  }

  for (i8 = 0; i8 < 3; i8++) {
    dxdot[i8 + 4][4] = 0.5F * b_fv[i8][0];
    dxdot[i8 + 4][5] = 0.5F * b_fv[i8][1];
    dxdot[i8 + 4][6] = 0.5F * b_fv[i8][2];
    dxdot[i8 + 7][0] = 0.0F;
    dxdot[i8 + 7][1] = 0.0F;
    dxdot[i8 + 7][2] = 0.0F;
    dxdot[i8 + 7][3] = 0.0F;
    dxdot[i8 + 7][4] = c_Js[i8][0];
    dxdot[i8 + 7][5] = c_Js[i8][1];
    dxdot[i8 + 7][6] = c_Js[i8][2];
  }
}

/*
 * Arguments    : float x0[5000][7]
 *                const float xg[7]
 *                float u0[4999][3]
 *                float u_lims[2][3]
 *                float B_ECI[5000][3]
 *                float Js[6][3]
 *                const float options[10]
 *                float x[5000][7]
 *                float u[4999][3]
 *                float K[4999][6][3]
 *                bool  *result
 * Return Type  : void
 */
void milqr(float x0[5000][7], const float xg[7], float u0[4999][3], float
           u_lims[2][3], float B_ECI[5000][3], float Js[6][3], const float
           options[10], float x[5000][7], float u[4999][3], float K[4999][6][3],
           bool *result)
{
  float lambda_max;
  float lambda_min;
  float lambda_scale;
  float alphas[11];
  float lambda;
  static float x_n[5000][7];
  int i;
  static float cx_n[5000][6];
  static float u_n[4999][3];
  int i1;
  static float cxx_n[5000][6][6];
  static float cu_n[4999][3];
  int i2;
  float cost_n;
  int i3;
  int i4;
  float dV[2];
  int i5;
  static float fx_n[4999][6][6];
  int i6;
  static float fu_n[4999][3][6];
  static float b_fv[4999][6][3];
  static float fx[4999][6][6];
  static float fu[4999][3][6];
  static float cx[5000][6];
  static float cu[4999][3];
  static float cxx[5000][6][6];
  static float cuu[4999][3][3];
  float cost;
  int i7;
  int i8;
  static float cuu_n[4999][3][3];
  static float l[4999][3];
  float iter;
  int b_iter;
  bool exitg1;
  bool backPassCheck;
  int exitg2;
  bool diverged;
  static float b_fv1[4999][3];
  static float b_fv2[4999][3];
  int j;
  float f;
  bool guard1 = false;
  float f1;
  float f2;
  float y;
  float maxval[4999];
  float f3;
  int k;
  float bsum;
  int b_k;
  bool lineSearchCheck;
  int alpha;
  bool exitg3;
  float b_cost_n;
  float expected_change;
  float c_ratio;
  float dcost;
  float b_u;
  int i9;
  int i10;
  int i11;
  int i12;
  int i13;
  int i14;
  int i15;
  int i16;

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
  /*  */
  /*  B_ECI - A time sequence of Earth magnetic field vectors in ECI (3, N) */
  /*  */
  /*  Js - Inertia Tensor and its inverse [J, Jinv] */
  /*  */
  /*  options - Array of solver options */
  /*  Outputs */
  /*  =========================================== */
  /*  x - Final nominal trajectory (n, N) */
  /*  */
  /*  u - Final open-loop controls (m, N-1) */
  /*  */
  /*  K - Feedback control gains (n-1, m, N-1)  */
  /*  */
  /*  result - Indicates convergence (boolean) */
  /*  */
  /* ============================================= */
  /*  Solver Options (passed in as array) */
  /*  Num simulation steps */
  /*  Timestep (Should be highest we can get away with) */
  /*  maximum iterations */
  /*  cost reduction exit tolerance */
  /*  feedforward control change exit criterion */
  /*  max regularizaion param allowed for exit */
  /*  minimum accepted cost reduction ratio */
  lambda_max = options[7];

  /*  maximum regularization parameter */
  lambda_min = options[8];

  /*  set lambda = 0 below this value */
  lambda_scale = options[9];

  /*  amount to scale dlambda by */
  /*  Initialize optimization params */
  power(alphas);

  /*  line search param vector */
  lambda = 1.0F;

  /*  error state size (3 param. error representation for attitude) */
  /*  Init matrices for update (otherwise MATLAB coder throws an error) */
  memset(&x_n[0][0], 0, 35000U * sizeof(float));
  for (i = 0; i < 4999; i++) {
    u_n[i][0] = 0.0F;
    u_n[i][1] = 0.0F;
    u_n[i][2] = 0.0F;
    for (i2 = 0; i2 < 6; i2++) {
      for (i4 = 0; i4 < 6; i4++) {
        fx_n[i][i2][i4] = 0.0F;
      }
    }

    for (i3 = 0; i3 < 3; i3++) {
      for (i5 = 0; i5 < 6; i5++) {
        fu_n[i][i3][i5] = 0.0F;
      }
    }
  }

  memset(&cx_n[0][0], 0, 30000U * sizeof(float));
  for (i1 = 0; i1 < 4999; i1++) {
    cu_n[i1][0] = 0.0F;
    cu_n[i1][1] = 0.0F;
    cu_n[i1][2] = 0.0F;
  }

  memset(&cxx_n[0][0][0], 0, 180000U * sizeof(float));
  cost_n = 0.0F;

  /*  Initial Forward rollout */
  dV[0] = 0.0F;
  dV[1] = 0.0F;
  for (i6 = 0; i6 < 4999; i6++) {
    for (i7 = 0; i7 < 3; i7++) {
      cuu_n[i6][i7][0] = 0.0F;
      cuu_n[i6][i7][1] = 0.0F;
      cuu_n[i6][i7][2] = 0.0F;
      l[i6][i7] = 0.0F;
    }

    for (i8 = 0; i8 < 6; i8++) {
      K[i6][i8][0] = 0.0F;
      K[i6][i8][1] = 0.0F;
      K[i6][i8][2] = 0.0F;
      b_fv[i6][i8][0] = 0.0F;
      b_fv[i6][i8][1] = 0.0F;
      b_fv[i6][i8][2] = 0.0F;
    }
  }

  forwardRollout(x0, xg, u0, b_fv, u_lims, options[1], B_ECI, Js, x, u, fx, fu,
                 cx, cu, cxx, cuu, &cost);

  /*  Convergence check params */
  /*  Expected cost change */
  /*  Ratio of cost change to expected cost change */
  *result = false;

  /*  Run MILQR Optimisation */
  iter = 1.0F;
  b_iter = 0;
  exitg1 = false;
  while ((!exitg1) && (b_iter <= (int)options[2] - 1)) {
    iter = (float)b_iter + 1.0F;

    /*  Backward Pass */
    /* ======================================= */
    backPassCheck = false;
    do {
      exitg2 = 0;
      if (!backPassCheck) {
        backwardPass(fx, fu, cx, cu, cxx, cuu, lambda, u_lims, u, l, K, dV,
                     &diverged);
        if (diverged) {
          /*  Warning: Cholesky factorizaton failed */
          /*  Increase regularization parameter (lambda) */
          /*  Increases or decreases the regularization parameter according */
          /*  to a non-linear scaling regime. */
          /*  increase lambda */
          lambda = fmaxf(lambda * fmaxf(lambda_scale, lambda_scale), lambda_min);
          if (lambda > lambda_max) {
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
    b_abs(l, b_fv1);
    b_abs(u, b_fv2);
    for (j = 0; j < 4999; j++) {
      f = b_fv1[j][0] / (b_fv2[j][0] + 1.0F);
      f1 = b_fv1[j][1] / (b_fv2[j][1] + 1.0F);
      f2 = b_fv1[j][2] / (b_fv2[j][2] + 1.0F);
      f3 = f;
      if (f < f1) {
        f3 = f1;
      }

      if (f3 < f2) {
        f3 = f2;
      }

      maxval[j] = f3;
    }

    /*  Avg over time of max change */
    guard1 = false;
    if (lambda < options[5]) {
      y = maxval[0];
      for (k = 0; k < 1023; k++) {
        y += maxval[k + 1];
      }

      bsum = maxval[1024];
      for (b_k = 2; b_k < 1025; b_k++) {
        bsum += maxval[b_k + 1023];
      }

      y += bsum;
      bsum = maxval[2048];
      for (b_k = 2; b_k < 1025; b_k++) {
        bsum += maxval[b_k + 2047];
      }

      y += bsum;
      bsum = maxval[3072];
      for (b_k = 2; b_k < 1025; b_k++) {
        bsum += maxval[b_k + 3071];
      }

      y += bsum;
      bsum = maxval[4096];
      for (b_k = 2; b_k < 904; b_k++) {
        bsum += maxval[b_k + 4095];
      }

      y += bsum;
      if (y / 4999.0F < options[4]) {
        /*  Success: Control change decreased below tolerance */
        *result = true;
        exitg1 = true;
      } else {
        guard1 = true;
      }
    } else {
      guard1 = true;
    }

    if (guard1) {
      /*  Forward Line-Search */
      /* =========================================== */
      lineSearchCheck = false;
      if (backPassCheck) {
        alpha = 0;
        exitg3 = false;
        while ((!exitg3) && (alpha < 11)) {
          b_forwardRollout(x, xg, u, l, K, alphas[alpha], u_lims, options[1],
                           B_ECI, Js, x_n, u_n, fx_n, fu_n, cx_n, cu_n, cxx_n,
                           cuu_n, &b_cost_n);
          cost_n = b_cost_n;
          expected_change = alphas[alpha] * dV[0] + powf(alphas[alpha], 2.0F) *
            dV[1];
          if (expected_change < 0.0F) {
            c_ratio = (b_cost_n - cost) / expected_change;
          } else {
            /*  Non-positive expected cost reduction */
            /*  actual cost change must be negative to accept the step */
            b_u = b_cost_n - cost;
            if (b_u < 0.0F) {
              c_ratio = 1.0F;
            } else if (b_u > 0.0F) {
              c_ratio = -1.0F;
            } else {
              c_ratio = -b_u;
            }
          }

          if (c_ratio > options[6]) {
            lineSearchCheck = true;
            exitg3 = true;
          } else {
            alpha++;
          }
        }
      }

      /*  Parameter Updates */
      /* ============================================= */
      if (lineSearchCheck) {
        /*  Decrease Lambda */
        /*  Increases or decreases the regularization parameter according */
        /*  to a non-linear scaling regime. */
        /*  decrease lambda */
        lambda = lambda * fminf(1.0F / lambda_scale, 1.0F / lambda_scale) *
          (float)(lambda > lambda_min);

        /*  set = 0 if lambda too small */
        dcost = cost - cost_n;

        /*  Update the trajectory and controls */
        memcpy(&x[0][0], &x_n[0][0], 35000U * sizeof(float));
        for (i9 = 0; i9 < 4999; i9++) {
          u[i9][0] = u_n[i9][0];
          u[i9][1] = u_n[i9][1];
          u[i9][2] = u_n[i9][2];
          for (i11 = 0; i11 < 6; i11++) {
            for (i14 = 0; i14 < 6; i14++) {
              fx[i9][i11][i14] = fx_n[i9][i11][i14];
            }
          }

          for (i13 = 0; i13 < 3; i13++) {
            for (i16 = 0; i16 < 6; i16++) {
              fu[i9][i13][i16] = fu_n[i9][i13][i16];
            }
          }
        }

        memcpy(&cx[0][0], &cx_n[0][0], 30000U * sizeof(float));
        for (i10 = 0; i10 < 4999; i10++) {
          cu[i10][0] = cu_n[i10][0];
          cu[i10][1] = cu_n[i10][1];
          cu[i10][2] = cu_n[i10][2];
        }

        memcpy(&cxx[0][0][0], &cxx_n[0][0][0], 180000U * sizeof(float));
        for (i12 = 0; i12 < 4999; i12++) {
          for (i15 = 0; i15 < 3; i15++) {
            cuu[i12][i15][0] = cuu_n[i12][i15][0];
            cuu[i12][i15][1] = cuu_n[i12][i15][1];
            cuu[i12][i15][2] = cuu_n[i12][i15][2];
          }
        }

        cost = cost_n;

        /*  Change in cost small enough to terminate? */
        if (dcost < options[3]) {
          *result = true;

          /*  Success: cost change < tolerance */
        }

        b_iter++;
      } else {
        /*  No cost reduction (based on cost change ratio) */
        /*  Increase lambda */
        /*  Increases or decreases the regularization parameter according */
        /*  to a non-linear scaling regime. */
        /*  increase lambda */
        lambda = fmaxf(lambda * fmaxf(lambda_scale, lambda_scale), lambda_min);
        if (lambda > lambda_max) {
          *result = false;

          /*  Diverged: new lambda > lambda_max */
          exitg1 = true;
        } else {
          b_iter++;
        }
      }
    }
  }

  if (iter == options[2]) {
    /*  Ddin't converge completely */
    *result = false;
  }
}
