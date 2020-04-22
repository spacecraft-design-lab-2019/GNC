/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: milqr_efficient.c
 *
 * MATLAB Coder version            : 4.1
 * C/C++ source code generated on  : 20-Apr-2020 22:26:32
 */

/* Include Files */
#include <math.h>
#include <string.h>
#include "milqr_efficient.h"

/* Variable Definitions */
static const signed char iv0[9] = { 1, 0, 0, 0, 1, 0, 0, 0, 1 };

static const float fv0[9] = { 0.1F, 0.0F, 0.0F, 0.0F, 0.1F, 0.0F, 0.0F, 0.0F,
  0.1F };

static const float fv1[9] = { 0.001F, 0.0F, 0.0F, 0.0F, 0.001F, 0.0F, 0.0F, 0.0F,
  0.001F };

/* Function Declarations */
static void b_abs(const float x[297], float y[297]);
static void b_forwardRollout(const float x[700], const float xg[7], const float
  u[297], const float l[297], const float K[1782], float alpha, const float
  u_lims[6], float dt, const float B_ECI[300], float xnew[700], float unew[297],
  float *cost);
static void b_sign(float *x);
static void backwardPass(float lambda, const float u_lims[6], const float u[297],
  const float x[700], const float xg[7], float dt, const float B_ECI[300], float
  l[297], float K[1782], float dV[2], boolean_T *diverged);
static void boxQPsolve(const float Quu[9], const float Qu[3], const float
  lower_lim[3], const float upper_lim[3], const float u0[3], float u[3], float
  *result, float Luu_data[], int Luu_size[2], boolean_T b_free[3]);
static void chol_free(const float A_data[], const int A_size[2], float L_data[],
                      int L_size[2], float *fail);
static void forwardRollout(const float x[700], const float xg[7], const float u
  [297], const float K[1782], const float u_lims[6], float dt, const float
  B_ECI[300], float xnew[700], float unew[297], float *cost);
static float mean(const float x[99]);
static void power(float y[11]);
static void satellite_dynamics(const float x[7], const float u[3], const float
  B_ECI[3], float xdot[7], float dxdot[70]);
static void updateLambda(float *lambda, float direction);

/* Function Definitions */

/*
 * Arguments    : const float x[297]
 *                float y[297]
 * Return Type  : void
 */
static void b_abs(const float x[297], float y[297])
{
  int k;
  for (k = 0; k < 297; k++) {
    y[k] = fabsf(x[k]);
  }
}

/*
 * Arguments    : const float x[700]
 *                const float xg[7]
 *                const float u[297]
 *                const float l[297]
 *                const float K[1782]
 *                float alpha
 *                const float u_lims[6]
 *                float dt
 *                const float B_ECI[300]
 *                float xnew[700]
 *                float unew[297]
 *                float *cost
 * Return Type  : void
 */
static void b_forwardRollout(const float x[700], const float xg[7], const float
  u[297], const float l[297], const float K[1782], float alpha, const float
  u_lims[6], float dt, const float B_ECI[300], float xnew[700], float unew[297],
  float *cost)
{
  int unew_tmp;
  float a;
  int k;
  int dx_tmp;
  float dx[6];
  int b_dx_tmp;
  int c_dx_tmp;
  float q_idx_0;
  int q_idx_1_tmp;
  int q_idx_2_tmp;
  int q_idx_3_tmp;
  float z1[3];
  float q[16];
  float f7;
  float q_error[4];
  int b_k;
  float b[7];
  float dxdot1[70];
  float b_xnew[7];
  float dxdot2[70];
  float x1[7];
  float f8;

  /*  get rid of everything after unew except cost */
  /*  Uses an rk method to roll out a trajectory */
  /*  Returns the new trajectory, cost and the derivatives along the trajectory */
  /*  If feed-forward and feed-back controls l, K are non-zero */
  /*  then the unew returned is the new control sequeunce. Otherwise unew = u */
  /*  Sizes */
  /*  Initialize outputs (Don't overwrite x or u) */
  memset(&xnew[0], 0, 700U * sizeof(float));
  memset(&unew[0], 0, 297U * sizeof(float));
  *cost = 0.0F;
  for (unew_tmp = 0; unew_tmp < 7; unew_tmp++) {
    xnew[unew_tmp] = x[unew_tmp];
  }

  a = 0.5F * dt;
  for (k = 0; k < 99; k++) {
    /*  Find the state error vector dx */
    dx_tmp = 7 * k + 4;
    dx[3] = xnew[dx_tmp] - x[dx_tmp];
    b_dx_tmp = 7 * k + 5;
    dx[4] = xnew[b_dx_tmp] - x[b_dx_tmp];
    c_dx_tmp = 7 * k + 6;
    dx[5] = xnew[c_dx_tmp] - x[c_dx_tmp];

    /*  Calculate error between qk and q_nom */
    /*  Defined as conj(q_nom)*qnew */
    /*  Returns error as Rodrigues parameters (3x1) */
    q_idx_0 = x[7 * k];
    q_idx_1_tmp = 1 + 7 * k;
    q_idx_2_tmp = 2 + 7 * k;
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
    q[4] = x[q_idx_1_tmp];
    q[8] = x[q_idx_2_tmp];
    q[12] = x[q_idx_3_tmp];
    q[1] = -x[q_idx_1_tmp];
    q[5] = q_idx_0;
    q[9] = x[q_idx_3_tmp];
    q[13] = -x[q_idx_2_tmp];
    q[2] = -x[q_idx_2_tmp];
    q[6] = -x[q_idx_3_tmp];
    q[10] = q_idx_0;
    q[14] = x[q_idx_1_tmp];
    q[3] = -x[q_idx_3_tmp];
    q[7] = x[q_idx_2_tmp];
    q[11] = -x[q_idx_1_tmp];
    q[15] = q_idx_0;
    q_idx_0 = 0.0F;
    for (unew_tmp = 0; unew_tmp < 4; unew_tmp++) {
      q_error[unew_tmp] = 0.0F;
      f7 = ((q[unew_tmp] * xnew[7 * k] + q[unew_tmp + 4] * xnew[q_idx_1_tmp]) +
            q[unew_tmp + 8] * xnew[q_idx_2_tmp]) + q[unew_tmp + 12] *
        xnew[q_idx_3_tmp];
      q_error[unew_tmp] = f7;
      q_idx_0 += f7 * f7;
    }

    q_idx_0 = sqrtf(q_idx_0);
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
      f7 = 0.0F;
      for (unew_tmp = 0; unew_tmp < 6; unew_tmp++) {
        f7 += K[(b_k + 3 * unew_tmp) + 18 * k] * dx[unew_tmp];
      }

      unew_tmp = b_k + 3 * k;
      unew[unew_tmp] = fminf(u_lims[3 + b_k], fmaxf(u_lims[b_k], (u[unew_tmp] -
        alpha * l[unew_tmp]) - f7));
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
    satellite_dynamics(*(float (*)[7])&xnew[7 * k], *(float (*)[3])&unew[3 * k],
                       *(float (*)[3])&B_ECI[3 * k], b, dxdot1);
    for (unew_tmp = 0; unew_tmp < 7; unew_tmp++) {
      b_xnew[unew_tmp] = xnew[unew_tmp + 7 * k] + a * b[unew_tmp];
    }

    satellite_dynamics(b_xnew, *(float (*)[3])&unew[3 * k], *(float (*)[3])&
                       B_ECI[3 * k], b, dxdot2);
    for (unew_tmp = 0; unew_tmp < 7; unew_tmp++) {
      x1[unew_tmp] = xnew[unew_tmp + 7 * k] + dt * b[unew_tmp];
    }

    /*  Re-normalize the quaternion */
    q_idx_0 = sqrtf(((x1[0] * x1[0] + x1[1] * x1[1]) + x1[2] * x1[2]) + x1[3] *
                    x1[3]);
    unew_tmp = 7 * (k + 1);
    xnew[unew_tmp] = x1[0] / q_idx_0;
    xnew[1 + unew_tmp] = x1[1] / q_idx_0;
    xnew[2 + unew_tmp] = x1[2] / q_idx_0;
    xnew[3 + unew_tmp] = x1[3] / q_idx_0;
    xnew[unew_tmp + 4] = x1[4];
    xnew[unew_tmp + 5] = x1[5];
    xnew[unew_tmp + 6] = x1[6];

    /*  Continuous time Jacobians */
    /*  Discrete time Jacobians */
    /*  First form the state error Jacobians E(x_k) & E(x_k+1) */
    /*  Calculate the cost */
    /*  make second function that's satellite cost derivatives */
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
    /*  cumulative angular velocity cost hessian */
    /*  cumulative geodesic cost weighting */
    /*  Finds the geodesic quaternion-error cost */
    /*  quat_cost = min(1+qd'q, 1-qd'q)) where qd is the desired attitude quaternion  */
    /*  Also records the sign (+ or -) which minimizes quat_cost */
    /*  this is used when calculating the Jacobain and Hessian */
    q_idx_0 = ((xg[0] * xnew[7 * k] + xg[1] * xnew[q_idx_1_tmp]) + xg[2] *
               xnew[q_idx_2_tmp]) + xg[3] * xnew[q_idx_3_tmp];
    z1[0] = xnew[dx_tmp] - xg[4];
    z1[1] = xnew[b_dx_tmp] - xg[5];
    z1[2] = xnew[c_dx_tmp] - xg[6];

    /*  State cost Hessian */
    /*  State cost Jacobian */
    /*  Control cost Hessian & Jacobian */
    f7 = 0.0F;
    for (unew_tmp = 0; unew_tmp < 3; unew_tmp++) {
      f7 += ((0.5F * z1[0] * fv0[3 * unew_tmp] + 0.5F * z1[1] * fv0[1 + 3 *
              unew_tmp]) + 0.5F * z1[2] * fv0[2 + 3 * unew_tmp]) * z1[unew_tmp];
    }

    f8 = 0.0F;
    for (unew_tmp = 0; unew_tmp < 3; unew_tmp++) {
      f8 += ((0.5F * unew[3 * k] * fv1[3 * unew_tmp] + 0.5F * unew[1 + 3 * k] *
              fv1[1 + 3 * unew_tmp]) + 0.5F * unew[2 + 3 * k] * fv1[2 + 3 *
             unew_tmp]) * unew[unew_tmp + 3 * k];
    }

    *cost += (fminf(1.0F + q_idx_0, 1.0F - q_idx_0) + f7) + f8;
  }

  /*  Final cost */
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
  q_idx_0 = ((xg[0] * xnew[693] + xg[1] * xnew[694]) + xg[2] * xnew[695]) + xg[3]
    * xnew[696];
  z1[0] = xnew[697] - xg[4];
  z1[1] = xnew[698] - xg[5];
  z1[2] = xnew[699] - xg[6];

  /*  State cost Hessian */
  /*  State cost Jacobian */
  /*  Control cost Hessian & Jacobian */
  f7 = 0.0F;
  for (unew_tmp = 0; unew_tmp < 3; unew_tmp++) {
    f7 += ((0.5F * z1[0] * (float)iv0[3 * unew_tmp] + 0.5F * z1[1] * (float)iv0
            [1 + 3 * unew_tmp]) + 0.5F * z1[2] * (float)iv0[2 + 3 * unew_tmp]) *
      z1[unew_tmp];
  }

  *cost += 10.0F * fminf(1.0F + q_idx_0, 1.0F - q_idx_0) + f7;
}

/*
 * Arguments    : float *x
 * Return Type  : void
 */
static void b_sign(float *x)
{
  if (*x < 0.0F) {
    *x = -1.0F;
  } else {
    if (*x > 0.0F) {
      *x = 1.0F;
    }
  }
}

/*
 * Arguments    : float lambda
 *                const float u_lims[6]
 *                const float u[297]
 *                const float x[700]
 *                const float xg[7]
 *                float dt
 *                const float B_ECI[300]
 *                float l[297]
 *                float K[1782]
 *                float dV[2]
 *                boolean_T *diverged
 * Return Type  : void
 */
static void backwardPass(float lambda, const float u_lims[6], const float u[297],
  const float x[700], const float xg[7], float dt, const float B_ECI[300], float
  l[297], float K[1782], float dV[2], boolean_T *diverged)
{
  float xg_tmp;
  float b_xg;
  int c_sign;
  int a;
  int i2;
  float Vxx[36];
  int k;
  static const signed char iv1[18] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1,
    0, 0, 0, 1 };

  float Quu[9];
  float b_a[12];
  int a_tmp;
  float Vx[6];
  float b_x[3];
  int b_k;
  boolean_T exitg1;
  float b[7];
  float dxdot1[70];
  float c_x[7];
  float dxdot2[70];
  float x1[7];
  int i3;
  float fx_tmp[49];
  float b_fx_tmp[42];
  float f3;
  int i4;
  float f4;
  float fv6[49];
  static const signed char iv2[49] = { 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
    0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
    1, 0, 0, 0, 0, 0, 0, 0, 1 };

  float f5;
  float d_x[42];
  float c_fx_tmp[42];
  float f6;
  float fx[36];
  float fu[18];
  float b_dt[21];
  int c_a;
  float Qu_tmp[18];
  float Qx_tmp[36];
  float b_Quu[9];
  float Qu[3];
  float b_u_lims[3];
  float Quu_tmp[18];
  float fv7[3];
  float Qux[18];
  float lk[3];
  float result;
  float Luu_data[9];
  int Luu_size[2];
  boolean_T b_free[3];
  float Kk[18];
  boolean_T out;
  boolean_T exitg2;
  signed char tmp_data[3];
  float b_Kk[18];
  float d_a[6];
  float e_a[6];
  float f_a[36];
  float b_Qx_tmp[36];
  int value;
  int n;
  int Y_size_idx_0;
  float Y_data[18];
  float coeffs_data[3];
  int a_size_idx_1;
  int i;
  int loop_ub;
  int b_loop_ub;
  int X_size_idx_0;
  float a_data[3];
  int j;
  int LT_size_idx_0;
  float b_data[3];
  float LT_data[9];
  int c_loop_ub;
  int d_loop_ub;
  int b_i;
  signed char b_tmp_data[3];

  /*  Perfoms the LQR backward pass to find the optimal controls */
  /*  Solves a quadratic program (QP) at each timestep for the optimal */
  /*  controls given the control limits */
  /*  insert section where I evaluate the cost and dynamics derivatives */
  /*  function [cx,cu,cxx,cuu,cxu] = cost_derivatives(x,u) */
  /*  function [fx,fu] = state_derivatives(x,u) */
  /*  Initialize matrices (for C code, not needed in MATLAB) */
  memset(&l[0], 0, 297U * sizeof(float));
  memset(&K[0], 0, 1782U * sizeof(float));

  /*  Change in cost */
  dV[0] = 0.0F;
  dV[1] = 0.0F;

  /*  Set cost-to-go Jacobian and Hessian equal to final costs */
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
  xg_tmp = xg[0] * x[693];
  b_xg = ((xg_tmp + xg[1] * x[694]) + xg[2] * x[695]) + xg[3] * x[696];
  if (1.0F + b_xg < 1.0F - b_xg) {
    c_sign = 1;
  } else {
    c_sign = -1;
  }

  /*  State cost Hessian */
  a = -c_sign * 10;
  b_xg = ((xg_tmp + xg[1] * x[694]) + xg[2] * x[695]) + xg[3] * x[696];
  for (i2 = 0; i2 < 3; i2++) {
    Vxx[6 * i2] = (float)(a * iv0[3 * i2]) * b_xg;
    k = 6 * (i2 + 3);
    Vxx[k] = 0.0F;
    Vxx[1 + 6 * i2] = (float)(a * iv0[1 + 3 * i2]) * b_xg;
    Vxx[1 + k] = 0.0F;
    Vxx[2 + 6 * i2] = (float)(a * iv0[2 + 3 * i2]) * b_xg;
    Vxx[2 + k] = 0.0F;
  }

  for (i2 = 0; i2 < 6; i2++) {
    Vxx[6 * i2 + 3] = iv1[3 * i2];
    Vxx[6 * i2 + 4] = iv1[1 + 3 * i2];
    Vxx[6 * i2 + 5] = iv1[2 + 3 * i2];
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
  Quu[0] = 0.0F;
  Quu[3] = -x[696];
  Quu[6] = x[695];
  Quu[1] = x[696];
  Quu[4] = 0.0F;
  Quu[7] = -x[694];
  Quu[2] = -x[695];
  Quu[5] = x[694];
  Quu[8] = 0.0F;
  for (i2 = 0; i2 < 3; i2++) {
    b_a[i2] = (float)a * -x[i2 + 694];
    a_tmp = 3 * (i2 + 1);
    b_a[a_tmp] = (float)a * (x[693] * (float)iv0[i2] + Quu[i2]);
    b_a[1 + a_tmp] = (float)a * (x[693] * (float)iv0[i2 + 3] + Quu[i2 + 3]);
    b_a[2 + a_tmp] = (float)a * (x[693] * (float)iv0[i2 + 6] + Quu[i2 + 6]);
    b_x[i2] = x[i2 + 697] - xg[4 + i2];
  }

  for (i2 = 0; i2 < 3; i2++) {
    Vx[i2] = ((b_a[i2] * xg[0] + b_a[i2 + 3] * xg[1]) + b_a[i2 + 6] * xg[2]) +
      b_a[i2 + 9] * xg[3];
    Vx[i2 + 3] = ((float)iv0[i2] * b_x[0] + (float)iv0[i2 + 3] * b_x[1]) +
      (float)iv0[i2 + 6] * b_x[2];
  }

  /*  Control cost Hessian & Jacobian */
  *diverged = false;
  b_k = 0;
  exitg1 = false;
  while ((!exitg1) && (b_k < 99)) {
    /*  calculate dynamics derivatives */
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
    satellite_dynamics(*(float (*)[7])&x[7 * (98 - b_k)], *(float (*)[3])&u[3 *
                       (98 - b_k)], *(float (*)[3])&B_ECI[3 * (98 - b_k)], b,
                       dxdot1);
    xg_tmp = 0.5F * dt;
    for (i2 = 0; i2 < 7; i2++) {
      c_x[i2] = x[i2 + 7 * (98 - b_k)] + xg_tmp * b[i2];
    }

    satellite_dynamics(c_x, *(float (*)[3])&u[3 * (98 - b_k)], *(float (*)[3])&
                       B_ECI[3 * (98 - b_k)], b, dxdot2);

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
    xg_tmp = 0.5F * dt * dt;
    for (i2 = 0; i2 < 7; i2++) {
      x1[i2] = x[i2 + 7 * (98 - b_k)] + dt * b[i2];
      for (i3 = 0; i3 < 7; i3++) {
        a_tmp = i3 + 7 * i2;
        fx_tmp[a_tmp] = xg_tmp * dxdot2[a_tmp];
      }
    }

    Quu[0] = 0.0F;
    Quu[3] = -x1[3];
    Quu[6] = x1[2];
    Quu[1] = x1[3];
    Quu[4] = 0.0F;
    Quu[7] = -x1[1];
    Quu[2] = -x1[2];
    Quu[5] = x1[1];
    Quu[8] = 0.0F;
    for (i2 = 0; i2 < 3; i2++) {
      b_fx_tmp[i2] = -x1[1 + i2];
      b_fx_tmp[i2 + 3] = 0.0F;
      a_tmp = 6 * (i2 + 1);
      b_fx_tmp[a_tmp] = x1[0] * (float)iv0[i2] + Quu[i2];
      b_fx_tmp[a_tmp + 3] = 0.0F;
      b_fx_tmp[1 + a_tmp] = x1[0] * (float)iv0[i2 + 3] + Quu[i2 + 3];
      b_fx_tmp[a_tmp + 4] = 0.0F;
      b_fx_tmp[2 + a_tmp] = x1[0] * (float)iv0[i2 + 6] + Quu[i2 + 6];
      b_fx_tmp[a_tmp + 5] = 0.0F;
      for (i3 = 0; i3 < 6; i3++) {
        b_fx_tmp[i3 + 6 * (i2 + 4)] = iv1[i2 + 3 * i3];
      }
    }

    for (i2 = 0; i2 < 7; i2++) {
      for (i3 = 0; i3 < 7; i3++) {
        f3 = 0.0F;
        for (i4 = 0; i4 < 7; i4++) {
          f3 += fx_tmp[i2 + 7 * i4] * dxdot1[i4 + 7 * i3];
        }

        i4 = i2 + 7 * i3;
        fv6[i4] = ((float)iv2[i4] + dt * dxdot2[i4]) + f3;
      }
    }

    xg_tmp = x[7 * (98 - b_k)];
    Quu[0] = 0.0F;
    f3 = x[3 + 7 * (98 - b_k)];
    Quu[3] = -f3;
    f4 = x[2 + 7 * (98 - b_k)];
    Quu[6] = f4;
    Quu[1] = f3;
    Quu[4] = 0.0F;
    f5 = x[1 + 7 * (98 - b_k)];
    Quu[7] = -f5;
    Quu[2] = -f4;
    Quu[5] = f5;
    Quu[8] = 0.0F;
    for (i2 = 0; i2 < 6; i2++) {
      for (i3 = 0; i3 < 7; i3++) {
        a_tmp = i2 + 6 * i3;
        c_fx_tmp[a_tmp] = 0.0F;
        f6 = 0.0F;
        for (i4 = 0; i4 < 7; i4++) {
          f6 += b_fx_tmp[i2 + 6 * i4] * fv6[i4 + 7 * i3];
        }

        c_fx_tmp[a_tmp] = f6;
      }
    }

    for (i2 = 0; i2 < 3; i2++) {
      d_x[7 * i2] = -x[(i2 + 7 * (98 - b_k)) + 1];
      k = 7 * (i2 + 3);
      d_x[k] = 0.0F;
      d_x[7 * i2 + 1] = xg_tmp * (float)iv0[3 * i2] + Quu[3 * i2];
      d_x[k + 1] = 0.0F;
      a_tmp = 1 + 3 * i2;
      d_x[7 * i2 + 2] = xg_tmp * (float)iv0[a_tmp] + Quu[a_tmp];
      d_x[k + 2] = 0.0F;
      a_tmp = 2 + 3 * i2;
      d_x[7 * i2 + 3] = xg_tmp * (float)iv0[a_tmp] + Quu[a_tmp];
      d_x[k + 3] = 0.0F;
    }

    for (i2 = 0; i2 < 6; i2++) {
      d_x[7 * i2 + 4] = iv1[3 * i2];
      d_x[7 * i2 + 5] = iv1[1 + 3 * i2];
      d_x[7 * i2 + 6] = iv1[2 + 3 * i2];
    }

    for (i2 = 0; i2 < 6; i2++) {
      for (i3 = 0; i3 < 6; i3++) {
        k = i2 + 6 * i3;
        fx[k] = 0.0F;
        f6 = 0.0F;
        for (i4 = 0; i4 < 7; i4++) {
          f6 += c_fx_tmp[i2 + 6 * i4] * d_x[i4 + 7 * i3];
        }

        fx[k] = f6;
      }
    }

    for (i2 = 0; i2 < 7; i2++) {
      for (i3 = 0; i3 < 3; i3++) {
        f6 = 0.0F;
        for (i4 = 0; i4 < 7; i4++) {
          f6 += fx_tmp[i2 + 7 * i4] * dxdot1[i4 + 7 * (7 + i3)];
        }

        b_dt[i2 + 7 * i3] = dt * dxdot2[i2 + 7 * (7 + i3)] + f6;
      }
    }

    for (i2 = 0; i2 < 6; i2++) {
      for (i3 = 0; i3 < 3; i3++) {
        k = i2 + 6 * i3;
        fu[k] = 0.0F;
        f6 = 0.0F;
        for (i4 = 0; i4 < 7; i4++) {
          f6 += b_fx_tmp[i2 + 6 * i4] * b_dt[i4 + 7 * i3];
        }

        fu[k] = f6;
      }
    }

    /*  Define cost gradients and hessians */
    /*  convert cost and dynamics derivatives into functions */
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
    xg_tmp = xg[0] * x[7 * (98 - b_k)];
    b_xg = ((xg_tmp + xg[1] * f5) + xg[2] * f4) + xg[3] * f3;
    if (1.0F + b_xg < 1.0F - b_xg) {
      c_sign = 1;
    } else {
      c_sign = -1;
    }

    /*  State cost Hessian */
    a = -c_sign * 10;
    b_xg = ((xg_tmp + xg[1] * x[1 + 7 * (98 - b_k)]) + xg[2] * x[2 + 7 * (98 -
             b_k)]) + xg[3] * x[3 + 7 * (98 - b_k)];

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
    c_a = c_sign * 10;

    /*  Control cost Hessian & Jacobian */
    for (i2 = 0; i2 < 6; i2++) {
      for (i3 = 0; i3 < 6; i3++) {
        Qx_tmp[i3 + 6 * i2] = fx[i2 + 6 * i3];
      }

      Qu_tmp[3 * i2] = fu[i2];
      Qu_tmp[1 + 3 * i2] = fu[i2 + 6];
      Qu_tmp[2 + 3 * i2] = fu[i2 + 12];
    }

    for (i2 = 0; i2 < 3; i2++) {
      a_tmp = 3 * (98 - b_k);
      b_x[i2] = 0.0F;
      for (i3 = 0; i3 < 6; i3++) {
        k = i2 + 3 * i3;
        Quu_tmp[k] = 0.0F;
        f6 = 0.0F;
        for (i4 = 0; i4 < 6; i4++) {
          f6 += Qu_tmp[i2 + 3 * i4] * Vxx[i4 + 6 * i3];
        }

        Quu_tmp[k] = f6;
        b_x[i2] += Qu_tmp[k] * Vx[i3];
      }

      Qu[i2] = ((fv1[i2] * u[a_tmp] + fv1[i2 + 3] * u[1 + a_tmp]) + fv1[i2 + 6] *
                u[2 + a_tmp]) + b_x[i2];
      for (i3 = 0; i3 < 3; i3++) {
        f6 = 0.0F;
        for (i4 = 0; i4 < 6; i4++) {
          f6 += Quu_tmp[i2 + 3 * i4] * fu[i4 + 6 * i3];
        }

        k = i2 + 3 * i3;
        b_Quu[k] = fv1[k] + f6;
      }

      for (i3 = 0; i3 < 6; i3++) {
        k = i2 + 3 * i3;
        Qux[k] = 0.0F;
        f6 = 0.0F;
        for (i4 = 0; i4 < 6; i4++) {
          f6 += Quu_tmp[i2 + 3 * i4] * fx[i4 + 6 * i3];
        }

        Qux[k] = f6;
      }
    }

    /*  Regularization (for Cholesky positive definiteness) */
    /*  Solve the Quadratic program with control limits */
    for (i2 = 0; i2 < 9; i2++) {
      Quu[i2] = b_Quu[i2] + (float)iv0[i2] * lambda;
    }

    b_x[0] = u_lims[0] - u[3 * (98 - b_k)];
    b_u_lims[0] = u_lims[3] - u[3 * (98 - b_k)];
    i2 = 3 * ((int)fminf(99.0F, (99.0F + -(float)b_k) + 1.0F) - 1);
    fv7[0] = -l[i2];
    b_x[1] = u_lims[1] - u[1 + 3 * (98 - b_k)];
    b_u_lims[1] = u_lims[4] - u[1 + 3 * (98 - b_k)];
    fv7[1] = -l[1 + i2];
    b_x[2] = u_lims[2] - u[2 + 3 * (98 - b_k)];
    b_u_lims[2] = u_lims[5] - u[2 + 3 * (98 - b_k)];
    fv7[2] = -l[2 + i2];
    boxQPsolve(Quu, Qu, b_x, b_u_lims, fv7, lk, &result, Luu_data, Luu_size,
               b_free);
    if (result < 2.0F) {
      *diverged = true;

      /*  fprintf('\nDiverged with lambda = %f\n',lambda); */
      exitg1 = true;
    } else {
      /*  Solve for feedback gains in non-clamped rows of u */
      /*  (using cholesky factor of Quu) */
      memset(&Kk[0], 0, 18U * sizeof(float));
      out = false;
      k = 0;
      exitg2 = false;
      while ((!exitg2) && (k < 3)) {
        if (b_free[k]) {
          out = true;
          exitg2 = true;
        } else {
          k++;
        }
      }

      if (out) {
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
        value = Luu_size[0];
        n = Luu_size[0];

        /*  Check sizes match and L is lower-triangular */
        /*  Solve LY = B for Y */
        /* ======================= */
        Y_size_idx_0 = Luu_size[0];
        a_tmp = Luu_size[0] * 6;
        if (0 <= a_tmp - 1) {
          memset(&Y_data[0], 0, (unsigned int)(a_tmp * (int)sizeof(float)));
        }

        if (0 <= Luu_size[0] - 1) {
          memset(&coeffs_data[0], 0, (unsigned int)(Luu_size[0] * (int)sizeof
                  (float)));
        }

        i2 = Luu_size[0];
        if (0 <= Luu_size[0] - 1) {
          a_size_idx_1 = Luu_size[0];
          loop_ub = Luu_size[0];
          b_loop_ub = Luu_size[0];
        }

        for (i = 0; i < i2; i++) {
          if (1 + i != 1) {
            for (i3 = 0; i3 < i; i3++) {
              coeffs_data[i3] = Luu_data[i + Luu_size[0] * i3];
            }
          }

          if (0 <= loop_ub - 1) {
            memcpy(&a_data[0], &coeffs_data[0], (unsigned int)(loop_ub * (int)
                    sizeof(float)));
          }

          for (j = 0; j < 6; j++) {
            for (i3 = 0; i3 < b_loop_ub; i3++) {
              b_data[i3] = Y_data[i3 + Y_size_idx_0 * j];
            }

            if ((a_size_idx_1 == 1) || (Y_size_idx_0 == 1)) {
              xg_tmp = 0.0F;
              for (i3 = 0; i3 < a_size_idx_1; i3++) {
                xg_tmp += a_data[i3] * b_data[i3];
              }
            } else {
              xg_tmp = 0.0F;
              for (i3 = 0; i3 < a_size_idx_1; i3++) {
                xg_tmp += a_data[i3] * b_data[i3];
              }
            }

            Y_data[i + Y_size_idx_0 * j] = (Qux[(tmp_data[i] + 3 * j) - 1] -
              xg_tmp) / Luu_data[i + Luu_size[0] * i];
          }
        }

        /*  Solve L'X = Y for X */
        /* ======================= */
        X_size_idx_0 = Luu_size[0];
        a_tmp = Luu_size[0] * 6;
        if (0 <= a_tmp - 1) {
          memset(&Qu_tmp[0], 0, (unsigned int)(a_tmp * (int)sizeof(float)));
        }

        LT_size_idx_0 = Luu_size[1];
        a_tmp = Luu_size[0];
        for (i2 = 0; i2 < a_tmp; i2++) {
          k = Luu_size[1];
          for (i3 = 0; i3 < k; i3++) {
            LT_data[i3 + LT_size_idx_0 * i2] = Luu_data[i2 + Luu_size[0] * i3];
          }
        }

        if (0 <= Luu_size[0] - 1) {
          memset(&coeffs_data[0], 0, (unsigned int)(Luu_size[0] * (int)sizeof
                  (float)));
        }

        i2 = (int)((1.0F + (-1.0F - (float)Luu_size[0])) / -1.0F);
        if (0 <= i2 - 1) {
          a_size_idx_1 = Luu_size[0];
          c_loop_ub = Luu_size[0];
          d_loop_ub = Luu_size[0];
        }

        for (i = 0; i < i2; i++) {
          b_i = (value - i) - 1;
          if (b_i + 1 != value) {
            if (b_i + 2 > n) {
              i3 = -1;
              i4 = 0;
              k = -1;
            } else {
              i3 = b_i;
              i4 = value;
              k = b_i;
            }

            a_tmp = i4 - i3;
            for (i4 = 0; i4 <= a_tmp - 2; i4++) {
              coeffs_data[(k + i4) + 1] = LT_data[b_i + LT_size_idx_0 * ((i3 +
                i4) + 1)];
            }
          }

          if (0 <= c_loop_ub - 1) {
            memcpy(&a_data[0], &coeffs_data[0], (unsigned int)(c_loop_ub * (int)
                    sizeof(float)));
          }

          for (j = 0; j < 6; j++) {
            for (i3 = 0; i3 < d_loop_ub; i3++) {
              b_data[i3] = Qu_tmp[i3 + X_size_idx_0 * j];
            }

            if ((a_size_idx_1 == 1) || (X_size_idx_0 == 1)) {
              xg_tmp = 0.0F;
              for (i3 = 0; i3 < a_size_idx_1; i3++) {
                xg_tmp += a_data[i3] * b_data[i3];
              }
            } else {
              xg_tmp = 0.0F;
              for (i3 = 0; i3 < a_size_idx_1; i3++) {
                xg_tmp += a_data[i3] * b_data[i3];
              }
            }

            Qu_tmp[b_i + X_size_idx_0 * j] = (Y_data[b_i + Y_size_idx_0 * j] -
              xg_tmp) / LT_data[b_i + LT_size_idx_0 * b_i];
          }
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

        a_tmp = Luu_size[0];
        for (i2 = 0; i2 < 6; i2++) {
          for (i3 = 0; i3 < a_tmp; i3++) {
            Kk[(b_tmp_data[i3] + 3 * i2) - 1] = -Qu_tmp[i3 + X_size_idx_0 * i2];
          }
        }
      }

      /*  Update Cost to Go Jacobian and Hessian */
      for (i2 = 0; i2 < 3; i2++) {
        for (i3 = 0; i3 < 6; i3++) {
          b_Kk[i3 + 6 * i2] = Kk[i2 + 3 * i3];
        }
      }

      memcpy(&Qu_tmp[0], &b_Kk[0], 18U * sizeof(float));
      for (i2 = 0; i2 < 6; i2++) {
        for (i3 = 0; i3 < 3; i3++) {
          k = i2 + 6 * i3;
          b_Kk[k] = 0.0F;
          b_Kk[k] = (Qu_tmp[i2] * b_Quu[3 * i3] + Qu_tmp[i2 + 6] * b_Quu[1 + 3 *
                     i3]) + Qu_tmp[i2 + 12] * b_Quu[2 + 3 * i3];
        }
      }

      memcpy(&Quu_tmp[0], &b_Kk[0], 18U * sizeof(float));
      xg_tmp = x[7 * (98 - b_k)];
      Quu[0] = 0.0F;
      Quu[3] = -x[3 + 7 * (98 - b_k)];
      Quu[6] = f4;
      Quu[1] = f3;
      Quu[4] = 0.0F;
      Quu[7] = -x[1 + 7 * (98 - b_k)];
      Quu[2] = -x[2 + 7 * (98 - b_k)];
      Quu[5] = f5;
      Quu[8] = 0.0F;
      for (i2 = 0; i2 < 3; i2++) {
        for (i3 = 0; i3 < 6; i3++) {
          b_Kk[i3 + 6 * i2] = Qux[i2 + 3 * i3];
        }

        a_tmp = i2 + 7 * (98 - b_k);
        b_a[i2] = (float)c_a * -x[a_tmp + 1];
        k = 3 * (i2 + 1);
        b_a[k] = (float)c_a * (xg_tmp * (float)iv0[i2] + Quu[i2]);
        b_a[1 + k] = (float)c_a * (xg_tmp * (float)iv0[i2 + 3] + Quu[i2 + 3]);
        b_a[2 + k] = (float)c_a * (xg_tmp * (float)iv0[i2 + 6] + Quu[i2 + 6]);
        b_x[i2] = x[a_tmp + 4] - xg[4 + i2];
      }

      for (i2 = 0; i2 < 3; i2++) {
        d_a[i2] = ((b_a[i2] * xg[0] + b_a[i2 + 3] * xg[1]) + b_a[i2 + 6] * xg[2])
          + b_a[i2 + 9] * xg[3];
        d_a[i2 + 3] = ((float)iv0[i2] * b_x[0] + (float)iv0[i2 + 3] * b_x[1]) +
          (float)iv0[i2 + 6] * b_x[2];
      }

      for (i2 = 0; i2 < 6; i2++) {
        xg_tmp = 0.0F;
        for (i3 = 0; i3 < 6; i3++) {
          xg_tmp += Qx_tmp[i2 + 6 * i3] * Vx[i3];
        }

        e_a[i2] = (d_a[i2] + xg_tmp) + ((Quu_tmp[i2] * lk[0] + Quu_tmp[i2 + 6] *
          lk[1]) + Quu_tmp[i2 + 12] * lk[2]);
      }

      for (i2 = 0; i2 < 6; i2++) {
        Vx[i2] = (e_a[i2] + ((Qu_tmp[i2] * Qu[0] + Qu_tmp[i2 + 6] * Qu[1]) +
                             Qu_tmp[i2 + 12] * Qu[2])) + ((b_Kk[i2] * lk[0] +
          b_Kk[i2 + 6] * lk[1]) + b_Kk[i2 + 12] * lk[2]);
        for (i3 = 0; i3 < 6; i3++) {
          k = i2 + 6 * i3;
          b_Qx_tmp[k] = 0.0F;
          f3 = 0.0F;
          for (i4 = 0; i4 < 6; i4++) {
            f3 += Qx_tmp[i2 + 6 * i4] * Vxx[i4 + 6 * i3];
          }

          b_Qx_tmp[k] = f3;
        }
      }

      for (i2 = 0; i2 < 3; i2++) {
        f_a[6 * i2] = (float)(a * iv0[3 * i2]) * b_xg;
        a_tmp = 6 * (i2 + 3);
        f_a[a_tmp] = 0.0F;
        f_a[1 + 6 * i2] = (float)(a * iv0[1 + 3 * i2]) * b_xg;
        f_a[1 + a_tmp] = 0.0F;
        f_a[2 + 6 * i2] = (float)(a * iv0[2 + 3 * i2]) * b_xg;
        f_a[2 + a_tmp] = 0.0F;
      }

      for (i2 = 0; i2 < 6; i2++) {
        f_a[6 * i2 + 3] = iv1[3 * i2];
        f_a[6 * i2 + 4] = iv1[1 + 3 * i2];
        f_a[6 * i2 + 5] = iv1[2 + 3 * i2];
        for (i3 = 0; i3 < 6; i3++) {
          k = i2 + 6 * i3;
          Qx_tmp[k] = 0.0F;
          f3 = 0.0F;
          for (i4 = 0; i4 < 6; i4++) {
            f3 += b_Qx_tmp[i2 + 6 * i4] * fx[i4 + 6 * i3];
          }

          Qx_tmp[k] = f3;
        }
      }

      for (i2 = 0; i2 < 6; i2++) {
        for (i3 = 0; i3 < 6; i3++) {
          f3 = Kk[3 * i3];
          i4 = 1 + 3 * i3;
          k = 2 + 3 * i3;
          a_tmp = i2 + 6 * i3;
          b_Qx_tmp[a_tmp] = (f_a[a_tmp] + Qx_tmp[a_tmp]) + ((Quu_tmp[i2] * f3 +
            Quu_tmp[i2 + 6] * Kk[i4]) + Quu_tmp[i2 + 12] * Kk[k]);
          Vxx[a_tmp] = (b_Qx_tmp[a_tmp] + ((Qu_tmp[i2] * Qux[3 * i3] + Qu_tmp[i2
            + 6] * Qux[i4]) + Qu_tmp[i2 + 12] * Qux[k])) + ((b_Kk[i2] * f3 +
            b_Kk[i2 + 6] * Kk[1 + 3 * i3]) + b_Kk[i2 + 12] * Kk[2 + 3 * i3]);
        }
      }

      for (i2 = 0; i2 < 6; i2++) {
        for (i3 = 0; i3 < 6; i3++) {
          k = i3 + 6 * i2;
          b_Qx_tmp[k] = 0.5F * (Vxx[k] + Vxx[i2 + 6 * i3]);
        }
      }

      memcpy(&Vxx[0], &b_Qx_tmp[0], 36U * sizeof(float));

      /*  Ensure Hessian is symmetric */
      /*  Record control cost change to check convergence */
      xg_tmp = 0.0F;
      for (i2 = 0; i2 < 3; i2++) {
        xg_tmp += ((0.5F * lk[0] * b_Quu[3 * i2] + 0.5F * lk[1] * b_Quu[1 + 3 *
                    i2]) + 0.5F * lk[2] * b_Quu[2 + 3 * i2]) * lk[i2];
      }

      dV[0] += (lk[0] * Qu[0] + lk[1] * Qu[1]) + lk[2] * Qu[2];
      dV[1] += xg_tmp;

      /*  Update Control Vectors */
      l[3 * (98 - b_k)] = -lk[0];
      l[1 + 3 * (98 - b_k)] = -lk[1];
      l[2 + 3 * (98 - b_k)] = -lk[2];
      for (i2 = 0; i2 < 6; i2++) {
        k = 3 * i2 + 18 * (98 - b_k);
        K[k] = -Kk[3 * i2];
        K[k + 1] = -Kk[1 + 3 * i2];
        K[k + 2] = -Kk[2 + 3 * i2];
      }

      b_k++;
    }
  }
}

/*
 * Arguments    : const float Quu[9]
 *                const float Qu[3]
 *                const float lower_lim[3]
 *                const float upper_lim[3]
 *                const float u0[3]
 *                float u[3]
 *                float *result
 *                float Luu_data[]
 *                int Luu_size[2]
 *                boolean_T b_free[3]
 * Return Type  : void
 */
static void boxQPsolve(const float Quu[9], const float Qu[3], const float
  lower_lim[3], const float upper_lim[3], const float u0[3], float u[3], float
  *result, float Luu_data[], int Luu_size[2], boolean_T b_free[3])
{
  boolean_T clamped[3];
  int i5;
  float old_cost;
  float scale;
  float z1[3];
  float absxk;
  float t;
  float cost;
  int iter;
  int b_iter;
  boolean_T exitg1;
  int i;
  boolean_T out;
  float grad[3];
  int k;
  boolean_T prev_clamped[3];
  boolean_T exitg2;
  boolean_T b_clamped;
  boolean_T factorize;
  boolean_T guard1 = false;
  int trueCount;
  signed char tmp_data[3];
  signed char b_tmp_data[3];
  int Quu_size[2];
  float Quu_data[9];
  float indef;
  int i6;
  float y;
  float grad_norm;
  float grad_clamped[3];
  signed char c_tmp_data[3];
  int value;
  int n;
  float Y_data[3];
  float coeffs_data[3];
  int a_size_idx_1;
  int loop_ub;
  int b_loop_ub;
  float a_data[3];
  int LT_size_idx_0;
  int c_loop_ub;
  float b_data[3];
  float LT_data[9];
  int d_loop_ub;
  int e_loop_ub;
  int b_i;
  float delta_u;
  float b_delta_u[3];
  float step;
  float u_c_idx_0;
  float u_c_idx_1;
  float u_c_idx_2;
  float cost_c;

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
  for (i5 = 0; i5 < 9; i5++) {
    Luu_data[i5] = 0.0F;
  }

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
  scale = fmaxf(lower_lim[0], fminf(upper_lim[0], u0[0]));
  z1[0] = scale;
  u[0] = scale;
  absxk = scale * Qu[0];
  scale = fmaxf(lower_lim[1], fminf(upper_lim[1], u0[1]));
  z1[1] = scale;
  u[1] = scale;
  absxk += scale * Qu[1];
  scale = fmaxf(lower_lim[2], fminf(upper_lim[2], u0[2]));
  z1[2] = scale;
  u[2] = scale;
  absxk += scale * Qu[2];
  t = 0.0F;
  for (i5 = 0; i5 < 3; i5++) {
    t += ((0.5F * z1[0] * Quu[3 * i5] + 0.5F * z1[1] * Quu[1 + 3 * i5]) + 0.5F *
          scale * Quu[2 + 3 * i5]) * z1[i5];
  }

  cost = absxk + t;

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
        for (i = 0; i < 3; i++) {
          scale = Qu[i] + ((Quu[i] * u[0] + Quu[i + 3] * u[1]) + Quu[i + 6] * u
                           [2]);
          grad[i] = scale;
          prev_clamped[i] = clamped[i];
          out = false;
          b_clamped = false;
          if ((u[i] == lower_lim[i]) && (scale > 0.0F)) {
            out = true;
            b_clamped = true;
          }

          if ((u[i] == upper_lim[i]) && (scale < 0.0F)) {
            out = true;
            b_clamped = true;
          }

          b_free[i] = !out;
          clamped[i] = b_clamped;
        }

        /*  Check if all controls clamped */
        out = true;
        k = 0;
        exitg2 = false;
        while ((!exitg2) && (k < 3)) {
          if (!clamped[k]) {
            out = false;
            exitg2 = true;
          } else {
            k++;
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
            for (i5 = 0; i5 < trueCount; i5++) {
              for (i6 = 0; i6 < trueCount; i6++) {
                Quu_data[i6 + trueCount * i5] = Quu[(tmp_data[i6] + 3 *
                  (tmp_data[i5] - 1)) - 1];
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
              grad_norm = 0.0F;
            } else {
              y = 0.0F;
              if (trueCount == 1) {
                grad_norm = fabsf(grad[b_tmp_data[0] - 1]);
              } else {
                scale = 1.29246971E-26F;
                for (k = 0; k < trueCount; k++) {
                  absxk = fabsf(grad[b_tmp_data[k] - 1]);
                  if (absxk > scale) {
                    t = scale / absxk;
                    y = 1.0F + y * t * t;
                    scale = absxk;
                  } else {
                    t = absxk / scale;
                    y += t * t;
                  }
                }

                grad_norm = scale * sqrtf(y);
              }
            }

            if (grad_norm < 1.0E-8F) {
              *result = 5.0F;
              exitg1 = true;
            } else {
              /*  get search direction */
              scale = u[0] * (float)clamped[0];
              absxk = u[1] * (float)clamped[1];
              t = u[2] * (float)clamped[2];
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
              value = Luu_size[0];
              n = Luu_size[0];

              /*  Check sizes match and L is lower-triangular */
              /*  Solve LY = B for Y */
              /* ======================= */
              if (0 <= Luu_size[0] - 1) {
                memset(&Y_data[0], 0, (unsigned int)(Luu_size[0] * (int)sizeof
                        (float)));
              }

              if (0 <= Luu_size[0] - 1) {
                memset(&coeffs_data[0], 0, (unsigned int)(Luu_size[0] * (int)
                        sizeof(float)));
              }

              i5 = Luu_size[0];
              if (0 <= Luu_size[0] - 1) {
                a_size_idx_1 = Luu_size[0];
                loop_ub = Luu_size[0];
                b_loop_ub = Luu_size[0];
              }

              for (i = 0; i < i5; i++) {
                if (1 + i != 1) {
                  for (i6 = 0; i6 < i; i6++) {
                    coeffs_data[i6] = Luu_data[i + Luu_size[0] * i6];
                  }
                }

                if (0 <= loop_ub - 1) {
                  memcpy(&a_data[0], &coeffs_data[0], (unsigned int)(loop_ub *
                          (int)sizeof(float)));
                }

                if (0 <= b_loop_ub - 1) {
                  memcpy(&b_data[0], &Y_data[0], (unsigned int)(b_loop_ub * (int)
                          sizeof(float)));
                }

                if ((a_size_idx_1 == 1) || (Luu_size[0] == 1)) {
                  scale = 0.0F;
                  for (i6 = 0; i6 < a_size_idx_1; i6++) {
                    scale += a_data[i6] * b_data[i6];
                  }
                } else {
                  scale = 0.0F;
                  for (i6 = 0; i6 < a_size_idx_1; i6++) {
                    scale += a_data[i6] * b_data[i6];
                  }
                }

                Y_data[i] = (grad_clamped[c_tmp_data[i] - 1] - scale) /
                  Luu_data[i + Luu_size[0] * i];
              }

              /*  Solve L'X = Y for X */
              /* ======================= */
              if (0 <= Luu_size[0] - 1) {
                memset(&z1[0], 0, (unsigned int)(Luu_size[0] * (int)sizeof(float)));
              }

              LT_size_idx_0 = Luu_size[1];
              c_loop_ub = Luu_size[0];
              for (i5 = 0; i5 < c_loop_ub; i5++) {
                k = Luu_size[1];
                for (i6 = 0; i6 < k; i6++) {
                  LT_data[i6 + LT_size_idx_0 * i5] = Luu_data[i5 + Luu_size[0] *
                    i6];
                }
              }

              if (0 <= Luu_size[0] - 1) {
                memset(&coeffs_data[0], 0, (unsigned int)(Luu_size[0] * (int)
                        sizeof(float)));
              }

              i5 = (int)((1.0F + (-1.0F - (float)Luu_size[0])) / -1.0F);
              if (0 <= i5 - 1) {
                a_size_idx_1 = Luu_size[0];
                d_loop_ub = Luu_size[0];
                e_loop_ub = Luu_size[0];
              }

              for (i = 0; i < i5; i++) {
                b_i = (value - i) - 1;
                if (b_i + 1 != value) {
                  if (b_i + 2 > n) {
                    i6 = -1;
                    k = 0;
                    trueCount = -1;
                  } else {
                    i6 = b_i;
                    k = value;
                    trueCount = b_i;
                  }

                  c_loop_ub = k - i6;
                  for (k = 0; k <= c_loop_ub - 2; k++) {
                    coeffs_data[(trueCount + k) + 1] = LT_data[b_i +
                      LT_size_idx_0 * ((i6 + k) + 1)];
                  }
                }

                if (0 <= d_loop_ub - 1) {
                  memcpy(&a_data[0], &coeffs_data[0], (unsigned int)(d_loop_ub *
                          (int)sizeof(float)));
                }

                if (0 <= e_loop_ub - 1) {
                  memcpy(&b_data[0], &z1[0], (unsigned int)(e_loop_ub * (int)
                          sizeof(float)));
                }

                if ((a_size_idx_1 == 1) || (Luu_size[0] == 1)) {
                  scale = 0.0F;
                  for (i6 = 0; i6 < a_size_idx_1; i6++) {
                    scale += a_data[i6] * b_data[i6];
                  }
                } else {
                  scale = 0.0F;
                  for (i6 = 0; i6 < a_size_idx_1; i6++) {
                    scale += a_data[i6] * b_data[i6];
                  }
                }

                z1[b_i] = (Y_data[b_i] - scale) / LT_data[b_i + LT_size_idx_0 *
                  b_i];
              }

              c_loop_ub = Luu_size[0];
              for (i5 = 0; i5 < c_loop_ub; i5++) {
                b_data[i5] = -z1[i5] - u[c_tmp_data[i5] - 1];
              }

              k = 0;

              /*  cholesky solver */
              /*  check projected change in cost is a reduction */
              delta_u = 0.0F;
              if (b_free[0]) {
                delta_u = b_data[0];
                k = 1;
              }

              z1[0] = delta_u * grad[0];
              b_delta_u[0] = delta_u;
              delta_u = 0.0F;
              if (b_free[1]) {
                delta_u = b_data[k];
                k++;
              }

              z1[1] = delta_u * grad[1];
              b_delta_u[1] = delta_u;
              delta_u = 0.0F;
              if (b_free[2]) {
                delta_u = b_data[k];
              }

              y = (z1[0] + z1[1]) + delta_u * grad[2];
              if (y >= 0.0F) {
                /*  (should not happen) */
                exitg1 = true;
              } else {
                /*  Armijo linesearch */
                step = 1.0F;

                /*  Returns array x with all values clamped between lower and upper */
                scale = fmaxf(lower_lim[0], fminf(upper_lim[0], u[0] +
                  b_delta_u[0]));
                z1[0] = scale;
                u_c_idx_0 = scale;
                absxk = scale * Qu[0];
                scale = fmaxf(lower_lim[1], fminf(upper_lim[1], u[1] +
                  b_delta_u[1]));
                z1[1] = scale;
                u_c_idx_1 = scale;
                absxk += scale * Qu[1];
                scale = fmaxf(lower_lim[2], fminf(upper_lim[2], u[2] + delta_u));
                z1[2] = scale;
                u_c_idx_2 = scale;
                absxk += scale * Qu[2];
                t = 0.0F;
                for (i5 = 0; i5 < 3; i5++) {
                  t += ((0.5F * z1[0] * Quu[3 * i5] + 0.5F * z1[1] * Quu[1 + 3 *
                         i5]) + 0.5F * scale * Quu[2 + 3 * i5]) * z1[i5];
                }

                cost_c = absxk + t;
                exitg2 = false;
                while ((!exitg2) && ((cost_c - cost) / (step * y) < 0.1F)) {
                  step *= 0.6F;

                  /*  Returns array x with all values clamped between lower and upper */
                  scale = fmaxf(lower_lim[0], fminf(upper_lim[0], u[0] + step *
                    b_delta_u[0]));
                  z1[0] = scale;
                  u_c_idx_0 = scale;
                  absxk = scale * Qu[0];
                  scale = fmaxf(lower_lim[1], fminf(upper_lim[1], u[1] + step *
                    b_delta_u[1]));
                  z1[1] = scale;
                  u_c_idx_1 = scale;
                  absxk += scale * Qu[1];
                  scale = fmaxf(lower_lim[2], fminf(upper_lim[2], u[2] + step *
                    delta_u));
                  z1[2] = scale;
                  u_c_idx_2 = scale;
                  absxk += scale * Qu[2];
                  t = 0.0F;
                  for (i5 = 0; i5 < 3; i5++) {
                    t += ((0.5F * z1[0] * Quu[3 * i5] + 0.5F * z1[1] * Quu[1 + 3
                           * i5]) + 0.5F * scale * Quu[2 + 3 * i5]) * z1[i5];
                  }

                  cost_c = absxk + t;
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
  int nmj;
  int k;
  float b_A_data[9];
  int n;
  int info;
  int b_n;
  int j;
  boolean_T exitg1;
  int jj;
  float ajj;
  int ix;
  int iy;
  int i7;
  int i8;
  float c_A_data[9];
  int jp1j;
  int iac;
  float c;
  int ia;

  /*  Wrapper for MATLAB chol for use with auto coder */
  /*  Inputs: */
  /* =========== */
  /*  A - positive semi-definite matrix */
  A_size_idx_0 = A_size[0];
  nmj = A_size[1];
  k = A_size[0] * A_size[1];
  if (0 <= k - 1) {
    memcpy(&b_A_data[0], &A_data[0], (unsigned int)(k * (int)sizeof(float)));
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
        ajj = 0.0F;
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
        if (ajj > 0.0F) {
          ajj = sqrtf(ajj);
          b_A_data[jj] = ajj;
          if (j + 1 < b_n) {
            nmj = (b_n - j) - 1;
            k = j + 2;
            jp1j = jj + 2;
            if ((nmj == 0) || (j == 0)) {
            } else {
              ix = j;
              i7 = (j + n * (j - 1)) + 2;
              for (iac = k; n < 0 ? iac >= i7 : iac <= i7; iac += n) {
                c = -b_A_data[ix];
                iy = jj + 1;
                i8 = (iac + nmj) - 1;
                for (ia = iac; ia <= i8; ia++) {
                  b_A_data[iy] += b_A_data[ia - 1] * c;
                  iy++;
                }

                ix += n;
              }
            }

            ajj = 1.0F / ajj;
            i7 = jj + nmj;
            for (k = jp1j; k <= i7 + 1; k++) {
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
        b_A_data[k + A_size_idx_0 * (j - 1)] = 0.0F;
      }
    }

    if (1 > nmj) {
      k = 0;
      nmj = 0;
    } else {
      k = nmj;
    }

    for (i7 = 0; i7 < nmj; i7++) {
      for (i8 = 0; i8 < k; i8++) {
        c_A_data[i8 + k * i7] = b_A_data[i8 + A_size_idx_0 * i7];
      }
    }

    A_size_idx_0 = k;
    k *= nmj;
    if (0 <= k - 1) {
      memcpy(&b_A_data[0], &c_A_data[0], (unsigned int)(k * (int)sizeof(float)));
    }
  }

  *fail = (float)info;
  L_size[0] = A_size_idx_0;
  L_size[1] = nmj;
  k = A_size_idx_0 * nmj;
  if (0 <= k - 1) {
    memcpy(&L_data[0], &b_A_data[0], (unsigned int)(k * (int)sizeof(float)));
  }
}

/*
 * Arguments    : const float x[700]
 *                const float xg[7]
 *                const float u[297]
 *                const float K[1782]
 *                const float u_lims[6]
 *                float dt
 *                const float B_ECI[300]
 *                float xnew[700]
 *                float unew[297]
 *                float *cost
 * Return Type  : void
 */
static void forwardRollout(const float x[700], const float xg[7], const float u
  [297], const float K[1782], const float u_lims[6], float dt, const float
  B_ECI[300], float xnew[700], float unew[297], float *cost)
{
  int unew_tmp;
  float a;
  int k;
  int dx_tmp;
  float dx[6];
  int b_dx_tmp;
  int c_dx_tmp;
  float q_idx_0;
  int q_idx_1_tmp;
  int q_idx_2_tmp;
  int q_idx_3_tmp;
  float z1[3];
  float q[16];
  float f1;
  float q_error[4];
  int b_k;
  float b[7];
  float dxdot1[70];
  float b_xnew[7];
  float dxdot2[70];
  float x1[7];
  float f2;

  /*  get rid of everything after unew except cost */
  /*  Uses an rk method to roll out a trajectory */
  /*  Returns the new trajectory, cost and the derivatives along the trajectory */
  /*  If feed-forward and feed-back controls l, K are non-zero */
  /*  then the unew returned is the new control sequeunce. Otherwise unew = u */
  /*  Sizes */
  /*  Initialize outputs (Don't overwrite x or u) */
  memset(&xnew[0], 0, 700U * sizeof(float));
  memset(&unew[0], 0, 297U * sizeof(float));
  *cost = 0.0F;
  for (unew_tmp = 0; unew_tmp < 7; unew_tmp++) {
    xnew[unew_tmp] = x[unew_tmp];
  }

  a = 0.5F * dt;
  for (k = 0; k < 99; k++) {
    /*  Find the state error vector dx */
    dx_tmp = 7 * k + 4;
    dx[3] = xnew[dx_tmp] - x[dx_tmp];
    b_dx_tmp = 7 * k + 5;
    dx[4] = xnew[b_dx_tmp] - x[b_dx_tmp];
    c_dx_tmp = 7 * k + 6;
    dx[5] = xnew[c_dx_tmp] - x[c_dx_tmp];

    /*  Calculate error between qk and q_nom */
    /*  Defined as conj(q_nom)*qnew */
    /*  Returns error as Rodrigues parameters (3x1) */
    q_idx_0 = x[7 * k];
    q_idx_1_tmp = 1 + 7 * k;
    q_idx_2_tmp = 2 + 7 * k;
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
    q[4] = x[q_idx_1_tmp];
    q[8] = x[q_idx_2_tmp];
    q[12] = x[q_idx_3_tmp];
    q[1] = -x[q_idx_1_tmp];
    q[5] = q_idx_0;
    q[9] = x[q_idx_3_tmp];
    q[13] = -x[q_idx_2_tmp];
    q[2] = -x[q_idx_2_tmp];
    q[6] = -x[q_idx_3_tmp];
    q[10] = q_idx_0;
    q[14] = x[q_idx_1_tmp];
    q[3] = -x[q_idx_3_tmp];
    q[7] = x[q_idx_2_tmp];
    q[11] = -x[q_idx_1_tmp];
    q[15] = q_idx_0;
    q_idx_0 = 0.0F;
    for (unew_tmp = 0; unew_tmp < 4; unew_tmp++) {
      q_error[unew_tmp] = 0.0F;
      f1 = ((q[unew_tmp] * xnew[7 * k] + q[unew_tmp + 4] * xnew[q_idx_1_tmp]) +
            q[unew_tmp + 8] * xnew[q_idx_2_tmp]) + q[unew_tmp + 12] *
        xnew[q_idx_3_tmp];
      q_error[unew_tmp] = f1;
      q_idx_0 += f1 * f1;
    }

    q_idx_0 = sqrtf(q_idx_0);
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
      f1 = 0.0F;
      for (unew_tmp = 0; unew_tmp < 6; unew_tmp++) {
        f1 += K[(b_k + 3 * unew_tmp) + 18 * k] * dx[unew_tmp];
      }

      unew_tmp = b_k + 3 * k;
      unew[unew_tmp] = fminf(u_lims[3 + b_k], fmaxf(u_lims[b_k], u[unew_tmp] -
        f1));
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
    satellite_dynamics(*(float (*)[7])&xnew[7 * k], *(float (*)[3])&unew[3 * k],
                       *(float (*)[3])&B_ECI[3 * k], b, dxdot1);
    for (unew_tmp = 0; unew_tmp < 7; unew_tmp++) {
      b_xnew[unew_tmp] = xnew[unew_tmp + 7 * k] + a * b[unew_tmp];
    }

    satellite_dynamics(b_xnew, *(float (*)[3])&unew[3 * k], *(float (*)[3])&
                       B_ECI[3 * k], b, dxdot2);
    for (unew_tmp = 0; unew_tmp < 7; unew_tmp++) {
      x1[unew_tmp] = xnew[unew_tmp + 7 * k] + dt * b[unew_tmp];
    }

    /*  Re-normalize the quaternion */
    q_idx_0 = sqrtf(((x1[0] * x1[0] + x1[1] * x1[1]) + x1[2] * x1[2]) + x1[3] *
                    x1[3]);
    unew_tmp = 7 * (k + 1);
    xnew[unew_tmp] = x1[0] / q_idx_0;
    xnew[1 + unew_tmp] = x1[1] / q_idx_0;
    xnew[2 + unew_tmp] = x1[2] / q_idx_0;
    xnew[3 + unew_tmp] = x1[3] / q_idx_0;
    xnew[unew_tmp + 4] = x1[4];
    xnew[unew_tmp + 5] = x1[5];
    xnew[unew_tmp + 6] = x1[6];

    /*  Continuous time Jacobians */
    /*  Discrete time Jacobians */
    /*  First form the state error Jacobians E(x_k) & E(x_k+1) */
    /*  Calculate the cost */
    /*  make second function that's satellite cost derivatives */
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
    /*  cumulative angular velocity cost hessian */
    /*  cumulative geodesic cost weighting */
    /*  Finds the geodesic quaternion-error cost */
    /*  quat_cost = min(1+qd'q, 1-qd'q)) where qd is the desired attitude quaternion  */
    /*  Also records the sign (+ or -) which minimizes quat_cost */
    /*  this is used when calculating the Jacobain and Hessian */
    q_idx_0 = ((xg[0] * xnew[7 * k] + xg[1] * xnew[q_idx_1_tmp]) + xg[2] *
               xnew[q_idx_2_tmp]) + xg[3] * xnew[q_idx_3_tmp];
    z1[0] = xnew[dx_tmp] - xg[4];
    z1[1] = xnew[b_dx_tmp] - xg[5];
    z1[2] = xnew[c_dx_tmp] - xg[6];

    /*  State cost Hessian */
    /*  State cost Jacobian */
    /*  Control cost Hessian & Jacobian */
    f1 = 0.0F;
    for (unew_tmp = 0; unew_tmp < 3; unew_tmp++) {
      f1 += ((0.5F * z1[0] * fv0[3 * unew_tmp] + 0.5F * z1[1] * fv0[1 + 3 *
              unew_tmp]) + 0.5F * z1[2] * fv0[2 + 3 * unew_tmp]) * z1[unew_tmp];
    }

    f2 = 0.0F;
    for (unew_tmp = 0; unew_tmp < 3; unew_tmp++) {
      f2 += ((0.5F * unew[3 * k] * fv1[3 * unew_tmp] + 0.5F * unew[1 + 3 * k] *
              fv1[1 + 3 * unew_tmp]) + 0.5F * unew[2 + 3 * k] * fv1[2 + 3 *
             unew_tmp]) * unew[unew_tmp + 3 * k];
    }

    *cost += (fminf(1.0F + q_idx_0, 1.0F - q_idx_0) + f1) + f2;
  }

  /*  Final cost */
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
  q_idx_0 = ((xg[0] * xnew[693] + xg[1] * xnew[694]) + xg[2] * xnew[695]) + xg[3]
    * xnew[696];
  z1[0] = xnew[697] - xg[4];
  z1[1] = xnew[698] - xg[5];
  z1[2] = xnew[699] - xg[6];

  /*  State cost Hessian */
  /*  State cost Jacobian */
  /*  Control cost Hessian & Jacobian */
  f1 = 0.0F;
  for (unew_tmp = 0; unew_tmp < 3; unew_tmp++) {
    f1 += ((0.5F * z1[0] * (float)iv0[3 * unew_tmp] + 0.5F * z1[1] * (float)iv0
            [1 + 3 * unew_tmp]) + 0.5F * z1[2] * (float)iv0[2 + 3 * unew_tmp]) *
      z1[unew_tmp];
  }

  *cost += 10.0F * fminf(1.0F + q_idx_0, 1.0F - q_idx_0) + f1;
}

/*
 * Arguments    : const float x[99]
 * Return Type  : float
 */
static float mean(const float x[99])
{
  float y;
  int k;
  y = x[0];
  for (k = 0; k < 98; k++) {
    y += x[k + 1];
  }

  y /= 99.0F;
  return y;
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
 *                float xdot[7]
 *                float dxdot[70]
 * Return Type  : void
 */
static void satellite_dynamics(const float x[7], const float u[3], const float
  B_ECI[3], float xdot[7], float dxdot[70])
{
  float c;
  float q_out[4];
  float B_B_idx_0;
  float B_B_idx_1;
  float B_B_idx_2;
  float qdot_tmp[9];
  int i0;
  float wdot_tmp[9];
  float B_B[3];
  int i1;
  float fv4[12];
  float q1[3];
  int a_tmp;
  static const signed char a[9] = { 100, 0, 0, 0, 100, 0, 0, 0, 100 };

  int wdot_tmp_tmp;
  float b_wdot_tmp[9];
  static const float b[9] = { 0.01F, 0.0F, 0.0F, 0.0F, 0.01F, 0.0F, 0.0F, 0.0F,
    0.01F };

  float fv5[9];
  float b_a[9];
  static const short c_a[9] = { -200, 0, 0, 0, -200, 0, 0, 0, -200 };

  int dxdot_tmp;
  int b_dxdot_tmp;
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
  /*  this function multiplies quaternions */
  c = (B_ECI[0] * x[1] + B_ECI[1] * x[2]) + B_ECI[2] * x[3];
  q_out[1] = x[0] * B_ECI[0] + (B_ECI[1] * x[3] - B_ECI[2] * x[2]);
  q_out[2] = x[0] * B_ECI[1] + (B_ECI[2] * x[1] - B_ECI[0] * x[3]);
  q_out[3] = x[0] * B_ECI[2] + (B_ECI[0] * x[2] - B_ECI[1] * x[1]);

  /*  this function multiplies quaternions */
  B_B_idx_0 = (x[0] * q_out[1] + (0.0F - c) * -x[1]) + (-x[2] * q_out[3] - -x[3]
    * q_out[2]);
  B_B_idx_1 = (x[0] * q_out[2] + (0.0F - c) * -x[2]) + (-x[3] * q_out[1] - -x[1]
    * q_out[3]);
  B_B_idx_2 = (x[0] * q_out[3] + (0.0F - c) * -x[3]) + (-x[1] * q_out[2] - -x[2]
    * q_out[1]);

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
  qdot_tmp[0] = 0.0F;
  qdot_tmp[3] = -x[3];
  qdot_tmp[6] = x[2];
  qdot_tmp[1] = x[3];
  qdot_tmp[4] = 0.0F;
  qdot_tmp[7] = -x[1];
  qdot_tmp[2] = -x[2];
  qdot_tmp[5] = x[1];
  qdot_tmp[8] = 0.0F;
  for (i0 = 0; i0 < 9; i0++) {
    qdot_tmp[i0] += x[0] * (float)iv0[i0];
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
  wdot_tmp[0] = 0.0F;
  wdot_tmp[3] = -x[6];
  wdot_tmp[6] = x[5];
  wdot_tmp[1] = x[6];
  wdot_tmp[4] = 0.0F;
  wdot_tmp[7] = -x[4];
  wdot_tmp[2] = -x[5];
  wdot_tmp[5] = x[4];
  wdot_tmp[8] = 0.0F;

  /*  with magnetorquers */
  B_B[0] = -(B_B_idx_1 * u[2] - B_B_idx_2 * u[1]);
  B_B[1] = -(B_B_idx_2 * u[0] - B_B_idx_0 * u[2]);
  B_B[2] = -(B_B_idx_0 * u[1] - B_B_idx_1 * u[0]);
  for (i0 = 0; i0 < 3; i0++) {
    i1 = i0 << 2;
    fv4[i1] = 0.5F * -x[1 + i0];
    q1[i0] = 0.0F;
    for (a_tmp = 0; a_tmp < 3; a_tmp++) {
      wdot_tmp_tmp = i0 + 3 * a_tmp;
      b_wdot_tmp[wdot_tmp_tmp] = 0.0F;
      c = (wdot_tmp[i0] * b[3 * a_tmp] + wdot_tmp[i0 + 3] * b[1 + 3 * a_tmp]) +
        wdot_tmp[i0 + 6] * b[2 + 3 * a_tmp];
      b_wdot_tmp[wdot_tmp_tmp] = c;
      fv4[(a_tmp + i1) + 1] = 0.5F * qdot_tmp[a_tmp + 3 * i0];
      q1[i0] += c * x[4 + a_tmp];
    }

    B_B[i0] -= q1[i0];
  }

  for (i0 = 0; i0 < 4; i0++) {
    q_out[i0] = 0.0F;
    q_out[i0] = (fv4[i0] * x[4] + fv4[i0 + 4] * x[5]) + fv4[i0 + 8] * x[6];
  }

  for (i0 = 0; i0 < 3; i0++) {
    q1[i0] = 0.0F;
    q1[i0] = ((float)a[i0] * B_B[0] + (float)a[i0 + 3] * B_B[1]) + (float)a[i0 +
      6] * B_B[2];
  }

  xdot[0] = q_out[0];
  xdot[1] = q_out[1];
  xdot[2] = q_out[2];
  xdot[3] = q_out[3];

  /*  Jacobians  */
  for (i0 = 0; i0 < 3; i0++) {
    xdot[i0 + 4] = q1[i0];
    q1[i0] = (b[i0] * x[4] + b[i0 + 3] * x[5]) + b[i0 + 6] * x[6];
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
  fv5[0] = 0.0F;
  fv5[3] = -q1[2];
  fv5[6] = q1[1];
  fv5[1] = q1[2];
  fv5[4] = 0.0F;
  fv5[7] = -q1[0];
  fv5[2] = -q1[1];
  fv5[5] = q1[0];
  fv5[8] = 0.0F;
  for (i0 = 0; i0 < 9; i0++) {
    b_wdot_tmp[i0] -= fv5[i0];
  }

  for (i0 = 0; i0 < 3; i0++) {
    for (i1 = 0; i1 < 3; i1++) {
      a_tmp = i0 + 3 * i1;
      b_a[a_tmp] = 0.0F;
      b_a[a_tmp] = ((float)c_a[i0] * b_wdot_tmp[3 * i1] + (float)c_a[i0 + 3] *
                    b_wdot_tmp[1 + 3 * i1]) + (float)c_a[i0 + 6] * b_wdot_tmp[2
        + 3 * i1];
    }
  }

  fv5[0] = 0.0F;
  fv5[3] = -B_B_idx_2;
  fv5[6] = B_B_idx_1;
  fv5[1] = B_B_idx_2;
  fv5[4] = 0.0F;
  fv5[7] = -B_B_idx_0;
  fv5[2] = -B_B_idx_1;
  fv5[5] = B_B_idx_0;
  fv5[8] = 0.0F;
  dxdot[0] = 0.0F;
  for (i0 = 0; i0 < 3; i0++) {
    c = x[4 + i0];
    dxdot_tmp = 7 * (i0 + 1);
    dxdot[dxdot_tmp] = 0.5F * -c;
    b_dxdot_tmp = 7 * (i0 + 4);
    dxdot[b_dxdot_tmp] = 0.5F * -x[1 + i0];
    dxdot[i0 + 1] = 0.5F * c;
    for (i1 = 0; i1 < 3; i1++) {
      wdot_tmp_tmp = i0 + 3 * i1;
      b_wdot_tmp[wdot_tmp_tmp] = 0.0F;
      b_wdot_tmp[wdot_tmp_tmp] = ((float)d_a[i0] * fv5[3 * i1] + (float)d_a[i0 +
        3] * fv5[1 + 3 * i1]) + (float)d_a[i0 + 6] * fv5[2 + 3 * i1];
      a_tmp = i1 + 3 * i0;
      dxdot[(i1 + dxdot_tmp) + 1] = 0.5F * -wdot_tmp[a_tmp];
      dxdot[(i1 + b_dxdot_tmp) + 1] = 0.5F * qdot_tmp[a_tmp];
    }
  }

  for (i0 = 0; i0 < 4; i0++) {
    dxdot[7 * i0 + 4] = 0.0F;
    dxdot[7 * i0 + 5] = 0.0F;
    dxdot[7 * i0 + 6] = 0.0F;
  }

  for (i0 = 0; i0 < 3; i0++) {
    dxdot_tmp = 7 * (i0 + 4);
    dxdot[dxdot_tmp + 4] = 0.5F * b_a[3 * i0];
    b_dxdot_tmp = 1 + 3 * i0;
    dxdot[dxdot_tmp + 5] = 0.5F * b_a[b_dxdot_tmp];
    a_tmp = 2 + 3 * i0;
    dxdot[dxdot_tmp + 6] = 0.5F * b_a[a_tmp];
    dxdot_tmp = 7 * (i0 + 7);
    dxdot[dxdot_tmp] = 0.0F;
    dxdot[1 + dxdot_tmp] = 0.0F;
    dxdot[2 + dxdot_tmp] = 0.0F;
    dxdot[3 + dxdot_tmp] = 0.0F;
    dxdot[dxdot_tmp + 4] = b_wdot_tmp[3 * i0];
    dxdot[dxdot_tmp + 5] = b_wdot_tmp[b_dxdot_tmp];
    dxdot[dxdot_tmp + 6] = b_wdot_tmp[a_tmp];
  }
}

/*
 * Arguments    : float *lambda
 *                float direction
 * Return Type  : void
 */
static void updateLambda(float *lambda, float direction)
{
  /*  Increases or decreases the regularization parameter according */
  /*  to a non-linear scaling regime. */
  if (direction == 1.0F) {
    /*  increase lambda */
    *lambda = fmaxf(*lambda * 1.6F, 1.0E-6F);
  } else {
    /*  decrease lambda */
    *lambda = *lambda * 0.625F * (float)(*lambda > 1.0E-6F);

    /*  set = 0 if lambda too small */
  }

/*
 * Arguments    : const float x0[700]
 *                const float xg[7]
 *                const float u0[297]
 *                const float u_lims[6]
 *                float dt
 *                const float B_ECI[300]
 *                float x[700]
 *                float u[297]
 *                float K[1782]
 *                boolean_T *result
 * Return Type  : void
 */
void milqr_efficient(const float x0[700], const float xg[7], const float u0[297],
                     const float u_lims[6], float dt, const float B_ECI[300],
                     float x[700], float u[297], float K[1782], boolean_T
                     *result)
{
  float alphas[11];
  float lambda;
  float cost_n;
  static float x_n[700];
  int j;
  float u_n[297];
  float dV[2];
  float l[297];
  static float fv2[1782];
  float cost;
  int iter;
  int b_iter;
  boolean_T exitg1;
  boolean_T backPassCheck;
  int exitg2;
  boolean_T diverged;
  float z[297];
  float fv3[297];
  float maxval[99];
  float f0;
  boolean_T lineSearchCheck;
  boolean_T exitg3;
  float expected_change;
  float dcost;
  float c_ratio;

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
  lambda = 1.0F;

  /*  error state size (3 param. error representation for attitude) */
  /*  variable initialization */
  cost_n = 0.0F;
  memcpy(&x_n[0], &x0[0], 700U * sizeof(float));

  /*  Initial Forward rollout */
  for (j = 0; j < 297; j++) {
    u_n[j] = u0[j];
    l[j] = 0.0F;
  }

  memset(&K[0], 0, 1782U * sizeof(float));
  dV[0] = 0.0F;
  dV[1] = 0.0F;

  /*  Throw out everything after x and u, then in backward pass evalute on the */
  /*  spot */
  memset(&fv2[0], 0, 1782U * sizeof(float));
  forwardRollout(x0, xg, u0, fv2, u_lims, dt, B_ECI, x, u, &cost);

  /*  Convergence check params */
  /*  Expected cost change */
  /*  Ratio of cost change to expected cost change */
  *result = false;

  /*  fprintf("\n=====Running MILQR Optimisation====\n"); */
  iter = 1;
  b_iter = 0;
  exitg1 = false;
  while ((!exitg1) && (b_iter < 300)) {
    iter = b_iter + 1;

    /*      fprintf("\n---New Iteration---"); */
    /*      fprintf("\n lambda: "); */
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
        if (diverged) {
          /*  fprintf("---Warning: Cholesky factorizaton failed---\n"); */
          /*  Increase regularization parameter (lambda) */
          updateLambda(&lambda, 1.0F);
          if (lambda > 1.0E+7F) {
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
    b_abs(l, z);
    b_abs(u, fv3);
    for (j = 0; j < 297; j++) {
      z[j] /= fv3[j] + 1.0F;
    }

    for (j = 0; j < 99; j++) {
      maxval[j] = z[3 * j];
      f0 = z[3 * j + 1];
      if (maxval[j] < f0) {
        maxval[j] = f0;
      }

      f0 = z[3 * j + 2];
      if (maxval[j] < f0) {
        maxval[j] = f0;
      }
    }

    /*  Avg over time of max change */
    if ((lambda < 1.0E-5F) && (mean(maxval) < 0.0001F)) {
      /*  fprintf("\n---Success: Control change decreased below tolerance---\n"); */
      *result = true;
      exitg1 = true;
    } else {
      /*  Forward Line-Search */
      /* =========================================== */
      lineSearchCheck = false;
      if (backPassCheck) {
        j = 0;
        exitg3 = false;
        while ((!exitg3) && (j < 11)) {
          b_forwardRollout(x, xg, u, l, K, alphas[j], u_lims, dt, B_ECI, x_n,
                           u_n, &cost_n);
          expected_change = alphas[j] * dV[0] + powf(alphas[j], 2.0F) * dV[1];
          if (expected_change < 0.0F) {
            c_ratio = (cost_n - cost) / expected_change;
          } else {
            /*  Non-positive expected cost reduction */
            /*  actual cost change must be negative to accept the step */
            f0 = cost_n - cost;
            b_sign(&f0);
            c_ratio = -f0;
          }

          if (c_ratio > 0.0F) {
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
        updateLambda(&lambda, -1.0F);
        dcost = cost - cost_n;

        /*  Update the trajectory and controls */
        memcpy(&x[0], &x_n[0], 700U * sizeof(float));
        memcpy(&u[0], &u_n[0], 297U * sizeof(float));
        cost = cost_n;

        /*  Change in cost small enough to terminate? */
        if (dcost < 1.0E-7F) {
          *result = true;

          /*  fprintf('\n---Success: cost change < tolerance---\n'); */
          exitg1 = true;
        } else {
          b_iter++;
        }
      } else {
        /*  No cost reduction (based on cost change ratio) */
        /*  Increase lambda */
        updateLambda(&lambda, 1.0F);
        if (lambda > 1.0E+7F) {
          /*  fprintf("\n---Diverged: new lambda > lambda_max---\n"); */
          exitg1 = true;
        } else {
          b_iter++;
        }
      }
    }
  }

  if (iter == 300) {
    /*  Ddin't converge completely */
    *result = false;

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
